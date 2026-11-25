# -----------------------------------------------------------------------------
# 02_fit_betabinomial_models.R
# Per-cluster beta-binomial generalized linear mixed model (glmmTMB):
#   cbind(k, outside_cluster) ~ Group + offset(background_logit) + (1 | SampleID)
# Group effect = log odds ratio (n-irAE vs Control) of an expanded clone's cells
# falling into the cluster, relative to the patient's own cluster frequency.
# Wald 95% CI; convergence = positive-definite Hessian; BH-FDR across clusters.
# -----------------------------------------------------------------------------

fit_one_cluster_model <- function(model_df, cluster_name) {
  dat <- model_df |>
    dplyr::filter(as.character(panel) == cluster_name) |>
    droplevels()

  empty_result <- function(status_text) {
    data.frame(
      panel = cluster_name,
      n_observations = nrow(dat),
      n_patients = dplyr::n_distinct(dat$SampleID),
      n_control_clones = sum(dat$Group == "Control"),
      n_nirae_clones = sum(dat$Group == "n-irAE"),
      estimate_log_odds = NA_real_, SE = NA_real_,
      odds_ratio = NA_real_, CI_low = NA_real_, CI_high = NA_real_,
      p_value = NA_real_, convergence = FALSE, status = status_text,
      stringsAsFactors = FALSE
    )
  }

  if (nrow(dat) == 0)
    return(list(model = NULL, data = dat, result = empty_result("No data")))
  if (dplyr::n_distinct(dat$Group) < 2)
    return(list(model = NULL, data = dat,
                result = empty_result("Only one group available")))
  if (all(dat$k == 0))
    return(list(model = NULL, data = dat,
                result = empty_result("All clone-cluster counts are zero")))
  if (all(dat$k == dat$n))
    return(list(model = NULL, data = dat,
                result = empty_result("All clone cells are inside the cluster")))

  fit <- tryCatch(
    glmmTMB::glmmTMB(
      cbind(k, outside_cluster) ~ Group + offset(background_logit) +
        (1 | SampleID),
      family = glmmTMB::betabinomial(link = "logit"),
      data = dat
    ),
    error = function(e) e
  )
  if (inherits(fit, "error"))
    return(list(model = NULL, data = dat,
                result = empty_result(paste("Model error:", fit$message))))

  coef_tab <- summary(fit)$coefficients$cond
  if (!"Groupn-irAE" %in% rownames(coef_tab))
    return(list(model = fit, data = dat,
                result = empty_result("Group coefficient not available")))

  est <- coef_tab["Groupn-irAE", "Estimate"]
  se  <- coef_tab["Groupn-irAE", "Std. Error"]
  p   <- coef_tab["Groupn-irAE", "Pr(>|z|)"]
  converged <- isTRUE(fit$sdr$pdHess)

  result <- data.frame(
    panel = cluster_name,
    n_observations = nrow(dat),
    n_patients = dplyr::n_distinct(dat$SampleID),
    n_control_clones = sum(dat$Group == "Control"),
    n_nirae_clones = sum(dat$Group == "n-irAE"),
    estimate_log_odds = est, SE = se,
    odds_ratio = exp(est),
    CI_low = exp(est - 1.96 * se), CI_high = exp(est + 1.96 * se),
    p_value = p, convergence = converged,
    status = ifelse(converged, "Model fitted", "Convergence warning"),
    stringsAsFactors = FALSE
  )
  list(model = fit, data = dat, result = result)
}

fit_all_clusters <- function(model_df, cluster_order) {
  panel_order <- paste0("C", cluster_order)
  model_list <- setNames(
    lapply(panel_order, function(cl) fit_one_cluster_model(model_df, cl)),
    panel_order
  )
  statistics <- dplyr::bind_rows(lapply(model_list, function(x) x$result)) |>
    dplyr::mutate(
      FDR_within_lineage = p.adjust(p_value, method = "BH"),
      direction = dplyr::case_when(
        is.na(odds_ratio) ~ "Not testable",
        odds_ratio > 1 ~ "Higher in n-irAE",
        odds_ratio < 1 ~ "Higher in Control",
        TRUE ~ "No estimated difference"
      )
    ) |>
    dplyr::arrange(p_value)
  list(
    statistics = statistics,
    models = lapply(model_list, function(x) x$model),
    data = lapply(model_list, function(x) x$data)
  )
}

# BH correction across all CD4 and CD8 clusters together (14 tests)
add_fdr_all_clusters <- function(cd4_stats, cd8_stats) {
  dplyr::bind_rows(
    dplyr::mutate(cd4_stats, lineage = "CD4"),
    dplyr::mutate(cd8_stats, lineage = "CD8")
  ) |>
    dplyr::mutate(FDR_all_clusters = p.adjust(p_value, method = "BH"))
}
