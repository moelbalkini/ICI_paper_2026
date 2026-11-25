# =============================================================================
# Complete analysis script: clonal beta-binomial mixed model + differential expression
# Run in the R session where gex_CD4_post, gex_CD8_post (re-clustered
# post-treatment T cells with TCR clonotypes) and GEX_mono_dc (monocytes/DCs)
# are loaded.
#
# Part 0  prepare_topclone_counts()
# Part 1  original clonal model, diagnostics and plots (unchanged)
# Part 2  FDR across all 14 clusters
# Part 3  cell-level differential expression in monocytes (FindMarkers: Wilcoxon, MAST; Fig. 3E)
# =============================================================================
library(dplyr)
library(ggplot2)

# ---------------------------------------------------------------- Part 0 -----
# -----------------------------------------------------------------------------
# 01_prepare_topclone_counts.R
# Build the per-clone x per-cluster count table used by the beta-binomial model.
#
# For every patient (SampleID), the `top_n` most expanded TCR clonotypes are
# selected. For each of these clones and each transcriptional cluster we record:
#   n                 cells of the clone
#   k                 cells of the clone inside the cluster
#   outside_cluster   n - k
#   N                 all TCR+ cells of the patient
#   K                 TCR+ cells of the patient inside the cluster
#   background_logit  logit(K / N), the patient-specific cluster frequency,
#                     used as offset so that the model tests over-representation
#                     of expanded clones relative to the patient's own background.
#
# Input: a Seurat object whose meta.data contains the columns
#   SampleID, Group ("Control" / "n-irAE"), clonotype, seurat_clusters
# -----------------------------------------------------------------------------

prepare_topclone_counts <- function(seu_obj, lineage, top_n = 10) {
  md <- seu_obj@meta.data |>
    dplyr::filter(!is.na(clonotype), clonotype != "") |>
    dplyr::mutate(
      Group = factor(Group, levels = c("Control", "n-irAE")),
      cluster_id = as.character(seurat_clusters)
    )

  cluster_order <- sort(unique(md$cluster_id))
  num <- suppressWarnings(as.numeric(cluster_order))
  if (!anyNA(num)) cluster_order <- cluster_order[order(num)]

  # patient background: all TCR+ cells per patient and cluster
  background <- md |>
    dplyr::count(SampleID, cluster_id, name = "K") |>
    tidyr::complete(SampleID, cluster_id = cluster_order,
                    fill = list(K = 0)) |>
    dplyr::group_by(SampleID) |>
    dplyr::mutate(N = sum(K)) |>
    dplyr::ungroup()

  # top-n expanded clonotypes per patient
  top_clones <- md |>
    dplyr::count(SampleID, Group, clonotype, name = "n") |>
    dplyr::arrange(SampleID, dplyr::desc(n)) |>
    dplyr::group_by(SampleID) |>
    dplyr::slice_head(n = top_n) |>
    dplyr::ungroup()

  # cells of each top clone per cluster (zeros included)
  clone_cluster <- md |>
    dplyr::semi_join(top_clones, by = c("SampleID", "clonotype")) |>
    dplyr::count(SampleID, clonotype, cluster_id, name = "k") |>
    tidyr::complete(tidyr::nesting(SampleID, clonotype),
                    cluster_id = cluster_order, fill = list(k = 0)) |>
    dplyr::left_join(top_clones, by = c("SampleID", "clonotype")) |>
    dplyr::left_join(background, by = c("SampleID", "cluster_id")) |>
    dplyr::mutate(
      lineage = lineage,
      panel = paste0("C", cluster_id),
      outside_cluster = n - k,
      background_logit = stats::qlogis(
        pmin(pmax(K / N, 1e-6), 1 - 1e-6)
      )
    )

  list(clone_cluster = clone_cluster, cluster_order = cluster_order)
}


# ---------------------------------------------------------------- Part 1 -----
cd4_prepared <- prepare_topclone_counts(
  seu_obj = gex_CD4_post,
  lineage = "CD4",
  top_n = 10
)

cd8_prepared <- prepare_topclone_counts(
  seu_obj = gex_CD8_post,
  lineage = "CD8",
  top_n = 10
)
##
cd4_model_df <- cd4_prepared$clone_cluster
cd8_model_df <- cd8_prepared$clone_cluster
##
cd4_model_df |>
  dplyr::summarise(
    patients = dplyr::n_distinct(SampleID),
    clones = dplyr::n_distinct(
      paste(SampleID, clonotype)
    ),
    clusters = dplyr::n_distinct(panel),
    rows = dplyr::n()
  )

cd8_model_df |>
  dplyr::summarise(
    patients = dplyr::n_distinct(SampleID),
    clones = dplyr::n_distinct(
      paste(SampleID, clonotype)
    ),
    clusters = dplyr::n_distinct(panel),
    rows = dplyr::n()
  )
table(cd4_model_df$panel, cd4_model_df$Group)
table(cd8_model_df$panel, cd8_model_df$Group)
fit_one_cluster_model <- function(
    model_df,
    cluster_name
) {
  
  dat <- model_df |>
    dplyr::filter(
      as.character(panel) == cluster_name
    ) |>
    droplevels()
  
  empty_result <- function(status_text) {
    data.frame(
      panel = cluster_name,
      n_observations = nrow(dat),
      n_patients = dplyr::n_distinct(dat$SampleID),
      n_control_clones = sum(dat$Group == "Control"),
      n_nirae_clones = sum(dat$Group == "n-irAE"),
      estimate_log_odds = NA_real_,
      SE = NA_real_,
      odds_ratio = NA_real_,
      CI_low = NA_real_,
      CI_high = NA_real_,
      p_value = NA_real_,
      convergence = FALSE,
      status = status_text,
      stringsAsFactors = FALSE
    )
  }
  
  if (nrow(dat) == 0) {
    return(
      list(
        model = NULL,
        data = dat,
        result = empty_result("No data")
      )
    )
  }
  
  if (dplyr::n_distinct(dat$Group) < 2) {
    return(
      list(
        model = NULL,
        data = dat,
        result = empty_result("Only one group available")
      )
    )
  }
  
  if (all(dat$k == 0)) {
    return(
      list(
        model = NULL,
        data = dat,
        result = empty_result(
          "All clone-cluster counts are zero"
        )
      )
    )
  }
  
  if (all(dat$k == dat$n)) {
    return(
      list(
        model = NULL,
        data = dat,
        result = empty_result(
          "All clone cells are inside the cluster"
        )
      )
    )
  }
  
  fit <- tryCatch(
    glmmTMB::glmmTMB(
      cbind(k, outside_cluster) ~
        Group +
        offset(background_logit) +
        (1 | SampleID),
      family = glmmTMB::betabinomial(
        link = "logit"
      ),
      data = dat
    ),
    error = function(e) e
  )
  
  if (inherits(fit, "error")) {
    return(
      list(
        model = NULL,
        data = dat,
        result = empty_result(
          paste("Model error:", fit$message)
        )
      )
    )
  }
  
  coefficient_table <- summary(fit)$coefficients$cond
  
  if (!"Groupn-irAE" %in% rownames(coefficient_table)) {
    return(
      list(
        model = fit,
        data = dat,
        result = empty_result(
          "Group coefficient not available"
        )
      )
    )
  }
  
  estimate <- coefficient_table[
    "Groupn-irAE",
    "Estimate"
  ]
  
  standard_error <- coefficient_table[
    "Groupn-irAE",
    "Std. Error"
  ]
  
  p_value <- coefficient_table[
    "Groupn-irAE",
    "Pr(>|z|)"
  ]
  
  converged <- isTRUE(fit$sdr$pdHess)
  
  result <- data.frame(
    panel = cluster_name,
    n_observations = nrow(dat),
    n_patients = dplyr::n_distinct(dat$SampleID),
    n_control_clones = sum(dat$Group == "Control"),
    n_nirae_clones = sum(dat$Group == "n-irAE"),
    estimate_log_odds = estimate,
    SE = standard_error,
    odds_ratio = exp(estimate),
    CI_low = exp(
      estimate - 1.96 * standard_error
    ),
    CI_high = exp(
      estimate + 1.96 * standard_error
    ),
    p_value = p_value,
    convergence = converged,
    status = ifelse(
      converged,
      "Model fitted",
      "Convergence warning"
    ),
    stringsAsFactors = FALSE
  )
  
  list(
    model = fit,
    data = dat,
    result = result
  )
}
fit_all_clusters <- function(
    model_df,
    cluster_order
) {
  
  panel_order <- paste0("C", cluster_order)
  
  model_list <- setNames(
    lapply(
      panel_order,
      function(cluster_name) {
        fit_one_cluster_model(
          model_df = model_df,
          cluster_name = cluster_name
        )
      }
    ),
    panel_order
  )
  
  statistics <- dplyr::bind_rows(
    lapply(
      model_list,
      function(x) x$result
    )
  ) |>
    dplyr::mutate(
      FDR = p.adjust(
        p_value,
        method = "BH"
      ),
      
      direction = dplyr::case_when(
        is.na(odds_ratio) ~ "Not testable",
        odds_ratio > 1 ~ "Higher in n-irAE",
        odds_ratio < 1 ~ "Higher in Control",
        TRUE ~ "No estimated difference"
      )
    ) |>
    dplyr::arrange(p_value)
  
  fitted_models <- lapply(
    model_list,
    function(x) x$model
  )
  
  model_data <- lapply(
    model_list,
    function(x) x$data
  )
  
  list(
    statistics = statistics,
    models = fitted_models,
    data = model_data
  )
}

cd4_mixed_results <- fit_all_clusters(
  model_df = cd4_model_df,
  cluster_order = cd4_prepared$cluster_order
)

cd8_mixed_results <- fit_all_clusters(
  model_df = cd8_model_df,
  cluster_order = cd8_prepared$cluster_order
)
cd4_mixed_stats <- cd4_mixed_results$statistics
cd8_mixed_stats <- cd8_mixed_results$statistics
cd4_mixed_stats
cd8_mixed_stats
cd4_mixed_stats |>
  dplyr::select(
    panel,
    n_observations,
    n_patients,
    odds_ratio,
    CI_low,
    CI_high,
    p_value,
    FDR,
    direction,
    convergence,
    status
  )

cd8_mixed_stats |>
  dplyr::select(
    panel,
    n_observations,
    n_patients,
    odds_ratio,
    CI_low,
    CI_high,
    p_value,
    FDR,
    direction,
    convergence,
    status
  )
sapply(
  cd4_mixed_results$models,
  function(model) {
    if (is.null(model)) return(NA)
    model$sdr$pdHess
  }
)

sapply(
  cd8_mixed_results$models,
  function(model) {
    if (is.null(model)) return(NA)
    model$sdr$pdHess
  }
)
###


compare_model_families <- function(
    model_df,
    cluster_name
) {
  
  dat <- model_df |>
    dplyr::filter(
      as.character(panel) == cluster_name
    ) |>
    droplevels()
  
  fit_binomial <- glmmTMB::glmmTMB(
    cbind(k, outside_cluster) ~
      Group +
      offset(background_logit) +
      (1 | SampleID),
    family = binomial(link = "logit"),
    data = dat
  )
  
  fit_betabinomial <- glmmTMB::glmmTMB(
    cbind(k, outside_cluster) ~
      Group +
      offset(background_logit) +
      (1 | SampleID),
    family = glmmTMB::betabinomial(
      link = "logit"
    ),
    data = dat
  )
  
  list(
    AIC = AIC(
      fit_binomial,
      fit_betabinomial
    ),
    binomial = fit_binomial,
    betabinomial = fit_betabinomial
  )
}
cd8_c2_family_check <- compare_model_families(
  cd8_model_df,
  "C2"
)

cd8_c2_family_check$AIC
run_model_diagnostics <- function(
    fitted_model,
    simulation_n = 1000,
    seed = 123
) {
  
  if (is.null(fitted_model)) {
    stop("The supplied model is NULL.")
  }
  
  set.seed(seed)
  
  simulated_residuals <- DHARMa::simulateResiduals(
    fittedModel = fitted_model,
    n = simulation_n
  )
  
  plot(simulated_residuals)
  
  list(
    residuals = simulated_residuals,
    uniformity = DHARMa::testUniformity(
      simulated_residuals
    ),
    dispersion = DHARMa::testDispersion(
      simulated_residuals
    ),
    zero_inflation = DHARMa::testZeroInflation(
      simulated_residuals
    )
  )
}
diag_cd8_c2 <- run_model_diagnostics(
  cd8_mixed_results$models[["C2"]]
)

diag_cd8_c2$uniformity
diag_cd8_c2$dispersion
diag_cd8_c2$zero_inflation

plot_mixed_model_forest <- function(
    stats_df,
    title_text
) {
  
  plot_df <- stats_df |>
    dplyr::filter(
      !is.na(odds_ratio),
      !is.na(CI_low),
      !is.na(CI_high),
      CI_low > 0,
      CI_high > 0
    ) |>
    dplyr::mutate(
      panel = factor(
        panel,
        levels = panel[
          order(odds_ratio)
        ]
      )
    )
  
  ggplot(
    plot_df,
    aes(
      x = odds_ratio,
      y = panel
    )
  ) +
    geom_vline(
      xintercept = 1,
      linetype = "dashed",
      linewidth = 0.45
    ) +
    geom_errorbar(
      aes(
        xmin = CI_low,
        xmax = CI_high
      ),
      orientation = "y",
      width = 0.15,
      linewidth = 0.6
    ) +
    geom_point(size = 3) +
    scale_x_log10() +
    theme_classic(base_size = 14) +
    labs(
      title = title_text,
      x = "Odds ratio: n-irAE versus Control\n(log scale)",
      y = NULL
    ) +
    theme(
      plot.title = element_text(
        face = "bold",
        hjust = 0.5
      ),
      axis.text.y = element_text(
        face = "bold"
      )
    )
}
p_cd4_forest <- plot_mixed_model_forest(
  cd4_mixed_stats,
  "CD4 top-clone cluster representation"
)

p_cd8_forest <- plot_mixed_model_forest(
  cd8_mixed_stats,
  "CD8 top-clone cluster representation"
)

p_cd4_forest
p_cd8_forest
make_patient_plot_data <- function(model_df) {
  
  model_df |>
    dplyr::group_by(
      lineage,
      SampleID,
      Group,
      cluster_id,
      panel
    ) |>
    dplyr::summarise(
      total_topclone_cells = sum(n),
      topclone_cluster_cells = sum(k),
      N = dplyr::first(N),
      K = dplyr::first(K),
      
      observed_fraction =
        topclone_cluster_cells /
        total_topclone_cells,
      
      expected_fraction = K / N,
      
      enrichment_ratio =
        observed_fraction /
        expected_fraction,
      
      .groups = "drop"
    ) |>
    dplyr::mutate(
      Group = factor(
        Group,
        levels = c("Control", "n-irAE")
      )
    )
}
cd4_patient_plot_df <- make_patient_plot_data(
  cd4_model_df
)

cd8_patient_plot_df <- make_patient_plot_data(
  cd8_model_df
)
plot_patient_ratios <- function(
    patient_df,
    stats_df,
    title_text
) {
  
  panel_order_df <- patient_df |>
    dplyr::distinct(
      panel,
      cluster_id
    ) |>
    dplyr::mutate(
      cluster_num = suppressWarnings(
        as.numeric(cluster_id)
      )
    ) |>
    dplyr::arrange(
      is.na(cluster_num),
      cluster_num,
      cluster_id
    )
  
  panel_order <- as.character(
    panel_order_df$panel
  )
  
  plot_df <- patient_df |>
    dplyr::mutate(
      panel = factor(
        panel,
        levels = panel_order
      )
    )
  
  panel_y <- plot_df |>
    dplyr::group_by(panel) |>
    dplyr::summarise(
      y_max = max(
        enrichment_ratio,
        na.rm = TRUE
      ),
      .groups = "drop"
    ) |>
    dplyr::mutate(
      y_max = ifelse(
        is.finite(y_max),
        y_max,
        1
      ),
      label_y = pmax(
        y_max * 1.18,
        1.3
      )
    )
  
  annotation_df <- stats_df |>
    dplyr::mutate(
      panel = factor(
        panel,
        levels = panel_order
      ),
      
      label = dplyr::case_when(
        is.na(p_value) ~ "not testable",
        
        p_value < 0.001 ~ paste0(
          "p = ",
          format(
            p_value,
            scientific = TRUE,
            digits = 2
          )
        ),
        
        TRUE ~ paste0(
          "p = ",
          signif(p_value, 3)
        )
      )
    ) |>
    dplyr::left_join(
      panel_y,
      by = "panel"
    )
  
  ggplot(
    plot_df,
    aes(
      x = Group,
      y = enrichment_ratio
    )
  ) +
    geom_boxplot(
      width = 0.48,
      outlier.shape = NA,
      fill = "white",
      color = "black",
      linewidth = 0.6
    ) +
    geom_jitter(
      aes(color = Group),
      width = 0.10,
      height = 0,
      size = 2.8,
      alpha = 0.9
    ) +
    geom_hline(
      yintercept = 1,
      linetype = "dashed",
      linewidth = 0.4
    ) +
    geom_text(
      data = annotation_df,
      aes(
        x = 1.5,
        y = label_y,
        label = label
      ),
      inherit.aes = FALSE,
      size = 3.5
    ) +
    facet_wrap(
      ~panel,
      scales = "free_y"
    ) +
    scale_color_manual(
      values = c(
        "Control" = "#8E7CC3",
        "n-irAE" = "#C55A11"
      )
    ) +
    scale_y_continuous(
      expand = expansion(
        mult = c(0.05, 0.25)
      )
    ) +
    theme_classic(base_size = 14) +
    labs(
      title = title_text,
      x = NULL,
      y = "Patient-level observed/expected\nclone representation"
    ) +
    theme(
      strip.background = element_blank(),
      strip.text = element_text(
        face = "bold"
      ),
      legend.position = "none",
      plot.title = element_text(
        face = "bold",
        hjust = 0.5
      )
    )
}
p_cd4_patient <- plot_patient_ratios(
  cd4_patient_plot_df,
  cd4_mixed_stats,
  "CD4 top-clone cluster representation"
)

p_cd8_patient <- plot_patient_ratios(
  cd8_patient_plot_df,
  cd8_mixed_stats,
  "CD8 top-clone cluster representation"
)

p_cd4_patient
p_cd8_patient


# ---------------------------------------------------------------- Part 2 -----
# BH correction across all CD4 (8) and CD8 (6) clusters together.
# The column "FDR" above is corrected within each lineage only.
all_cluster_stats <- dplyr::bind_rows(
  dplyr::mutate(cd4_mixed_stats, lineage = "CD4"),
  dplyr::mutate(cd8_mixed_stats, lineage = "CD8")
) |>
  dplyr::mutate(FDR_all_clusters = p.adjust(p_value, method = "BH"))

all_cluster_stats |>
  dplyr::select(lineage, panel, odds_ratio, CI_low, CI_high,
                p_value, FDR, FDR_all_clusters, convergence)
write.csv(all_cluster_stats, "betabinomial_cluster_statistics.csv", row.names = FALSE)

# CHECK: CD4 TEMRA (cluster 8) should give OR 3.13, p = 0.027 as in the paper.

# ---------------------------------------------------------------- Part 3 -----
# Cell-level differential expression in monocytes/DCs (Figure 3E)

library(Seurat)
library(ggplot2)
library(ggpubr)
library(ggrepel)

genes <- c("IL18", "CCL2", "IL18R1", "IL18RAP", "IL18BP", "IL1B", "TNF", "IL6",
           "CXCL8", "CCL3", "CCL4", "NLRP3", "CXCL2")

GEX_mono_dc_nirAE <- subset(GEX_mono_dc, subset = Dataset == "nirAE")
unique(GEX_mono_dc_nirAE$Group_Timepoint)

DotPlot(GEX_mono_dc_nirAE, features = genes, group.by = "Group_Timepoint") + RotatedAxis()

Idents(GEX_mono_dc_nirAE) <- "Group_Timepoint"

# ---- 1. targeted gene panel: Wilcoxon --------------------------------------
de_results <- FindMarkers(
  object = GEX_mono_dc_nirAE,
  ident.1 = "n-irAE_post",
  ident.2 = "Control_post",
  features = genes,
  test.use = "wilcox",
  min.pct = 0.1,
  logfc.threshold = 0
)
print(de_results[de_results$p_val < 0.05, ])

for (g in genes) {
  plot_data <- FetchData(GEX_mono_dc_nirAE, vars = c(g, "Group_Timepoint"))
  colnames(plot_data) <- c("Expression", "Group")
  plot_data <- plot_data[plot_data$Group %in% c("Control_post", "n-irAE_post"), ]
  plot_data$Group <- factor(plot_data$Group, levels = c("Control_post", "n-irAE_post"))

  p_val <- de_results[g, "p_val"]
  p_star <- ifelse(p_val < 0.001, "***",
                   ifelse(p_val < 0.01, "**",
                          ifelse(p_val < 0.05, "*", "ns")))

  p <- ggplot(plot_data, aes(x = Group, y = Expression)) +
    geom_violin(fill = "white", color = "black", alpha = 0.1) +
    geom_jitter(aes(shape = Group), width = 0.2, size = 0.8, alpha = 0.3) +
    stat_summary(fun.data = "mean_sdl", geom = "errorbar", width = 0.1, color = "black") +
    stat_summary(fun = "mean", geom = "crossbar", width = 0.3, color = "red") +
    annotate("text", x = 1.5, y = max(plot_data$Expression) * 1.1,
             label = p_star, size = 6, fontface = "bold") +
    theme_classic() +
    scale_shape_manual(values = c("Control_post" = 15, "n-irAE_post" = 17)) +
    labs(title = g, subtitle = paste("Wilcoxon p =", format.pval(p_val, digits = 2)), x = NULL) +
    theme(legend.position = "none")
  print(p)
}

# ---- 2. transcriptome-wide: MAST + volcano (Figure 3E) ----------------------
all_de_results <- FindMarkers(
  object = GEX_mono_dc_nirAE,
  ident.1 = "n-irAE_post",
  ident.2 = "Control_post",
  test.use = "MAST",
  min.pct = 0.1,
  logfc.threshold = 0
)
all_de_results$gene <- rownames(all_de_results)
write.csv(all_de_results, "DE_MAST_monocytes_nirAE_post_vs_Control_post.csv", row.names = FALSE)

all_de_results$diffexpressed <- "Not Significant"
all_de_results$diffexpressed[all_de_results$avg_log2FC >  0.25 & all_de_results$p_val_adj < 0.05] <- "Upregulated"
all_de_results$diffexpressed[all_de_results$avg_log2FC < -0.25 & all_de_results$p_val_adj < 0.05] <- "Downregulated"

ggplot(all_de_results, aes(x = avg_log2FC, y = -log10(p_val_adj))) +
  geom_point(aes(color = diffexpressed), alpha = 0.4, size = 1) +
  scale_color_manual(values = c("Downregulated" = "#377eb8",
                                "Not Significant" = "grey80",
                                "Upregulated" = "#e41a1c")) +
  geom_text_repel(data = subset(all_de_results, gene %in% genes),
                  aes(label = gene), box.padding = 0.5, point.padding = 0.5,
                  force = 5, segment.color = "grey50", fontface = "bold") +
  geom_vline(xintercept = c(-0.25, 0.25), linetype = "dashed", alpha = 0.5) +
  geom_hline(yintercept = -log10(0.05), linetype = "dashed", alpha = 0.5) +
  theme_classic() +
  labs(title = "MAST Differential Expression: n-irAE_post vs Control_post",
       x = "log2(Fold Change)", y = "-log10(Adj. p-value)")
