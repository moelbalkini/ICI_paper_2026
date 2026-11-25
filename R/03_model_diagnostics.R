# -----------------------------------------------------------------------------
# 03_model_diagnostics.R
# Binomial vs beta-binomial comparison (AIC) and DHARMa residual diagnostics.
# -----------------------------------------------------------------------------

compare_model_families <- function(model_df, cluster_name) {
  dat <- model_df |>
    dplyr::filter(as.character(panel) == cluster_name) |>
    droplevels()
  f <- cbind(k, outside_cluster) ~ Group + offset(background_logit) +
    (1 | SampleID)
  fit_bin  <- glmmTMB::glmmTMB(f, family = binomial(link = "logit"), data = dat)
  fit_bbin <- glmmTMB::glmmTMB(f, family = glmmTMB::betabinomial(link = "logit"),
                               data = dat)
  list(AIC = AIC(fit_bin, fit_bbin), binomial = fit_bin, betabinomial = fit_bbin)
}

run_model_diagnostics <- function(fitted_model, simulation_n = 1000, seed = 123) {
  if (is.null(fitted_model)) stop("The supplied model is NULL.")
  set.seed(seed)
  res <- DHARMa::simulateResiduals(fittedModel = fitted_model, n = simulation_n)
  plot(res)
  list(
    residuals = res,
    uniformity = DHARMa::testUniformity(res),
    dispersion = DHARMa::testDispersion(res),
    zero_inflation = DHARMa::testZeroInflation(res)
  )
}
