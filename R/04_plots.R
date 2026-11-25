# -----------------------------------------------------------------------------
# 04_plots.R
# Forest plot of odds ratios (Figure 3B) and patient-level observed/expected
# clone representation per cluster (Supplementary Figure 7).
# -----------------------------------------------------------------------------
library(ggplot2)

plot_mixed_model_forest <- function(stats_df, title_text) {
  plot_df <- stats_df |>
    dplyr::filter(!is.na(odds_ratio), !is.na(CI_low), !is.na(CI_high),
                  CI_low > 0, CI_high > 0) |>
    dplyr::mutate(panel = factor(panel, levels = panel[order(odds_ratio)]))
  ggplot(plot_df, aes(x = odds_ratio, y = panel)) +
    geom_vline(xintercept = 1, linetype = "dashed", linewidth = 0.45) +
    geom_errorbar(aes(xmin = CI_low, xmax = CI_high), orientation = "y",
                  width = 0.15, linewidth = 0.6) +
    geom_point(size = 3) +
    scale_x_log10() +
    theme_classic(base_size = 14) +
    labs(title = title_text,
         x = "Odds ratio: n-irAE versus Control\n(log scale)", y = NULL) +
    theme(plot.title = element_text(face = "bold", hjust = 0.5),
          axis.text.y = element_text(face = "bold"))
}

make_patient_plot_data <- function(model_df) {
  model_df |>
    dplyr::group_by(lineage, SampleID, Group, cluster_id, panel) |>
    dplyr::summarise(
      total_topclone_cells = sum(n),
      topclone_cluster_cells = sum(k),
      N = dplyr::first(N), K = dplyr::first(K),
      observed_fraction = topclone_cluster_cells / total_topclone_cells,
      expected_fraction = K / N,
      enrichment_ratio = observed_fraction / expected_fraction,
      .groups = "drop"
    ) |>
    dplyr::mutate(Group = factor(Group, levels = c("Control", "n-irAE")))
}

plot_patient_ratios <- function(patient_df, stats_df, title_text) {
  panel_order <- patient_df |>
    dplyr::distinct(panel, cluster_id) |>
    dplyr::mutate(cluster_num = suppressWarnings(as.numeric(cluster_id))) |>
    dplyr::arrange(is.na(cluster_num), cluster_num, cluster_id) |>
    dplyr::pull(panel) |>
    as.character()
  plot_df <- dplyr::mutate(patient_df, panel = factor(panel, levels = panel_order))
  panel_y <- plot_df |>
    dplyr::group_by(panel) |>
    dplyr::summarise(y_max = max(enrichment_ratio, na.rm = TRUE), .groups = "drop") |>
    dplyr::mutate(y_max = ifelse(is.finite(y_max), y_max, 1),
                  label_y = pmax(y_max * 1.18, 1.3))
  annotation_df <- stats_df |>
    dplyr::mutate(
      panel = factor(panel, levels = panel_order),
      label = dplyr::case_when(
        is.na(p_value) ~ "not testable",
        p_value < 0.001 ~ paste0("p = ", format(p_value, scientific = TRUE, digits = 2)),
        TRUE ~ paste0("p = ", signif(p_value, 3))
      )
    ) |>
    dplyr::left_join(panel_y, by = "panel")
  ggplot(plot_df, aes(x = Group, y = enrichment_ratio)) +
    geom_boxplot(width = 0.48, outlier.shape = NA, fill = "white",
                 color = "black", linewidth = 0.6) +
    geom_jitter(aes(color = Group), width = 0.10, height = 0, size = 2.8, alpha = 0.9) +
    geom_hline(yintercept = 1, linetype = "dashed", linewidth = 0.4) +
    geom_text(data = annotation_df, aes(x = 1.5, y = label_y, label = label),
              inherit.aes = FALSE, size = 3.5) +
    facet_wrap(~panel, scales = "free_y") +
    scale_color_manual(values = c("Control" = "#8E7CC3", "n-irAE" = "#C55A11")) +
    scale_y_continuous(expand = expansion(mult = c(0.05, 0.25))) +
    theme_classic(base_size = 14) +
    labs(title = title_text, x = NULL,
         y = "Patient-level observed/expected\nclone representation") +
    theme(strip.background = element_blank(),
          strip.text = element_text(face = "bold"),
          legend.position = "none",
          plot.title = element_text(face = "bold", hjust = 0.5))
}
