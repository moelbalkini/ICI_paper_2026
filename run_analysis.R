# -----------------------------------------------------------------------------
# run_analysis.R
# Patient-level beta-binomial mixed-effects model of clonal representation
# across CD4+ and CD8+ T-cell clusters (Giguelay et al. 2025, GSE278119,
# post-treatment samples; control n = 9, n-irAE n = 8).
#
# Requires two Seurat objects with re-clustered post-treatment T cells and TCR
# clonotypes in meta.data (columns: SampleID, Group, clonotype, seurat_clusters):
#   data/gex_CD4_post.rds, data/gex_CD8_post.rds
# These are generated from the public GEO data with the preprocessing described
# in the Methods of the manuscript (SoupX, Seurat, Harmony, Azimuth mapping).
# -----------------------------------------------------------------------------

source("R/01_prepare_topclone_counts.R")
source("R/02_fit_betabinomial_models.R")
source("R/03_model_diagnostics.R")
source("R/04_plots.R")

gex_CD4_post <- readRDS("data/gex_CD4_post.rds")
gex_CD8_post <- readRDS("data/gex_CD8_post.rds")

# 1. per-clone x per-cluster counts (top 10 clones per patient)
cd4_prepared <- prepare_topclone_counts(gex_CD4_post, lineage = "CD4", top_n = 10)
cd8_prepared <- prepare_topclone_counts(gex_CD8_post, lineage = "CD8", top_n = 10)
cd4_model_df <- cd4_prepared$clone_cluster
cd8_model_df <- cd8_prepared$clone_cluster

# 2. beta-binomial GLMM per cluster
cd4_res <- fit_all_clusters(cd4_model_df, cd4_prepared$cluster_order)
cd8_res <- fit_all_clusters(cd8_model_df, cd8_prepared$cluster_order)
all_stats <- add_fdr_all_clusters(cd4_res$statistics, cd8_res$statistics)

dir.create("results", showWarnings = FALSE)
write.csv(all_stats, "results/betabinomial_cluster_statistics.csv", row.names = FALSE)
print(all_stats[, c("lineage", "panel", "odds_ratio", "CI_low", "CI_high",
                    "p_value", "FDR_within_lineage", "FDR_all_clusters",
                    "convergence")])

# 3. diagnostics (example: CD8 cluster C2)
print(compare_model_families(cd8_model_df, "C2")$AIC)
diag_cd8_c2 <- run_model_diagnostics(cd8_res$models[["C2"]])

# 4. figures
ggsave("results/Fig3B_forest_CD4.pdf",
       plot_mixed_model_forest(cd4_res$statistics, "CD4 top-clone cluster representation"),
       width = 6, height = 5)
ggsave("results/Fig3B_forest_CD8.pdf",
       plot_mixed_model_forest(cd8_res$statistics, "CD8 top-clone cluster representation"),
       width = 6, height = 4)
ggsave("results/SuppFig7A_CD4_patient_ratios.pdf",
       plot_patient_ratios(make_patient_plot_data(cd4_model_df), cd4_res$statistics,
                           "CD4 top-clone cluster representation"),
       width = 10, height = 8)
ggsave("results/SuppFig7B_CD8_patient_ratios.pdf",
       plot_patient_ratios(make_patient_plot_data(cd8_model_df), cd8_res$statistics,
                           "CD8 top-clone cluster representation"),
       width = 10, height = 7)

# 5. cell-level differential expression in monocytes/DCs (Figure 3E):
#    requires data/gex_mono_dc.rds (monocytes/DCs with Dataset and
#    Group_Timepoint in meta.data)
if (file.exists("data/gex_mono_dc.rds")) {
  GEX_mono_dc <- readRDS("data/gex_mono_dc.rds")
  source("R/05_DE_monocytes_FindMarkers.R")
}

writeLines(capture.output(sessionInfo()), "results/sessionInfo.txt")
