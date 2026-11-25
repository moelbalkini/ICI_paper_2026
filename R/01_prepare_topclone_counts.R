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
