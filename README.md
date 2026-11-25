# Clonal representation and differential expression in ICI-induced n-irAE

Code for the patient-level beta-binomial mixed-effects model of clonal representation
in:

> Elbalkini M, Thomas J, et al. *Neurotoxic immune-related adverse events following
> immune checkpoint blockade are associated with aberrant cytotoxic TEMRA T cells and
> serum cytokine responses.* Neuro-Oncology Advances (under revision).

## What the model does

Data: post-treatment PBMC scRNA-seq + TCR data from Giguelay et al. 2025 (GEO
GSE278119; controls n = 9, n-irAE n = 8).

1. For each patient, the 10 most expanded TCR clonotypes are selected.
2. For each clone and each CD4⁺ (C1–C8) or CD8⁺ (C1–C6) cluster, the cells of the
   clone inside (`k`) and outside the cluster are counted.
3. Per cluster, a beta-binomial generalized linear mixed model is fitted with glmmTMB:

   ```
   cbind(k, outside_cluster) ~ Group + offset(logit(K/N)) + (1 | SampleID)
   ```

   `K/N` is the patient's own frequency of that cluster, so the group effect is the
   odds ratio (n-irAE vs control) that an expanded clone is over-represented in the
   cluster relative to the patient background. Patient is a random intercept;
   the beta-binomial family accounts for overdispersion between clones.
4. Wald 95% confidence intervals; Benjamini-Hochberg FDR within each lineage and
   across all 14 clusters. Convergence is checked via a positive-definite Hessian;
   residuals are checked with DHARMa.

The CD4⁺ TEMRA comparison was pre-specified; all other clusters are exploratory.

## Differential expression

Differential expression in monocytes/DCs between n-irAE and control post-treatment
samples (Fig. 3E) is tested at the cell level with Seurat `FindMarkers`: a Wilcoxon
test for a pre-defined inflammatory gene panel, and MAST transcriptome-wide
(Bonferroni-adjusted p-values; |log2FC| > 0.25 and adjusted p < 0.05 shown as
significant in the volcano plot).

## Files

| File | Content |
|---|---|
| `run_analysis.R` | runs the full analysis and writes `results/` |
| `clonal_and_DE_complete.R` | the complete analysis as one script, as run for the paper |
| `R/01_prepare_topclone_counts.R` | builds the clone × cluster count table |
| `R/02_fit_betabinomial_models.R` | per-cluster beta-binomial GLMM, odds ratios, FDR |
| `R/03_model_diagnostics.R` | binomial vs beta-binomial AIC, DHARMa diagnostics |
| `R/04_plots.R` | forest plot (Fig. 3B) and patient-level ratios (Supp. Fig. 7) |
| `R/05_DE_monocytes_FindMarkers.R` | Wilcoxon gene panel, violin plots, MAST and volcano plot (Fig. 3E) |

## Input

Two Seurat objects (`data/gex_CD4_post.rds`, `data/gex_CD8_post.rds`) with re-clustered
post-treatment CD4⁺ / CD8⁺ T cells, and in `meta.data` the columns `SampleID`,
`Group` (`"Control"` / `"n-irAE"`), `clonotype` and `seurat_clusters`. They are
derived from the public GEO data using the preprocessing in the manuscript Methods
(SoupX, Seurat, DoubletFinder, Harmony, Azimuth reference mapping). No patient data
from the authors' own cohort are included.

For Fig. 3E, a Seurat object of monocytes/DCs (`data/gex_mono_dc.rds`) with the
`meta.data` columns `Dataset` and `Group_Timepoint` (`"Control_post"` / `"n-irAE_post"`).

## Requirements

R ≥ 4.5 with `Seurat`, `dplyr`, `tidyr`, `ggplot2`, `glmmTMB`, `DHARMa`, `ggpubr`, `ggrepel`, `MAST` (Bioconductor).

```r
install.packages(c("dplyr", "tidyr", "ggplot2", "glmmTMB", "DHARMa", "Seurat", "ggpubr", "ggrepel"))
BiocManager::install("MAST")
source("run_analysis.R")
```

## License

MIT (see `LICENSE`).
