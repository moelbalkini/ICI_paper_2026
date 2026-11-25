# -----------------------------------------------------------------------------
# 05_DE_monocytes_FindMarkers.R
# Cell-level differential expression in monocytes/DCs of the n-irAE dataset
# (Giguelay et al. 2025), n-irAE_post vs Control_post (Figure 3E).
#   - targeted panel: Wilcoxon test (Seurat::FindMarkers) + violin plots
#   - transcriptome-wide: MAST (Seurat::FindMarkers) + volcano plot
# Input: Seurat object GEX_mono_dc (monocytes/DCs) with meta.data columns
#        Dataset and Group_Timepoint.
# -----------------------------------------------------------------------------

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
