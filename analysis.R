# Loading libraries
library(DESeq2)
library(airway)
library(ggplot2)
library(pheatmap)
library(ggrepel)
library(dplyr)
library(tibble)

# Loading airway dataset
# Human airway cells treated with dexamethasone (a steroid) vs untreated control cells
data("airway")

# Converting to DESeq2 object
dds <- DESeqDataSet(airway, design = ~ dex)

# Removing genes with very low counts
dds <- dds[rowSums(counts(dds)) >= 10, ]

cat("Dataset loaded:", nrow(dds), "genes across", ncol(dds), "samples\n")

# Running DESeq2
dds <- DESeq(dds)

# Untreated vs treated with dexamethasone
res <- results(dds, contrast = c("dex", "trt", "untrt"))

# Ordering by adjusted p-value
res_ordered <- res[order(res$padj), ]

# Saving as csv
write.csv(as.data.frame(res_ordered), file = "deseq2_results.csv")
summary(res)

# Plot 1: Volcano
res_df <- as.data.frame(res) %>%
  rownames_to_column("gene") %>%
  filter(!is.na(padj)) %>%
  mutate(
    significance = case_when(
      padj < 0.05 & log2FoldChange > 1  ~ "Up",
      padj < 0.05 & log2FoldChange < -1 ~ "Down",
      TRUE ~ "Not significant"
    )
  )

volcano <- ggplot(res_df, aes(x = log2FoldChange, y = -log10(padj), color = significance)) +
  geom_point(alpha = 0.6, size = 1.5) +
  scale_color_manual(values = c("Up" = "#D85A30", "Down" = "#378ADD", "Not significant" = "#888780")) +
  geom_hline(yintercept = -log10(0.05), linetype = "dashed", color = "gray50") +
  geom_vline(xintercept = c(-1, 1), linetype = "dashed", color = "gray50") +
  labs(title = "Differential Gene Expression: Treated vs Untreated",
       x = "Log2 Fold Change", y = "-Log10 Adjusted P-value",
       color = "Gene status") +
  theme_minimal()

ggsave("volcano_plot.png", volcano, width = 8, height = 6, dpi = 150)
cat("Saved volcano_plot.png\n")

# Plot 2: Heatmap Visualisation

vsd <- vst(dds, blind = FALSE)  # variance-stabilising transform

top_genes <- head(order(res$padj, na.last = TRUE), 30)
mat <- assay(vsd)[top_genes, ]
mat <- mat - rowMeans(mat)  # centre each gene

png("heatmap.png", width = 800, height = 900, res = 120)
pheatmap(mat,
         annotation_col = as.data.frame(colData(vsd)[, "dex"]),
         main = "Top 30 Differentially Expressed Genes",
         fontsize_row = 8)
dev.off()
cat("Saved heatmap.png\n")

# Plot 3: PCA

pca_data <- plotPCA(vsd, intgroup = "dex", returnData = TRUE)
percent_var <- round(100 * attr(pca_data, "percentVar"))

pca_plot <- ggplot(pca_data, aes(x = PC1, y = PC2, color = dex)) +
  geom_point(size = 4) +
  xlab(paste0("PC1: ", percent_var[1], "% variance")) +
  ylab(paste0("PC2: ", percent_var[2], "% variance")) +
  labs(title = "PCA of RNA-seq Samples", color = "Treatment") +
  scale_color_manual(values = c("trt" = "#D85A30", "untrt" = "#378ADD")) +
  theme_minimal()

ggsave("pca_plot.png", pca_plot, width = 7, height = 5, dpi = 150)
cat("Saved pca_plot.png\n")
