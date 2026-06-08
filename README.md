# RNA-seq Differential Expression Analysis

Analysing the airway dataset in R using DESeq2

## Brief Summary

This short project compares the gene expression in human airway muscle cells in cells that have been treated with dexamethasone against a control group that has not been treated. This analysis identifies up- and downregulated genes by creating various visualisations: volcano plot, heatmap, and PCA plot.

DESeq was used for the differential expression, the ggplot2 package was used for the volcano and PCA plots, and the pheatmap package was used to create the heatmap.

## Outputs

- `deseq2_results.csv` — full results table
- `volcano_plot.png` — volcano plot of all genes
- `heatmap.png` — top 30 DE genes
- `pca_plot.png` — sample clustering by treatment

## Reference for Data Used
`airway` Bioconductor package (Himes et al. 2014, PMID 24926665)

