## CANADIAN BIOINFORMATICS WORKSHOPS 2021 - CANCER ANALYSIS - MODULE 7 (SINGLE CELL RNA)

This tutorial guides the user to a full single cell RNA-seq (scRNA-seq) analysis, starting from raw cell-vs-gene UMI counts matrices from Cell Ranger. It includes four hands-on labs:
1) Data QC and normalization
2) Dataset integration and batch effect correction
3) Dimension reduction, PCA, UMAP and cell clustering
4) Differential gene expression
These four labs use mainly the R packages Seurat and STACAS

Additionally, two tutorial inputs and code are included to do:
a) cell-cluster labeling tutorial with GSVA
b) Copy number variant (CNV) identification with InferCNV

For this module, we use seven glioblastoma scRNA-seq samples from [Richards et al 2021 (Nature Cancer), authors Figure 7](https://www.nature.com/articles/s43018-020-00154-9) as an example.

Raw scRNA-seq measurements are provided in MTX format, metadata, and code needed for this tutorial can be obtained from Zenodo: https://doi.org/10.5281/zenodo.4913828

