## CANADIAN BIOINFORMATICS WORKSHOPS 2021 - CANCER ANALYSIS - MODULE 7 (SINGLE CELL RNA)

This tutorial guides the user to a full single cell RNA-seq (scRNA-seq) analysis, starting from raw cell-vs-gene UMI counts matrices from Cell Ranger. It includes four hands-on labs (1 to 4) provided as Rmd files, and two extra tutorials (5 and 6) provided as R scripts:
1) Data QC and normalization
2) Dataset integration and batch effect correction
3) Dimension reduction, PCA, UMAP and cell clustering
4) Differential gene expression
5) cell-cluster labeling tutorial with GSVA
6) Copy number variant (CNV) identification with InferCNV

The first four labs use mainly the R packages Seurat and STACAS; whereas the two extra tutorials use GSVA and InferCNV.
Package versions are provided in the Git pages: https://jdime.github.io/CBW_2021_CAN_M7

For this module, we use seven glioblastoma scRNA-seq samples from [Richards et al 2021 (Nature Cancer), authors Figure 7](https://www.nature.com/articles/s43018-020-00154-9) as an example.

Raw scRNA-seq measurements are provided in MTX format, metadata, and code needed for this tutorial can be obtained from Zenodo:
https://zenodo.org/record/4913828#.YMAr3TZKjUI

