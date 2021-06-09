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

Raw scRNA-seq measurements are provided in MTX format, metadata, relevant intermediate outfiles, and code needed for this tutorial can be obtained from Zenodo:
https://zenodo.org/record/4918244#.YMEL9zZKjUI

### How to run scripts

*Scripts for Labs 1 to 4 are provide in two versions, pick the version that you feel more comfortable with:*<br />
a) \*Rmd files can be rendered using Rstudio<br />
b) R_PROGRAMS/\*R files can be run in a console/terminal like:<br />
`Rscript ~/path_to/CBW_CAN_2021_Module7_Lab1_QC_Normalization.R`<br />
`Rscript ~/path_to/CBW_CAN_2021_Module7_Lab2_Integration.R`<br />
`Rscript ~/path_to/CBW_CAN_2021_Module7_Lab3_PCA_Clustering_DimReduction.R`<br />
`Rscript ~/path_to/CBW_CAN_2021_Module7_Lab4_DGE.R`<br />

*Scripts for Labs 5 and 6 must be run in a console/terminal using one-line-commands, like:*<br />
`Rscript /path_to/CBW_CAN_SingleCell_5_GSVA.R -i /path_to/Richards_NatCancer_2021_AverageGeneExpression_GlobalClustering.tsv.bz2 -t DGE -c /path_to/LM22_signature.cutoff3000.UsedF1000paper.renamed.gmt -o /path_to/WITHOUT_G945_I_T -p Richards_NatCancer_2021_AverageGeneExpression_GlobalClustering -e 0.05 -f 0.1`<br />

`Rscript /path_to/CBW_CAN_SingleCell_6_InferCNV.R -i /path_to/ -t MTX -j /path_to/Richards_NatCancer_2021_GlobalClustering_CellClusters_G983_A_T.tsv.bz2 -k 4 -g /path_to/gencode_v19_gene_pos.txt -m 0.1 -n 0.1 -s 0.15 -o /path_to/ -p Richards_NatCancer_2021_GlobalClustering_CellClusters_G983_A_T -u MAX -l N -a 10000`<br />
<br />
*For help with Labs 5 and 6 scripts use:*<br />
`Rscript /path_to/CBW_CAN_SingleCell_5_GSVA.R -h`<br />
`Rscript /path_to/CBW_CAN_SingleCell_6_InferCNV.R -h`<br />
