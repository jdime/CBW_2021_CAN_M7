####################################
### Javier Diaz - javier.diazmejia@gmail.com
### Script based on:
### https://satijalab.org/seurat/articles/integration_introduction.html
####################################

####################################
### GENERAL OVERVIEW OF THIS SCRIPT
### 1) Loads scRNA-seq data from Cell Ranger (MTX files)
### 2) Generates QC plots for each dataset
### 3) Normalizes each dataset
### 4) Saves R objects with each normalized dataset
####################################

####################################
### HOW TO RUN THIS SCRIPT 
### 'Rscript ~/path_to/CBW_CAN_SingleCell_1_QC_Normalization.R'
####################################

####################################
###Install and check that Seurat v4.0 is available
####################################
# install.packages("Seurat")
library(Seurat)   # (CRAN) main scRNA-seq analysis package
library(bookdown) # (CRAN) to create a book of this script run

####################################
### Define input parameters
####################################
UserHomeDirectory <- Sys.getenv("HOME")[[1]]

StopWatchStart <- list()
StopWatchEnd   <- list()

StopWatchStart$Overall  <- Sys.time()

### Run specific inputs and parameters

### Select glioblastoma datasets from Richards_NatCancer_2021
DatasetIds <- list( ## List of input IDs to process
  `1` = "G1003_A_T",
  `2` = "G620_T",
  `3` = "G910_A_T",
  `4` = "G945_I_T",
  `5` = "G946_I_T",
  `6` = "G967_A_T",
  `7` = "G983_A_T"
)

### /path_to/MTX_DIRECTORIES from 10X Cell Ranger (i.e. 'filtered_feature_bc_matrix' directories
### with [features.tsv.gz, barcodes.tsv.gz and matrix.mtx.gz] files)
PathToDatasets <- list(
  "G1003_A_T" = "~/CourseData/CAN_data/Module7/MTX_FILES/Richards_NatCancer_2021/G1003_A_T/",
  "G620_T"    = "~/CourseData/CAN_data/Module7/MTX_FILES/Richards_NatCancer_2021/G620_T/",
  "G910_A_T"  = "~/CourseData/CAN_data/Module7/MTX_FILES/Richards_NatCancer_2021/G910_A_T/",
  "G945_I_T"  = "~/CourseData/CAN_data/Module7/MTX_FILES/Richards_NatCancer_2021/G945_I_T/",
  "G946_I_T"  = "~/CourseData/CAN_data/Module7/MTX_FILES/Richards_NatCancer_2021/G946_I_T/",
  "G967_A_T"  = "~/CourseData/CAN_data/Module7/MTX_FILES/Richards_NatCancer_2021/G967_A_T/",
  "G983_A_T"  = "~/CourseData/CAN_data/Module7/MTX_FILES/Richards_NatCancer_2021/G983_A_T/"
)

### Outputs
PathForOutfiles <- "~/CourseData/CAN_data/Module7/OUTFILES/Richards_NatCancer_2021" ## /path_to/out_directory
PrefixOutfiles  <- "Richards_NatCancer_2021"  ## Prefix for outfiles
PathForOutfiles <- gsub("^~/",paste0(UserHomeDirectory,"/"), PathForOutfiles)
dir.create(path = PathForOutfiles, recursive = T, showWarnings = F)

### Default parameters. Either suggested by Seurat developers, or tailored empirically.
DefaultParameters <- list(
  
  ### Parameters for QC plots
  CellPropertiesToQC = c("nFeature_RNA", "nCount_RNA", "mito.fraction"),
  
  ### Parameters for Seurat filters
  MinCells = 1,
  MinGenes = 50,
  MaxReads = 80000,
  MaxMitoFraction = 0.2
)

####################################
### Report R sessionInfo
####################################
writeLines("\n*** Report R sessionInfo ***\n")

OutfileRSessionInfo<-paste0(PathForOutfiles, "/", PrefixOutfiles, "_SingleCell_1_QC_Normalization_RSessionInfo.txt")
writeLines(capture.output(sessionInfo()), OutfileRSessionInfo)
capture.output(sessionInfo())

################################################################################################################################################
################################################################################################################################################
### HERE ARE THE FUNCTIONS TO LOAD AND QC DATASETS
################################################################################################################################################
################################################################################################################################################

writeLines("\n**** LOAD AND QC DATASETS ****\n")

####################################
###Load datasets from local (MTX files produced by Cell Ranger) and create Seurat objects
####################################

SeuratObjectsUnfiltered <-list()

for (datasetNumber in c(1:length(DatasetIds))) {
  dataset <- DatasetIds[[as.character(datasetNumber)]]
  DatasetPath <- PathToDatasets[[dataset]]
  DatasetPath<-gsub("^~/",paste0(UserHomeDirectory,"/"), DatasetPath)
  print(paste("Loading:", datasetNumber, dataset, DatasetPath, sep = " "))
  
  ### Create sparse matrix
  expression_matrix.mat <- Read10X(data.dir = DatasetPath, strip.suffix = T)
  colnames(expression_matrix.mat) <- paste(dataset, colnames(expression_matrix.mat), sep = "_") ### Add datasetID to cell-barcode
  dim(expression_matrix.mat)
  
  ### Create Seurat object
  seurat.object.u <- CreateSeuratObject(counts = expression_matrix.mat,
                                        min.cells = DefaultParameters$MinCells, 
                                        min.features = DefaultParameters$MinGenes, 
                                        project = paste0(PrefixOutfiles, "_", dataset)
                                        )
  seurat.object.u[['dataset.label']] <- dataset
  SeuratObjectsUnfiltered[[datasetNumber]]  <- seurat.object.u
}
print(SeuratObjectsUnfiltered)

####################################
###Filter cells based on number of genes, number of reads, and mitochondrial representation
####################################
SeuratObjectsFiltered <-list()

for (datasetNumber in c(1:length(DatasetIds))) {
  dataset <- DatasetIds[[as.character(datasetNumber)]]
  DatasetPath <- PathToDatasets[[dataset]]
  seurat.object.u <- SeuratObjectsUnfiltered[[datasetNumber]]
  print(paste("Filtering:", datasetNumber, dataset, DatasetPath, sep = " "))
  
  mitoRegExpressions <- paste(c("^MT-"), collapse = "|")
  mito.features <- grep(pattern = mitoRegExpressions, 
                        ignore.case = T, x = rownames(x = SeuratObjectsUnfiltered[[datasetNumber]]), value = T)

  ####################################
  ### Get mitochondrial genes
  ####################################
  writeLines(paste0("\n*** Get  mitochondrial genes for ", dataset, " ***\n"))

  if (length(mito.features)[[1]] > 0) {
    mito.fraction <- Matrix::colSums(x = GetAssayData(object = seurat.object.u, slot = 'counts')[mito.features, ]) / Matrix::colSums(x = GetAssayData(object = seurat.object.u, slot = 'counts'))
    seurat.object.u[['mito.fraction']] <- mito.fraction
  }else{
    mito.fraction <- 0 / Matrix::colSums(x = GetAssayData(object = seurat.object.u, slot = 'counts'))
    seurat.object.u[['mito.fraction']] <- mito.fraction
  }
  
  ####################################
  ### Generate violin plots
  ####################################
  writeLines(paste0("\n*** Generate violin plots ***\n"))
  
  VlnPlotPdf<-paste0(PathForOutfiles, "/", PrefixOutfiles, "_QC_ViolinPlots_", dataset, ".pdf")
  pdf(file=VlnPlotPdf, width = 10, height = 5)
  print(VlnPlot(seurat.object.u, features = c("nFeature_RNA", "nCount_RNA", "mito.fraction"), ncol = 3))
  dev.off()

  ####################################
  ### Filter cells based gene counts, number of reads, ribosomal and mitochondrial representation
  ####################################
  writeLines(paste0("\n*** Filter cells based on number of genes, number of reads, and mitochondrial representation for ", dataset, " ***\n"))
  
  if (length(mito.features)[[1]] > 0) {
    seurat.object.f<-subset(x = seurat.object.u, subset = 
                              nFeature_RNA  >= DefaultParameters$MinGenes
                            & nCount_RNA    <= DefaultParameters$MaxReads
                            & mito.fraction <= DefaultParameters$MaxMitoFraction)
    
  }else{
    seurat.object.f<-subset(x = seurat.object.u, subset = 
                              nFeature_RNA  >= DefaultParameters$MinGenes
                            & nCount_RNA    <= DefaultParameters$MaxReads)
  }
  
  SeuratObjectsFiltered[[datasetNumber]]  <- seurat.object.f

}

####################################
### Merge Seurat objects RNA assay
####################################
writeLines("\n*** Merge Seurat objects RNA assay ***\n")

FirstSeuratObject   <- SeuratObjectsFiltered[[1]]
RestOfSeuratObjectsFiltered <- SeuratObjectsFiltered[c(2:length(DatasetIds))]
RestOfDatasetsIds    <- unlist(DatasetIds[c(2:length(DatasetIds))])

seurat.object.merged <- merge(FirstSeuratObject, y = RestOfSeuratObjectsFiltered, project = PrefixOutfiles)
seurat.object.merged <- AddMetaData(object = seurat.object.merged, metadata = seurat.object.merged@meta.data$dataset.label, col.name = "dataset")
seurat.object.list <- SplitObject(seurat.object.merged, split.by = "dataset.label")

print(seurat.object.list)

################################################################################################################################################
################################################################################################################################################
### HERE ARE THE FUNCTIONS TO NORMALIZE DATASETS
################################################################################################################################################
################################################################################################################################################
writeLines("\n**** NORMALIZE DATASETS ****\n")

####################################
### Running SCTransform
####################################
writeLines("\n*** Running SCTransform ***\n")

for (i in 1:length(seurat.object.list)) {
  dataset <- names(seurat.object.list)[[i]]
  print(paste0("Normalizing: ", dataset))
  seurat.object.list[[i]] <- SCTransform(seurat.object.list[[i]], verbose = T)
}

################################################################################################################################################
################################################################################################################################################
### HERE ARE THE FUNCTIONS TO SAVE EACH DATASET R_OBJECT
################################################################################################################################################
################################################################################################################################################

writeLines("\n**** SAVE EACH DATASET R_OBJECT ****\n")

####################################
### Saving each dataset R object
####################################

writeLines("\n*** Saving each dataset R object ***\n")

for (i in 1:length(seurat.object.list)) {
  dataset <- names(seurat.object.list)[[i]]
  OutfileRDS<-paste0(PathForOutfiles, "/", PrefixOutfiles, "_", dataset , "_QC_Normalization.rds")
  saveRDS(seurat.object.list[[i]], file = OutfileRDS)
}

####################################
### Obtain computing time used
####################################
writeLines("\n*** Obtain computing time used ***\n")

StopWatchEnd$Overall  <- Sys.time()

OutfileCPUtimes<-paste0(PathForOutfiles, "/", PrefixOutfiles, "_QC_Normalization.txt")

Headers<-paste("Step", "Time(minutes)", sep="\t")
write.table(Headers, file = OutfileCPUtimes, row.names = F, col.names = F, sep="\t", quote = F, append = T)

lapply(names(StopWatchStart), function(stepToClock) {
  if (regexpr("POSIXct", class(StopWatchStart[[stepToClock]]), ignore.case = T)[1] == 1) {
    TimeStart <- StopWatchStart[[stepToClock]]
    TimeEnd   <- StopWatchEnd[[stepToClock]]
    TimeDiff <- format(difftime(TimeEnd, TimeStart, units = "min"))
    ReportTime<-c(paste(stepToClock, TimeDiff, sep = "\t", collapse = ""))
    write(file = OutfileCPUtimes, x=gsub(pattern = " mins", replacement = "", x = ReportTime), append = T)
  }
})

####################################
### Finish
####################################
writeLines(paste0("\nEND - Check:\n", PathForOutfiles, "\nFor outfiles\n\n"))

quit()
