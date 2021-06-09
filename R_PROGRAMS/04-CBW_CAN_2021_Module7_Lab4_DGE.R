####################################
### Javier Diaz - javier.diazmejia@gmail.com
### Script based on:
### https://satijalab.org/seurat/articles/integration_introduction.html
####################################

####################################
### GENERAL OVERVIEW OF THIS SCRIPT
### 1) Loads integrated datasets R object produced by script `CBW_CAN_SingleCell_3_PCA_Clustering_DimReduction.R`
### 2) Computes differential gene expression (DGE) for each cell-class vs. the rest of cells
### 3) Computes DGE for one cell-class vs. another cell-class
####################################

####################################
### HOW TO RUN THIS SCRIPT
### Using one-line-commands in a console or terminal type:
### 'Rscript ~/path_to/CBW_CAN_SingleCell_4_DGE.R'
####################################

####################################
### Tested with R v4.0.2
####################################

####################################
###Install and check that Seurat v4.0 is available
####################################
library(Seurat)  # (CRAN) main scRNA-seq analysis package
library(future)  # (CRAN) to run parallel processes

####################################
### Define inputs, parameters and outputs
####################################
UserHomeDirectory <- Sys.getenv("HOME")[[1]]

StopWatchStart <- list()
StopWatchEnd   <- list()

StopWatchStart$Overall  <- Sys.time()

### Inputs
AnchorFinder        <- "STACAS"
PathToIntegratedRds <- paste0("~/CourseData/CAN_data/Module7/OUTFILES/Richards_NatCancer_2021/", AnchorFinder, "/WITHOUT_G945_I_T", "/Richards_NatCancer_2021_PCA_Clustering_DimReduction.rds") ### /path_to/*rds (R object) from CBW_CAN_SingleCell_3_PCA_Clustering_DimReduction.R
InfileMetadata      <- "~/CourseData/CAN_data/Module7/METADATA/Richards_NatCancer_2021_without_G945_I_T.metadata.tsv"
PathToIntegratedRds <- gsub("^~/",paste0(UserHomeDirectory,"/"), PathToIntegratedRds)
ASSAY               <- "SCT" ### Note, there is debate in the field about whether DGE should be calculated using RNA or SCT

### Outputs
PathForOutfiles     <- paste0("~/CourseData/CAN_data/Module7/OUTFILES/Richards_NatCancer_2021/", AnchorFinder, "/WITHOUT_G945_I_T") ## /path_to/out_directory
PrefixOutfiles      <- "Richards_NatCancer_2021"  ## Prefix for outfiles
PathForOutfiles     <- gsub("^~/",paste0(UserHomeDirectory,"/"), PathForOutfiles)
dir.create(path = PathForOutfiles, recursive = T, showWarnings = F)

### Default parameters. Either suggested by Seurat developers, or tailored empirically.
DefaultParameters <- list(
  
  ### Parameters for DGE calculation OneVsRest
  Test_OneVsRest               = "wilcox",
  OnlyPos_OneVsRest            = F,
  ReturnMinPctThresh_OneVsRest = 0,  # Use 0 if all cells should be included
  ReturnPvalThresh_OneVsRest   = 1,  # Use 1 if NO p-value cutoff should be used
  ReturnLogFcThresh_OneVsRest  = 0,  # Use 0 if NO LogFc cutoff should be used

  ### Parameters for DGE calculation Paired
  Test_Paired                  = "wilcox",
  OnlyPos_Paired               = F,
  ReturnMinPctThresh_Paired    = 0,  # Use 0 if all cells should be included
  ReturnPvalThresh_Paired      = 1,  # Use 1 if NO p-value cutoff should be used
  ReturnLogFcThresh_Paired     = 0,  # Use 0 if NO LogFc cutoff should be used
  subclass1_Paired             = "malignant",
  subclass2_Paired             = "oligodendrocyte"
)

####################################
### Report R sessionInfo
####################################
writeLines("\n*** Report R sessionInfo ***\n")

OutfileRSessionInfo<-paste0(PathForOutfiles, "/", PrefixOutfiles, "_SingleCell_4_DGE_RSessionInfo.txt")
writeLines(capture.output(sessionInfo()), OutfileRSessionInfo)
capture.output(sessionInfo())

####################################
### Define number of cores and RAM for parallelization
####################################

NumbCores <- "MAX"
MaxGlobalVariables <- 30000

if (regexpr("^MAX$", NumbCores, ignore.case = T)[1] == 1) {
  NumbCoresToUse <- availableCores()[[1]]
}else if (regexpr("^[0-9]+$", NumbCores, ignore.case = T)[1] == 1) {
  NumbCoresToUse <- as.numeric(NumbCores)
}else{
  stop(paste0("Unexpected format for NumbCores: ", NumbCores, "\n"))
}

if (NumbCoresToUse == 1) {
  plan(strategy = "sequential")
  writeLines(paste0("\n", "*** Running: in 'sequential' mode with ", NumbCoresToUse, " core ***", "\n"))
}else if (NumbCoresToUse > 1) {
  plan(strategy = "multicore", workers = NumbCoresToUse)
  writeLines(paste0("\n", "*** Running: in 'multicore' mode with ", NumbCoresToUse, " cores ***", "\n"))
}else{
  stop(paste0("Unexpected NumbCores = ", NumbCoresToUse))
}

options(future.globals.maxSize = MaxGlobalVariables * 1024^2)

################################################################################################################################################
################################################################################################################################################
### HERE ARE THE FUNCTIONS TO LOAD DATASETS
################################################################################################################################################
################################################################################################################################################

writeLines("\n**** LOAD DATASETS ****\n")

####################################
### Load integrated datasets R object
####################################
writeLines("\n*** Load integrated datasets R object ***\n")

seurat.object.integrated <- readRDS(PathToIntegratedRds)

####################################
### Loading metadata from --infile_metadata
####################################
writeLines("\n*** Load cell-level metadata ***\n")

CellPropertiesFromMetadata <- data.frame(read.table(InfileMetadata, header = T, row.names = 1, check.names = F))
head(CellPropertiesFromMetadata)

################################################################################################################################################
################################################################################################################################################
### HERE ARE THE FUNCTIONS TO COMPUTE DGE USING GLOBAL CLUSTERS
################################################################################################################################################
################################################################################################################################################

writeLines("\n**** COMPUTE DGE USING GLOBAL CLUSTERS ****\n")

####################################
### Finding differentially expressed genes using global cell clusters, compares each cell cluster vs. the rest of cells
####################################

writeLines("\n*** Finding differentially expressed genes using global cell clusters, compares each cell cluster vs. the rest of cells ***\n")

Idents(object = seurat.object.integrated) <- "seurat_clusters"

### For details on this approach to use a pseudocount, see https://f1000research.com/articles/7-1522
FindMarkers.Pseudocount  <- 1/length(rownames(seurat.object.integrated@meta.data))

###  base = exp(1) is needed in Seurat v4 because the default uses Log2FC instead of LogNatFC, like v3 does
seurat.object.integrated.markers <- FindAllMarkers(object = seurat.object.integrated,
                                                   assay = ASSAY,
                                                   test.use = DefaultParameters$Test_OneVsRest,
                                                   only.pos = DefaultParameters$OnlyPos_OneVsRest,
                                                   min.pct =  DefaultParameters$ReturnMinPctThresh_OneVsRest,
                                                   return.thresh =   DefaultParameters$ReturnPvalThresh_OneVsRest,
                                                   logfc.threshold = DefaultParameters$ReturnLogFcThresh_OneVsRest,
                                                   pseudocount.use = FindMarkers.Pseudocount,
                                                   base = exp(1))

### IMPORTANT: use column 'gene' instead of rownames
head(seurat.object.integrated.markers)

SimplifiedDiffExprGenes.df <- seurat.object.integrated.markers[,c("cluster","gene","p_val","p_val_adj","avg_logFC","pct.1","pct.2")]

### Add -Log10(Pval) * sign of FC
SimplifiedDiffExprGenes.df[["mLog10Pval_FCsign"]] <- (log10(SimplifiedDiffExprGenes.df[["p_val"]])*-1) * sign(SimplifiedDiffExprGenes.df[,"avg_logFC"])
### Replace -Inf and Inf in mLog10Pval_FCsign
SimplifiedDiffExprGenes.df[,"mLog10Pval_FCsign"][SimplifiedDiffExprGenes.df[,"mLog10Pval_FCsign"] ==  Inf]  <- log10(min(SimplifiedDiffExprGenes.df[,"p_val"][SimplifiedDiffExprGenes.df[,"p_val"] > 0]))*-1
SimplifiedDiffExprGenes.df[,"mLog10Pval_FCsign"][SimplifiedDiffExprGenes.df[,"mLog10Pval_FCsign"] ==  -Inf] <- log10(max(SimplifiedDiffExprGenes.df[,"p_val"][SimplifiedDiffExprGenes.df[,"p_val"] > 0]))*-1
head(SimplifiedDiffExprGenes.df)

OutfileDGE <- paste0(PathForOutfiles, "/", PrefixOutfiles, "_DGE_GlobalClustering_", ASSAY, "_", DefaultParameters$Test_OneVsRest, ".tsv.bz2")
Outfile.con <- bzfile(OutfileDGE, "w")
write.table(x = data.frame(SimplifiedDiffExprGenes.df), file = Outfile.con, row.names = F, sep="\t", quote = F)
close(Outfile.con)

####################################
### Finding differentially expressed genes using metadata annotations, compares two cell classes vs. each other
####################################
writeLines("\n*** Finding differentially expressed genes using metadata annotations, compares two cell classes vs. each other ***\n")

Idents(seurat.object.integrated) <- "cell_type"
FindMarkers.Pseudocount  <- 1/length(rownames(seurat.object.integrated@meta.data))
###  base = exp(1) is needed in Seurat v4 because the default uses Log2FC instead of LogNatFC, like v3 does
seurat.object.integrated.markers <- data.frame(FindMarkers(object = seurat.object.integrated,
                                                           assay = ASSAY,
                                                           test.use = DefaultParameters$Test_Paired,
                                                           only.pos = DefaultParameters$OnlyPos_Paired,
                                                           ident.1  = DefaultParameters$subclass1_Paired,
                                                           ident.2  = DefaultParameters$subclass2_Paired,
                                                           min.pct  = DefaultParameters$ReturnMinPctThresh_Paired,
                                                           return.thresh   = DefaultParameters$ReturnPvalThresh_Paired,
                                                           logfc.threshold = DefaultParameters$ReturnLogFcThresh_Paired,
                                                           pseudocount.use = FindMarkers.Pseudocount,
                                                           base = exp(1)))
seurat.object.integrated.markers$class1 <- DefaultParameters$subclass1_Paired
seurat.object.integrated.markers$class2 <- DefaultParameters$subclass2_Paired
seurat.object.integrated.markers$gene   <- rownames(seurat.object.integrated.markers)

SimplifiedDiffExprGenes.df <- seurat.object.integrated.markers[,c("class1","class2","gene","p_val","p_val_adj","avg_logFC","pct.1","pct.2")]

### Add -Log10(Pval) * sign of FC
SimplifiedDiffExprGenes.df[["mLog10Pval_FCsign"]] <- (log10(SimplifiedDiffExprGenes.df[["p_val"]])*-1) * sign(SimplifiedDiffExprGenes.df[,"avg_logFC"])
### Replace -Inf and Inf in mLog10Pval_FCsign
SimplifiedDiffExprGenes.df[,"mLog10Pval_FCsign"][SimplifiedDiffExprGenes.df[,"mLog10Pval_FCsign"] ==  Inf]  <- log10(min(SimplifiedDiffExprGenes.df[,"p_val"][SimplifiedDiffExprGenes.df[,"p_val"] > 0]))*-1
SimplifiedDiffExprGenes.df[,"mLog10Pval_FCsign"][SimplifiedDiffExprGenes.df[,"mLog10Pval_FCsign"] ==  -Inf] <- log10(max(SimplifiedDiffExprGenes.df[,"p_val"][SimplifiedDiffExprGenes.df[,"p_val"] > 0]))*-1
head(SimplifiedDiffExprGenes.df)

OutfileDGE <- paste0(PathForOutfiles, "/", PrefixOutfiles, "_DGE_Paired_", ASSAY, "_", DefaultParameters$Test_Paired, ".tsv.bz2")
Outfile.con <- bzfile(OutfileDGE, "w")
write.table(x = data.frame(SimplifiedDiffExprGenes.df), file = Outfile.con, row.names = F, sep="\t", quote = F)
close(Outfile.con)

####################################
### Obtain computing time used
####################################
writeLines("\n*** Obtain computing time used ***\n")

StopWatchEnd$Overall  <- Sys.time()

OutfileCPUtimes<-paste0(PathForOutfiles, "/", PrefixOutfiles, "_DGE_CPUtimes.txt")
write(file = OutfileCPUtimes, x = paste("#Number_of_cores_used", NumbCoresToUse, sep = "\t", collapse = ""))
write(file = OutfileCPUtimes, x = paste("#MaxGlobalVariables", MaxGlobalVariables, sep = "\t", collapse = ""), append = T)

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

