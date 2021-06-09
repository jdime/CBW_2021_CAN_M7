####################################
### Javier Diaz - javier.diazmejia@gmail.com
### Script based on:
### https://satijalab.org/seurat/articles/integration_introduction.html
####################################

####################################
### GENERAL OVERVIEW OF THIS SCRIPT
### 1) Loads each normalized dataset R object produced by script `CBW_CAN_2021_Module7_Lab1__QC_Normalization'
### 2) Merges and integrates datasets correcting batch effects
### 3) Saves R object with integrated datasets
####################################

####################################
### HOW TO RUN THIS SCRIPT 
### 'Rscript ~/path_to/CBW_CAN_SingleCell_2_Integration.R'
####################################

####################################
###Install and check that Seurat v4.0 is available
####################################
library(Seurat) # (CRAN) main scRNA-seq analysis package
library(future) # (CRAN) to run parallel processes
library(STACAS) # (GitHub carmonalab/STACAS) for STACAS-based anchor detection

####################################
### Define inputs, parameters and outputs
####################################
UserHomeDirectory <- Sys.getenv("HOME")[[1]]

StopWatchStart <- list()
StopWatchEnd   <- list()

StopWatchStart$Overall  <- Sys.time()

### Run specific inputs and parameters

AnchorFinder        <- "STACAS" ## Either 'STACAS' or 'SEURAT'

### /path_to/*rds infiles (R objects) from CBW_CAN_2021_Module7_Lab1__QC_Normalization
PathToDatasets <- list(
  "G1003_A_T" = "~/CourseData/CAN_data/Module7/OUTFILES/Richards_NatCancer_2021/Richards_NatCancer_2021_G1003_A_T_QC_Normalization.rds",
  # "G945_I_T"  = "~/CourseData/CAN_data/Module7/OUTFILES/Richards_NatCancer_2021/Richards_NatCancer_2021_G945_I_T_QC_Normalization.rds",
  "G910_A_T"  = "~/CourseData/CAN_data/Module7/OUTFILES/Richards_NatCancer_2021/Richards_NatCancer_2021_G910_A_T_QC_Normalization.rds",
  "G983_A_T"  = "~/CourseData/CAN_data/Module7/OUTFILES/Richards_NatCancer_2021/Richards_NatCancer_2021_G983_A_T_QC_Normalization.rds",
  "G620_T"    = "~/CourseData/CAN_data/Module7/OUTFILES/Richards_NatCancer_2021/Richards_NatCancer_2021_G620_T_QC_Normalization.rds",
  "G946_I_T"  = "~/CourseData/CAN_data/Module7/OUTFILES/Richards_NatCancer_2021/Richards_NatCancer_2021_G946_I_T_QC_Normalization.rds",
  "G967_A_T"  = "~/CourseData/CAN_data/Module7/OUTFILES/Richards_NatCancer_2021/Richards_NatCancer_2021_G967_A_T_QC_Normalization.rds"
)

### Outputs
PathForOutfiles     <- paste0("~/CourseData/CAN_data/Module7/OUTFILES/Richards_NatCancer_2021/", AnchorFinder, "/WITHOUT_G945_I_T") ## /path_to/out_directory
PrefixOutfiles      <- "Richards_NatCancer_2021"  ## Prefix for outfiles
PathForOutfiles     <- gsub("^~/",paste0(UserHomeDirectory,"/"), PathForOutfiles)
dir.create(path = PathForOutfiles, recursive = T, showWarnings = F)

### Default parameters. Either suggested by Seurat developers, or tailored empirically.

### NOTE: when datasets being integrated have fewer anchors than needed for the integration,
### step IntegrateData() will produce an error:
### `Error in nn2(data = c(15.4534704633918, 15.4534704633918, 15.4534704633918,  : 
### Cannot find more nearest neighbours than there are points`
### To fix this using Seurat's anchor finder in step FindIntegrationAnchors() lower the k.filter value (e.g. 150, default = 200)
### To fix this using STACAS anchor finder, in step FilterAnchors.STACAS() increase the dist.pct and dist.thr values (e.g. 0.9, default = 0.8)
### More info: https://github.com/satijalab/seurat/issues/2056
DefaultParameters <- list(
  
  ### Parameters for dataset integration
  IntegrationNFeatures = 3000,
  StacasVarGenesIntegratedN = 500, #only needed if using STACAS
  PcaDimsUse = c(1:20),
  SeuratFindAnchorsKFilter = 200,
  StacasDistPct = 0.9,
  StacasDistThr = 0.9,
  
  ### ReferenceDatasets = either a <comma> delimited list of dataset ID(s) to be used as
  ### references for anchors or 'NA' to compare all-vs-all
  ReferenceDatasets = "NA"
  
)

####################################
### Report R sessionInfo
####################################
writeLines("\n*** Report R sessionInfo ***\n")

OutfileRSessionInfo<-paste0(PathForOutfiles, "/", PrefixOutfiles, "_SingleCell_2_Integration_RSessionInfo.txt")
writeLines(capture.output(sessionInfo()), OutfileRSessionInfo)
capture.output(sessionInfo())

####################################
### Define number of cores and RAM for parallelization
####################################

NumbCores <- 1
MaxGlobalVariables <- 10000

### Using plan(strategy = "multisession", workers = NumbCoresToUse)
### works out in Javier's laptop and HPC systems, but it fails in the CBW AWS instance,
### even setting up `MaxGlobalVariables <- 50000`, the run produced an error in the IntegrateData() step, like:
###   Failed to retrieve the result of MulticoreFuture (future_lapply-1) from the forked worker (on localhost; PID 194469).
###   Post-mortem diagnostic: No process exists with this PID, i.e. the forked localhost worker is no longer alive.
### This seems to be related to not having enough RAM in the AWS instance allocated for the course
### To fix this for the CBW CAN, we are:
###   a) removing from the analysis a dataset (G945_I_T) with high mito.fraction
###   b) Using `plan(strategy = "sequential")` in script CBW_CAN_SingleCell_2_Integration.R (i.e. no parallelization)
### 
### Notes: `options(future.globals.maxSize = MaxGlobalVariables)` sets the maximum allowed size in bytes.
###         So to set it to 10GB, you would run options(future.globals.maxSize = 10000 * 1024^2).
###         Increasing the value of future.globals.maxSize will increase your RAM usage so set this number mindfully.
###         https://satijalab.org/seurat/archive/v3.1/future_vignette.html
###         https://github.com/satijalab/seurat/issues/3249

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
  # plan(strategy = "multicore", workers = NumbCoresToUse)
  plan(strategy = "multisession", workers = NumbCoresToUse)
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
### Load R Objects
####################################

writeLines("\n*** Load R Objects ***\n")

InputsTable <- data.frame(matrix(unlist(PathToDatasets), nrow=length(PathToDatasets), byrow=TRUE))
rownames(InputsTable) <- names(PathToDatasets)
colnames(InputsTable)<-c("PathToRObject")
NumberOfDatasets <- nrow(InputsTable)

seurat.object.list <- list()
for (dataset in rownames(InputsTable)) {
  print(dataset)
  DatasetIndexInInputsTable <- which(x = rownames(InputsTable) == dataset)
  InputRobject <- InputsTable[dataset,"PathToRObject"]
  seurat.object.list[[DatasetIndexInInputsTable]] <- readRDS(InputRobject)
}

################################################################################################################################################
################################################################################################################################################
### HERE ARE THE FUNCTIONS TO INTEGRATE DATASETS
################################################################################################################################################
################################################################################################################################################

writeLines("\n**** INTEGRATE DATASETS ****\n")

####################################
### Get reference datasets
####################################

writeLines("\n*** Get reference datasets ***\n")

if (regexpr("^NA$", DefaultParameters$ReferenceDatasets , ignore.case = T)[1] == 1) {
  writeLines("\n*** Will compare all-vs-all datasets to get anchors ***\n")
  ReferenceDatasets.indices <- c(1:nrow(InputsTable))
}else{
  writeLines("\n*** Determine reference dataset indices ***\n")
  ReferenceDatasets.list <- unlist(strsplit(DefaultParameters$ReferenceDatasets, ","))
  NumberOfFoundReferenceDatasetIDs <- sum(ReferenceDatasets.list %in% rownames(InputsTable) == T)
  if (NumberOfFoundReferenceDatasetIDs == length(ReferenceDatasets.list)) {
    ReferenceDatasets.indices <- match(ReferenceDatasets.list, rownames(InputsTable))
  }else{
    stop(paste0("Requested ", length(ReferenceDatasets.list), " datasets as references, but found ", NumberOfFoundReferenceDatasetIDs))
  }
  print(paste0("Will use datasets: ", paste(as.character(ReferenceDatasets.indices), sep = "", collapse = ",")))
}

####################################
### Get integration anchors
####################################
writeLines("\n*** Get integration anchors ***\n")

if (regexpr("^STACAS$", AnchorFinder, ignore.case = T)[1] == 1) {
  
  ####################################
  ### Get anchors with STACAS
  ####################################
  writeLines("\n*** Get anchors with STACAS ***\n")

  seurat.object.anchors.unfiltered <- FindAnchors.STACAS(object.list = seurat.object.list, dims=DefaultParameters$PcaDimsUse, anchor.features=DefaultParameters$StacasVarGenesIntegratedN, reference = ReferenceDatasets.indices, verbose = T)
  seurat.object.anchors <- FilterAnchors.STACAS(seurat.object.anchors.unfiltered, dist.thr = DefaultParameters$StacasDistPct, DefaultParameters$StacasDistThr)
  SampleTree <- SampleTree.STACAS(seurat.object.anchors)

} else if (regexpr("^Seurat$", AnchorFinder , ignore.case = T)[1] == 1) {
  
  ####################################
  ### Get anchors with Seurat
  ####################################
  writeLines("\n*** Get anchors with Seurat ***\n")
  
  seurat.object.integratedfeatures <- SelectIntegrationFeatures(object.list = seurat.object.list, nfeatures = DefaultParameters$IntegrationNFeatures)
  seurat.object.list <- PrepSCTIntegration(object.list = seurat.object.list, anchor.features = seurat.object.integratedfeatures, verbose = T)
  seurat.object.anchors <- FindIntegrationAnchors(object.list = seurat.object.list, k.filter = DefaultParameters$SeuratFindAnchorsKFilter, normalization.method = "SCT", dims=DefaultParameters$PcaDimsUse, anchor.features = seurat.object.integratedfeatures, reference = ReferenceDatasets.indices, verbose = T)
  SampleTree <- NULL
  
} else {
  stop("ERROR: unexpected anchors function")
}

####################################
### Integrating datasets
####################################
writeLines("\n*** Integrating datasets ***\n")

seurat.object.integrated <- IntegrateData(anchorset = seurat.object.anchors, normalization.method = "SCT", sample.tree = SampleTree, preserve.order = T, verbose = T)

################################################################################################################################################
################################################################################################################################################
### HERE ARE THE FUNCTIONS TO SAVE THE R_OBJECT AND LOG FILES
################################################################################################################################################
################################################################################################################################################

####################################
### Saving the Integration R object
####################################

writeLines("\n*** Saving the Integration R object ***\n")

OutfileRDS<-paste0(PathForOutfiles, "/", PrefixOutfiles, "_Integration.rds")
saveRDS(seurat.object.integrated, file = OutfileRDS)

####################################
### Obtain computing time used
####################################
writeLines("\n*** Obtain computing time used ***\n")

StopWatchEnd$Overall  <- Sys.time()

OutfileCPUtimes<-paste0(PathForOutfiles, "/", PrefixOutfiles, "_Integration_CPUtimes.txt")
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
