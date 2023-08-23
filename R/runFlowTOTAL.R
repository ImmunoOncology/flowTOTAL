library(flowCore)

source("/Volumes/TRANSCEND/GitHub/flowTOTAL/R/helperFunctions.R")
source("/Volumes/TRANSCEND/GitHub/flowTOTAL/R/runPreprocessing.R")
source("/Volumes/TRANSCEND/GitHub/flowTOTAL/R/runBackgating.R")
source("/Volumes/TRANSCEND/GitHub/flowTOTAL/R/estimateProportion.R")

fcs_path <- "/Volumes/ExtremeSSD/IBIMA/Citometria-IBIMA/fcs_raw-civil--11-07-2023"
output <- "/Volumes/TRANSCEND/IBIMA/Citometria/flowTOTAL-pruebas"

metadata <- data.frame(filename=list.files(fcs_path, full.names = T))
metadata$file <-sapply(strsplit(metadata$filename, "/"), function(x) x[length(x)]) 
metadata$filename_clean <- paste0(output, "/fcs_clean/", metadata$file)

message("Step1: runPreprocessing()")
# runPreprocessing(metadata, output = output, report = T, ncores = NULL)

message("Step2: runDensityBackgating()")
metadata <- metadata[grep("Panel1", metadata$filename), ]
# runDensityBackgating(metadata = metadata, output = output, channel_bg = c("PE-Cy7-A+", "SSC-A-"))

message("Step3: runDensityBackgating()")
info_panel <- read.delim("/Volumes/TRANSCEND/IBIMA/Citometria/flowTOTAL-pruebas/Info_panel.txt")
log_file_track <- "/Volumes/TRANSCEND/IBIMA/Citometria/flowTOTAL-pruebas/Log_file_track.txt"
runEstimateProprotion(log_file_track, info_panel, output, ncores=NULL)

  

  
runFlowTOTAL <- function(fcs_path, output, panel_backgating=NULL, panel_estimate=NULL){

  metadata <- data.frame(filename=list.files(fcs_path, full.names = T))
  metadata$file <-sapply(strsplit(metadata$filename, "/"), function(x) x[length(x)]) 
  metadata$filename_clean <- paste0(output, "/fcs_clean/", metadata$file)
  
  message("Step1: runPreprocessing()")
  runPreprocessing(metadata, output = output, report = T, ncores = NULL)
  
  message("Step2: runDensityBackgating()")
  runDensityBackgating(metadata = metadata, output = output, channel_bg = panel_backgating)
  
  message("Step3: runDensityBackgating()")
  info_panel <- read.delim(panel_estimate)
  log_file_track <- paste0(output, "/Log_file_track.txt")
  runEstimateProprotion(log_file_track, info_panel, output, ncores=NULL)
  
}
  
  