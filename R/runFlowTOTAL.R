

runFlowTOTAL <- function(fcs_path, output, panel_backgating=NULL, panel_estimate=NULL){

  metadata <- data.frame(filename=list.files(fcs_path, full.names = T))
  metadata$file <-sapply(strsplit(metadata$filename, "/"), function(x) x[length(x)]) 
  metadata$filename_clean <- paste0(output, "/fcs_clean/", metadata$file)
  
  message("Step1: runPreprocessing()")
  runPreprocessing(metadata=metadata, output = output, report = T, ncores = NULL)
  
  message("Step2: runDensityBackgating()")
  runDensityBackgating(metadata = metadata, output = output, channel_bg = panel_backgating)
  
  message("Step3: runDensityBackgating()")
  info_panel <- read.delim(panel_estimate)
  log_file_track <- paste0(output, "/Log_file_track.txt")
  runEstimateProprotion(log_file_track, info_panel, output, ncores=NULL)
  
}
  

  