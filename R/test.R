# # # library(flowCore)
# # #
# # # source("/Volumes/TRANSCEND/GitHub/flowTOTAL/R/helperFunctions.R")
# # # source("/Volumes/TRANSCEND/GitHub/flowTOTAL/R/runPreprocessing.R")
# # # source("/Volumes/TRANSCEND/GitHub/flowTOTAL/R/runBackgating.R")
# # # source("/Volumes/TRANSCEND/GitHub/flowTOTAL/R/estimateProportion.R")
# # #
# # # fcs_path <- "/Volumes/ExtremeSSD/IBIMA/Citometria-IBIMA/fcs_raw-civil--11-07-2023"
# # # output <- "/Volumes/TRANSCEND/IBIMA/Citometria/flowTOTAL-pruebas"
# # #
# # # metadata <- data.frame(filename=list.files(fcs_path, full.names = T))
# # # metadata$file <-sapply(strsplit(metadata$filename, "/"), function(x) x[length(x)])
# # # metadata$filename_clean <- paste0(output, "/fcs_clean/", metadata$file)
# # #
# # # message("Step1: runPreprocessing()")
# # # # runPreprocessing(metadata, output = output, report = T, ncores = NULL)
# # #
# # # message("Step2: runDensityBackgating()")
# # # metadata <- metadata[grep("Panel1", metadata$filename), ]
# # # # runDensityBackgating(metadata = metadata, output = output, channel_bg = c("PE-Cy7-A+", "SSC-A-"))
# # #
# # # message("Step3: runDensityBackgating()")
# # # info_panel <- read.delim("/Volumes/TRANSCEND/IBIMA/Citometria/flowTOTAL-pruebas/Info_panel.txt")
# # # log_file_track <- "/Volumes/TRANSCEND/IBIMA/Citometria/flowTOTAL-pruebas/Log_file_track.txt"
# # # runEstimateProprotion(log_file_track, info_panel, output, ncores=NULL)
# # #
# # #
# # # remotes::install_github(repo = "ImmunoOncology/flowTOTAL", ref = "dev-1.0.0")
# # # library(flowTOTAL)
# # #
# # # flowTOTAL::runFlowTOTAL("/Volumes/TRANSCEND/IBIMA/Citometria/flowTOTAL-pruebas/fcs_raw",
# # #              "/Volumes/TRANSCEND/IBIMA/Citometria/flowTOTAL-pruebas",
# # #              panel_backgating=c("PE-Cy7-A+", "SSC-A-"),
# # #              panel_estimate="/Volumes/TRANSCEND/IBIMA/Citometria/flowTOTAL-pruebas/Info_panel.txt"
# # #             )
# # #
# #
# #
# #
# #
#
# remotes::install_github(repo = "ImmunoOncology/flowTOTAL", ref = "dev-1.0.0")
# library(flowCore)
# library(flowTOTAL)
# flowTOTAL::runFlowTOTAL("/Users/juanlu/Downloads/fcs_prueba",
#                         "/Users/juanlu/Downloads/fcs_prueba-res-master",
#                         panel_backgating=c("CD3+", "SSC-A-"),
#                         panel_estimate="/Users/juanlu/Downloads/fcs_prueba-res/Info_panel-Civil-Panel1.txt", cluster = NULL
# )
# #
# #
# # library(foreach)
# # library(parallel)
# #
# # cl <- parallel::makeCluster(ncores)
# # doParallel::registerDoParallel(cl)
# #
# #
# # flowTOTAL::runFlowTOTAL("/Users/juanlu/Downloads/fcs_prueba",
# #                         "/Users/juanlu/Downloads/fcs_prueba-res",
# #                         panel_backgating=c("CD3+", "SSC-A-"),
# #                         panel_estimate="Info_panel-Civil-Panel1.txt", cluster = NULL
# # )
# #
# # parallel::stopCluster(cl)
# #
# # flowTOTAL::runFlowTOTAL("fcs_raw_Civil_Panel1",
# #                         "flowTOTAL-Civil-Panel1",
# #                         panel_backgating=c("PE-Cy7-A+", "SSC-A-"),
# #                         panel_estimate="Info_panel-Civil-Panel1.txt", cluster = NULL
# # )
# #
# #
# # runEstimateProprotion(log_file_track = "pruebas-flowTOTAL/Log_file_track.txt", info_panel = read.delim("Info_panel-pruebas.txt"), output = "pruebas-flowTOTAL", ncores = 2)
# #
# #
# # runEstimateProprotion(log_file_track = paste0(file2, "/Log_file_track.txt"),
# #                       info_panel = read.delim(panel_estimate), output = file2, cluster=cl)
# #
