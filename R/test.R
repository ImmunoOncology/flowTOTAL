
setwd("/Volumes/TRANSCEND/IBIMA/flowTOTAL/flowTOTAL")
source("R/fn-helper.R")
source("R/preprocessing_function.R")

all_fcs <- list.files("/Volumes/TRANSCEND/IBIMA/Citometria/fcs_clean/", full.names = T)

output <- "/Volumes/TRANSCEND/IBIMA/Citometria/pruebas"
filename <- "mi_fcs.fcs"
file <- all_fcs[1]
run_Preprocessing(file = all_fcs[1], filename = filename, output = output)

ff <- flowCore::read.FCS(paste0(output, "/", filename))

