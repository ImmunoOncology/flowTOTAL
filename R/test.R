
setwd("/Volumes/TRANSCEND/IBIMA/flowTOTAL/flowTOTAL")
source("R/fn-helper.R")
source("R/preprocessing_function.R")
library(flowCore)

metadata <- read.delim("/Volumes/ExtremeSSD/IBIMA/Citometria/metadata.txt")
all_fcs <- list.files("/Users/juanlu/Desktop/fcs_clean", full.names = T, pattern = ".fcs$")
metadata$filename_clean <- gsub("/Volumes/ExtremeSSD/IBIMA/Citometria/fcs_raw/", "/Users/juanlu/Desktop/fcs_clean/", metadata$filename)
metadata <- metadata[match(all_fcs, metadata$filename_clean), ]

output <- "/Volumes/TRANSCEND/IBIMA/Citometria/pruebas/"
metadata_panel1 <- metadata[metadata$panel%in%"Panel1", ]

set.seed(123)
idt <- sample(1:nrow(metadata_panel1), 50)
metadata_panel1 <- metadata_panel1[idt, ]

rownames(metadata_panel1) <- 1:nrow(metadata_panel1)
for(i in 1:nrow(metadata_panel1)){
  filename_clean <- metadata_panel1$filename_clean[i]
  ff <- flowCore::read.FCS(filename_clean)
  marker_bg <- c("CD3+")
  dir_plot <- paste0(output, gsub("/Users/juanlu/Desktop/fcs_clean/", "", filename_clean))
  dir_plot <- gsub(".fcs$", ".pdf", dir_plot)
  dummy <- do_backgating(ff, marker_bg = marker_bg, dir_plot=dir_plot)
}





dummy$plot
dummy$plot_back
dummy$plot_filter




