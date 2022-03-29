
setwd("/Volumes/TRANSCEND/IBIMA/flowTOTAL/flowTOTAL")
source("R/fn-helper.R")
source("R/preprocessing_function.R")
library(flowCore)

metadata <- read.delim("/Volumes/ExtremeSSD/IBIMA/Citometria/metadata.txt")
all_fcs <- list.files("/Volumes/ExtremeSSD/IBIMA/Citometria/fcs_clean", full.names = T, pattern = ".fcs$")
metadata$filename_clean <- gsub("/Volumes/ExtremeSSD/IBIMA/Citometria/fcs_raw/", "/Volumes/ExtremeSSD/IBIMA/Citometria/fcs_clean/", metadata$filename)
metadata <- metadata[match(all_fcs, metadata$filename_clean), ]

output <- "/Volumes/TRANSCEND/IBIMA/Citometria/pruebas/"
metadata_panel1 <- metadata[metadata$panel%in%"Panel1", ]

library(flowCore)
library(ggplot2)
library(openCyto)
library(ggcyto)


log_file <- "/Volumes/ExtremeSSD/IBIMA/Citometria/Log_file_error.txt"
log_file_error <- function(messages){
  sink(log_file, append = T)
  for(i in messages){
    cat(i, "\n")
  }
  sink()
}

log_file_traditional <- "/Volumes/ExtremeSSD/IBIMA/Citometria/Log_file_error_traditional.txt"
log_file_error_traditional <- function(messages){
  sink(log_file_traditional, append = T)
  for(i in messages){
    cat(i, "\n")
  }
  sink()
}

track_file <- "/Volumes/ExtremeSSD/IBIMA/Citometria/Log_file_track.txt"
log_file_track <- function(messages){
  sink(track_file, append = T)
  for(i in messages){
    cat(i, "\n")
  }
  sink()
}


# # Panel5 - Civil ----------------------------------------------------------
# panel5.civil <- metadata[metadata$Location=="Civil" & metadata$panel=="Panel3", ]
# panel5.civil$panel <- "Panel5"
# metadata <- rbind(metadata, panel5.civil)

channel <- list()
channel[["Panel1"]] <- c("CD3+", "SSC-A-")
channel[["Panel2"]] <- c("CD4+", "SSC-A-")
channel[["Panel3"]] <- c("CD16+", "SSC-A-")
channel[["Panel4"]] <- c("CD19+", "SSC-A-")
channel[["Panel5"]] <- c("CD14+", "SSC-A+")

channel_lg_list <- list()
channel_lg_list[["Panel1"]] <- c(1)
channel_lg_list[["Panel2"]] <- c(1)
channel_lg_list[["Panel3"]] <- c(1)
channel_lg_list[["Panel4"]] <- c(1)
channel_lg_list[["Panel5"]] <- c(1)

#metadata <- metadata[sample(1:nrow(metadata), 300), ]
#metadata <- metadata[metadata$Panel=="Panel1", ]
rownames(metadata) <- 1:nrow(metadata)

cl <- parallel::makeCluster(6, setup_strategy = "sequential")
doParallel::registerDoParallel(cl)
library(doParallel)
system.time(clean <- foreach(i = 1:nrow(metadata)) %dopar% {

  library(flowCore)
  library(ggcyto)
  library(ggplot2)
  library(openCyto)

  library(trip)
  library(sp)

  filename_clean <- metadata$filename_clean[i]
  panel <- metadata$Panel[i]
  filename <- unlist(strsplit(filename_clean, "/"))
  filename <- filename[length(filename)]

  output.dir <- paste0("/Volumes/TRANSCEND/IBIMA/Citometria/Paper_results/", panel)
  if(!dir.exists(output.dir)) dir.create(output.dir)
  if(!dir.exists(paste0(output.dir, "/PDF"))) dir.create(paste0(output.dir, "/PDF"))
  if(!dir.exists(paste0(output.dir, "/FCS"))) dir.create(paste0(output.dir, "/FCS"))

  tryCatch({

    ff <- flowCore::read.FCS(filename_clean)

    location <- metadata$Location[i]
    my_channel <- channel[[panel]]
    channel_sign <- sapply(my_channel, function(x) substr(x, nchar(x), nchar(x)))
    my_channel <-  sapply(my_channel, function(x) substr(x, 1, nchar(x)-1))
    ff@parameters@data$desc[is.na(ff@parameters@data$desc)] <- ff@parameters@data$name[is.na(ff@parameters@data$desc)]
    channel_bg <- as.vector(ff@parameters@data$name[match(my_channel, ff@parameters@data$desc)])
    channel_bg <- paste0(channel_bg,  channel_sign)

    cutpoint_min <- 1000 # 5000
    chnl = c("FSC-A", "SSC-A")
    cutpoint_max = 50000
    logicle_chnls <- channel_lg_list[[panel]]
    res_bg <- do_backgating(ff, filename = filename, channel_bg = channel_bg, logicle_chnls = logicle_chnls, output.dir = output.dir, cutpoint_min = cutpoint_min)
    log_file_track(paste(filename_clean, filename, paste0(output.dir, "/FCS/", filename), panel, nrow(ff@exprs), res_bg, "Completed", i, sep = "\t"))

        },
  error=function(e) {
    log_file_error(paste(e, "Iter-->", i, "\n", "File: ", metadata$filename[i]))
    log_file_track(paste(filename_clean, filename, "\t", NA, NA, "Error", i, sep = "\t"))
    return(e)
  })

})

parallel::stopCluster(cl)



# Error -------------------------------------------------------------------

error <- read.delim("/Volumes/ExtremeSSD/IBIMA/Citometria/Log_file_track.txt", header = F)
error <- error[error$V7=="Error", ]

error_list <- list()
error_plot <- list()
for(i in error$V1){
  id <- unlist(strsplit(i, "/"))[7]
  error_list[[id]] <- as.vector(na.omit(read.FCS(i)@parameters@data$desc))
  error_plot[[id]] <- ggcyto::autoplot(read.FCS(i))
}

# metadata$filename[c(100, 99, 122, 134, 957)]
# ggcyto::autoplot(read.FCS(metadata$filename[100]), "CD19", "SSC-A", 100)
# ggcyto::autoplot(read.FCS(metadata$filename[99]), "CD19", "SSC-A", 100)
# ggcyto::autoplot(read.FCS(metadata$filename[122]), "CD3", "SSC-A", 100)
# ggcyto::autoplot(read.FCS(metadata$filename[134]), "CD3", "SSC-A", 100)
# ggcyto::autoplot(read.FCS(metadata$filename[957]), "CD14", "SSC-A", 100)
#
#
# table(unlist(lapply(lapply(raw, read.FCS), function(a) a@parameters@data$desc)))
#
# panel3_civil_to_panel5_civil <- raw[grepl("Panel3", raw) & grepl("Civil", raw)]
# panel3_civil_to_panel5_civil_new <- gsub("_Panel3_", "_Panel5_", panel3_civil_to_panel5_civil)
# file.copy(panel3_civil_to_panel5_civil, panel3_civil_to_panel5_civil_new)
#
# panel4_marbella_to_panel5_marbella <- raw[grepl("Panel4", raw) & grepl("Marbella", raw)]
# panel4_marbella_to_panel5_marbella_new <- gsub("_Panel4_", "_Panel5_", panel4_marbella_to_panel5_marbella)
# file.copy(panel4_marbella_to_panel5_marbella, panel4_marbella_to_panel5_marbella_new)
# file.remove(panel4_marbella_to_panel5_marbella)
#
# panel3_marbella_to_panel4_marbella <- raw[grepl("Panel3", raw) & grepl("Marbella", raw)]
# panel3_marbella_to_panel4_marbella_new <- gsub("_Panel3_", "_Panel4_", panel3_marbella_to_panel4_marbella)
# file.copy(panel3_marbella_to_panel4_marbella, panel3_marbella_to_panel4_marbella_new)


# Tradicional ------------------------------------------------------------

info_panel <- read.delim("/Volumes/ExtremeSSD/IBIMA/Citometria/Info_panel.txt")
output.dir <- "/Volumes/TRANSCEND/IBIMA/Citometria/IMK_results/results"
if(!dir.exists(output.dir)) dir.create(output.dir)

files <- list.files("/Volumes/TRANSCEND/IBIMA/Citometria/IMK_results/", full.names = T, pattern = ".fcs$", recursive = T)
files <- files[grep("FCS_raw", files)]

batch <- rep("Marbella", length(files))
batch[grepl("Civil",files)] <- "Civil"

panel <- rep("Panel1", length(files))
panel[grepl("Panel2",files)] <- "Panel2"
panel[grepl("Panel3",files)] <- "Panel3"
panel[grepl("Panel4",files)] <- "Panel4"
panel[grepl("Panel5",files)] <- "Panel5"


cl <- parallel::makeCluster(6, setup_strategy = "sequential")
doParallel::registerDoParallel(cl)
library(doParallel)
system.time(clean <- foreach(i = 1:length(files)) %dopar% {
  library(flowCore)
  library(ggcyto)
  library(ggplot2)
  library(openCyto)

  tryCatch({

    id <- gsub(".fcs$", "", unlist(strsplit(files[i], "/"))[[length(unlist(strsplit(files[i], "/")))]])
    do_traditional(filename = files[i], id = id, panel = panel[i],  batch = batch[i], info_panel = info_panel, output.dir = output.dir, cutpoint_min = 0, cutpoint_max = 30000)

  },
  error=function(e) {
    log_file_error_traditional(paste(e, "Iter-->", i, "\n", "File: ", files[i]))
    return(e)
  })

})


counts_panel1 <- reshape2::melt(read.delim("/Volumes/TRANSCEND/IBIMA/Citometria/IMK_results/results/Traditional_counts_Panel1.txt"))
counts_panel2 <- reshape2::melt(read.delim("/Volumes/TRANSCEND/IBIMA/Citometria/IMK_results/results/Traditional_counts_Panel2.txt"))
counts_panel3 <- reshape2::melt(read.delim("/Volumes/TRANSCEND/IBIMA/Citometria/IMK_results/results/Traditional_counts_Panel3.txt"))
counts_panel4 <- reshape2::melt(read.delim("/Volumes/TRANSCEND/IBIMA/Citometria/IMK_results/results/Traditional_counts_Panel4.txt"))
counts_panel5 <- reshape2::melt(read.delim("/Volumes/TRANSCEND/IBIMA/Citometria/IMK_results/results/Traditional_counts_Panel5.txt"))
counts_panel3 <- counts_panel3[counts_panel3$ID!="", ]

counts_panel1$variable <- as.character(counts_panel1$variable)
counts_panel1$variable[counts_panel1$variable=="N.total"] <- "N.total.panel1"
counts_panel2$variable <- as.character(counts_panel2$variable)
counts_panel2$variable[counts_panel2$variable=="N.total"] <- "N.total.panel2"
counts_panel3$variable <- as.character(counts_panel3$variable)
counts_panel3$variable[counts_panel3$variable=="N.total"] <- "N.total.panel3"
counts_panel4$variable <- as.character(counts_panel4$variable)
counts_panel4$variable[counts_panel4$variable=="N.total"] <- "N.total.panel4"
counts_panel5$variable <- as.character(counts_panel5$variable)
counts_panel5$variable[counts_panel5$variable=="N.total"] <- "N.total.panel5"

counts <- rbind(counts_panel1, counts_panel2, counts_panel3, counts_panel4, counts_panel5)
counts <- counts[, c("ID", "variable", "value")]
counts$ID <- unlist(lapply(lapply(strsplit(counts$ID, "_"), `[`, c(3,4)), paste, sep = "", collapse = "-"))

counts$ID[duplicated(paste0(counts$ID, counts$variable))]

metadata <- read.delim("/Volumes/ExtremeSSD/IBIMA/Citometria/metadata_complete.txt")
#metadata <- unique(metadata[, c("Response_3_months", "Response_6_months", "Time", "IMK", "Diagnosis", "Location", "Diagnosis")])
metadata <- unique(metadata[, -which(colnames(metadata)%in%c("Panel", "ID", "filename", "TP_is_valid"))])
metadata$ID <- paste( metadata$IMK, metadata$Time, sep = "-")

counts$ID[!counts$ID%in%metadata$ID]

counts <- cbind(counts, metadata[match(counts$ID, metadata$ID), ])
counts <- reshape2::dcast(counts, ID~variable, value.var = "value")
counts <- cbind(counts, metadata[match(counts$ID, metadata$ID), ])

## Checks NA
View(counts[rowSums(is.na(counts[, -6]))>0, ])

# Panel1
ids <- na.omit(counts$ID[is.na(counts[, c(3, 5, 7, 13, 14)])])
ids <- paste0("Panel1_", gsub("-", "_", ids))
ids <- unlist(sapply(ids, grep, files))
table(unlist(lapply(ids, function(x){
  ff <- read.FCS(files[x])
  return(unlist(ff@parameters@data$desc))
})))

# Panel2
ids <- na.omit(counts$ID[is.na(counts[, c(4, 8, 15)])])
ids <- paste0("Panel2_", gsub("-", "_", ids))
ids <- unlist(sapply(ids, grep, files))


# Panel3
ids <- na.omit(counts$ID[is.na(counts[, c(9, 12)])])

# Panel4
ids <- na.omit(counts$ID[is.na(counts[, c(2, 10)])])
ids <- paste0("Panel4_", gsub("-", "_", ids))
ids <- unlist(sapply(ids, grep, files))

# Panel5
ids <- na.omit(counts$ID[is.na(counts[, c(6, 11)])])
ids <- paste0("Panel5_", gsub("-", "_", ids))
ids <- unlist(sapply(ids, grep, files))


error <- read.delim("/Volumes/ExtremeSSD/IBIMA/Citometria/Log_file_track.txt", header = F)
error$ID <- unlist(lapply(lapply(strsplit(error$V2, "_"), `[`, c(3, 4)), paste0, collapse="-"))

error.p1 <- error[grep("Panel1", error$V2), ]
error.p2 <- error[grep("Panel2", error$V2), ]
error.p3 <- error[grep("Panel3", error$V2), ]
error.p4 <- error[grep("Panel4", error$V2), ]
error.p5 <- error[grep("Panel5", error$V2), ]

counts$N.event.panel1 <- error.p1$V5[match(counts$ID, error.p1$ID)]
counts$N.event.panel2 <- error.p2$V5[match(counts$ID, error.p2$ID)]
counts$N.event.panel3 <- error.p3$V5[match(counts$ID, error.p3$ID)]
counts$N.event.panel4 <- error.p4$V5[match(counts$ID, error.p4$ID)]
counts$N.event.panel5 <- error.p5$V5[match(counts$ID, error.p5$ID)]

write.table(counts, "/Volumes/TRANSCEND/IBIMA/Citometria/IMK_results/results/IMK-Counts.txt", quote = F, sep = "\t", col.names = T, row.names = F)


# metadata <- read.delim("/Volumes/ExtremeSSD/IBIMA/Citometria/metadata_complete.txt")
# metadata$clean <- gsub("fcs_raw", "fcs_clean", metadata$filename)
# track <- read.delim("/Volumes/ExtremeSSD/IBIMA/Citometria/Log_file_track.txt", header = F)
# track <- cbind(track, metadata[match(track$V1, metadata$clean), ])
# track$ID <- paste(track$Location, track$panel, track$IMK, track$Time, track$Time_date, sep = "_")
#
# counts <- read.delim("/Volumes/TRANSCEND/IBIMA/Citometria/IMK_results/results/Traditional_counts.txt")
#
# counts[grepl("Panel5", counts$File), "ID"] <- gsub("_Panel3_", "_Panel5_", counts[grepl("Panel5", counts$File), "ID"])
# track[grepl("Panel5", track$V3), "ID"] <- gsub("_Panel3_", "_Panel5_", track[grepl("Panel5", track$V3), "ID"])
#
# all(counts$ID%in%track$ID)
#
# counts <- cbind(track, counts[match(track$ID, counts$ID), -c(1, 2)])


# simplyfi ----------------------------------------------------------------


track <- read.delim("/Volumes/ExtremeSSD/IBIMA/Citometria/Log_file_track.txt", header = F)
track <- track[track$V7=="Completed", ]

track$batch <- "Marbella"
track$batch[grepl("Civil", track$V3)] <- "Civil"

keep <- c("Time", "SSC-H", "SSC-A", "Original_ID", "FSC-A", "FSC-H")
info_experiment <- read.delim("/Volumes/ExtremeSSD/IBIMA/Citometria/Info_experiment.txt")

check <- list()

track$V3[!file.exists(track$V3)]

output.dir <- paste0("/Volumes/TRANSCEND/IBIMA/Citometria/IMK_results/")
output.panel <- list.dirs(output.dir, recursive = F)
output.panel <- output.panel[grep("Panel", output.panel)]
lapply(output.panel, function(file_panel){
  R.utils::copyDirectory(file.path(file_panel, "FCS"), file.path(file_panel, "FCS_raw"))
})


for(i in unique(track$batch)){
  check[[i]] <- list()
  for(j in unique(track$V4)){
    check[[i]][[j]] <- list()
    files <- track$V3[track$batch==i&track$V4==j]
    channel <- info_experiment$Channels[info_experiment$Batch==i&info_experiment$Panel==j]
    channel <- unlist(strsplit(channel, "[.]"))
    check[[i]][[j]] <- sapply(files, simplify_flowCore, keep = c(keep, channel))
  }
  message("Track ", i)
  message("Panel ", j)
}


# SCE ---------------------------------------------------------------------

if(!dir.exists("/Volumes/TRANSCEND/IBIMA/Citometria/IMK_results/RDS"))
  dir.create("/Volumes/TRANSCEND/IBIMA/Citometria/IMK_results/RDS")

metadata <- read.delim("/Volumes/ExtremeSSD/IBIMA/Citometria/metadata_complete.txt")
metadata$clean <- gsub("fcs_raw", "fcs_clean", metadata$filename)
track <- read.delim("/Volumes/ExtremeSSD/IBIMA/Citometria/Log_file_track.txt", header = F)
track <- track[track$V7=="Completed", ]
track <- track[track$V1%in%metadata$clean, ]
metadata <- metadata[match(track$V1, metadata$clean), ]
metadata$FCS_final <- track[match(metadata$clean, track$V1), "V3"]
all(file.exists(metadata$FCS_final))
metadata$FCS_final[!file.exists(metadata$FCS_final)]

# Panel1 ------------------------------------------------------------------

panel1 <- metadata$FCS_final[metadata$Panel=="Panel1"]
metadata_panel1 <-  metadata[metadata$Panel=="Panel1", ]
metadata_panel1 <- metadata_panel1[!metadata_panel1$FCS_final%in%c(
  "/Volumes/TRANSCEND/IBIMA/Citometria/IMK_results/Panel1/FCS/Marbella_Panel1_062IMK_T2_2019-06-12.fcs",
  "/Volumes/TRANSCEND/IBIMA/Citometria/IMK_results/Panel1/FCS/Marbella_Panel1_064IMK_T1_2019-06-25.fcs",
  "/Volumes/TRANSCEND/IBIMA/Citometria/IMK_results/Panel1/FCS/Marbella_Panel1_065IMK_T1_2019-06-20.fcs",
  "/Volumes/TRANSCEND/IBIMA/Citometria/IMK_results/Panel1/FCS/Marbella_Panel1_067IMK_T1_2019-06-26.fcs"
), ]

rownames(metadata_panel1) <- metadata_panel1$FCS_final
panel1 <- metadata_panel1$FCS_final

fcs <- CytoTree::runExprsExtract(panel1[1])
exclusions <- c("Original_ID<Original_ID>", "Time")
exclusions <- c("Original_ID<Original_ID>", "Time", "FSC-A<FSC-A>" ,"FSC-H<FSC-H>", "SSC-A<SSC-A>", "SSC-H<SSC-H>")
inclusions <- colnames(fcs)[!colnames(fcs)%in%exclusions]
inclusions <- unlist(lapply(strsplit(inclusions, "<"), `[`, 1))
markernames <- colnames(fcs)[!colnames(fcs)%in%exclusions]
markernames <- gsub(">", "", markernames)
markernames <- unlist(lapply(strsplit(markernames, "<"), `[`, 2))
names(markernames) <- inclusions

sce = scDataviz::processFCS(
  files = panel1,
  metadata = metadata_panel1,
  transformation = TRUE, # Transform data
  transFun = function (x) asinh(x), # Using asinh function (cofactor 5 is used by default)
  downsampleVar = NULL, # No downsampling by variance
  downsample = NULL, # Downsample to n cells per downsample_grouping
  colsRetain = inclusions,
  asinhFactor = 150,
  filter = TRUE,
  bgNoiseThreshold = -Inf,
  colsDiscard = exclusions,
  newColnames = markernames) # Discard columns not selected for downstream analysis
saveRDS(sce, "/Volumes/TRANSCEND/IBIMA/Citometria/IMK_results/RDS/panel1.RDS")

idt <- metadata_panel1$Time=="T1"
sce_T1 = scDataviz::processFCS(
  files = panel1[idt],
  metadata = metadata_panel1[idt, ],
  transformation = TRUE, # Transform data
  transFun = function (x) asinh(x), # Using asinh function (cofactor 5 is used by default)
  downsampleVar = NULL, # No downsampling by variance
  downsample = NULL, # Downsample to n cells per downsample_grouping
  colsRetain = inclusions,
  asinhFactor = 150,
  bgNoiseThreshold = -Inf,
  colsDiscard = exclusions,
  newColnames = markernames) # Discard columns not selected for downstream analysis
saveRDS(sce_T1, "/Volumes/TRANSCEND/IBIMA/Citometria/IMK_results/RDS/panel1_T1.RDS")

idt <- metadata_panel1$Diagnosis=="Lung"
sce_lung = scDataviz::processFCS(
  files = panel1[idt],
  metadata = metadata_panel1[idt, ],
  transformation = TRUE, # Transform data
  transFun = function (x) asinh(x), # Using asinh function (cofactor 5 is used by default)
  downsampleVar = NULL, # No downsampling by variance
  downsample = NULL, # Downsample to n cells per downsample_grouping
  colsRetain = inclusions,
  asinhFactor = 150,
  bgNoiseThreshold = -Inf,
  colsDiscard = exclusions,
  newColnames = markernames) # Discard columns not selected for downstream analysis
saveRDS(sce_lung, "/Volumes/TRANSCEND/IBIMA/Citometria/IMK_results/RDS/panel1_Lung.RDS")

idt <- metadata_panel1$Diagnosis=="Melanoma"
sce_melanoma = scDataviz::processFCS(
  files = panel1[idt],
  metadata = metadata_panel1[idt, ],
  transformation = TRUE, # Transform data
  transFun = function (x) asinh(x), # Using asinh function (cofactor 5 is used by default)
  downsampleVar = NULL, # No downsampling by variance
  downsample = NULL, # Downsample to n cells per downsample_grouping
  colsRetain = inclusions,
  asinhFactor = 150,
  bgNoiseThreshold = -Inf,
  colsDiscard = exclusions,
  newColnames = markernames) # Discard columns not selected for downstream analysis
saveRDS(sce_melanoma, "/Volumes/TRANSCEND/IBIMA/Citometria/IMK_results/RDS/panel1_Melanoma.RDS")


idt_Time_T1 <- metadata_panel1$Time=="T1"

idt <- metadata_panel1$Diagnosis=="Lung" & idt_Time_T1
sce_lung = scDataviz::processFCS(
  files = panel1[idt],
  metadata = metadata_panel1[idt, ],
  transformation = TRUE, # Transform data
  transFun = function (x) asinh(x), # Using asinh function (cofactor 5 is used by default)
  downsampleVar = NULL, # No downsampling by variance
  downsample = NULL, # Downsample to n cells per downsample_grouping
  colsRetain = inclusions,
  asinhFactor = 150,
  bgNoiseThreshold = -Inf,
  colsDiscard = exclusions,
  newColnames = markernames) # Discard columns not selected for downstream analysis
saveRDS(sce_lung, "/Volumes/TRANSCEND/IBIMA/Citometria/IMK_results/RDS/panel1_Lung_T1.RDS")

idt <- metadata_panel1$Diagnosis=="Melanoma"& idt_Time_T1
sce_melanoma = scDataviz::processFCS(
  files = panel1[idt],
  metadata = metadata_panel1[idt, ],
  transformation = TRUE, # Transform data
  transFun = function (x) asinh(x), # Using asinh function (cofactor 5 is used by default)
  downsampleVar = NULL, # No downsampling by variance
  downsample = NULL, # Downsample to n cells per downsample_grouping
  colsRetain = inclusions,
  asinhFactor = 150,
  bgNoiseThreshold = -Inf,
  colsDiscard = exclusions,
  newColnames = markernames) # Discard columns not selected for downstream analysis
saveRDS(sce_melanoma, "/Volumes/TRANSCEND/IBIMA/Citometria/IMK_results/RDS/panel1_Melanoma_T1.RDS")

idt <- metadata_panel1$Time=="T2"
sce_T2 = scDataviz::processFCS(
  files = panel1[idt],
  metadata = metadata_panel1[idt, ],
  transformation = TRUE, # Transform data
  transFun = function (x) asinh(x), # Using asinh function (cofactor 5 is used by default)
  downsampleVar = NULL, # No downsampling by variance
  downsample = NULL, # Downsample to n cells per downsample_grouping
  colsRetain = inclusions,
  bgNoiseThreshold = -Inf,
  asinhFactor = 150,
  colsDiscard = exclusions,
  newColnames = markernames) # Discard columns not selected for downstream analysis
saveRDS(sce_T2, "/Volumes/TRANSCEND/IBIMA/Citometria/IMK_results/RDS/panel1_T2.RDS")


idt_Time_T2 <- metadata_panel1$Time=="T2"

idt <- metadata_panel1$Diagnosis=="Lung" & idt_Time_T2
sce_lung = scDataviz::processFCS(
  files = panel1[idt],
  metadata = metadata_panel1[idt, ],
  transformation = TRUE, # Transform data
  transFun = function (x) asinh(x), # Using asinh function (cofactor 5 is used by default)
  downsampleVar = NULL, # No downsampling by variance
  downsample = NULL, # Downsample to n cells per downsample_grouping
  colsRetain = inclusions,
  asinhFactor = 150,
  bgNoiseThreshold = -Inf,
  colsDiscard = exclusions,
  newColnames = markernames) # Discard columns not selected for downstream analysis
saveRDS(sce_lung, "/Volumes/TRANSCEND/IBIMA/Citometria/IMK_results/RDS/panel1_Lung_T2.RDS")

idt <- metadata_panel1$Diagnosis=="Melanoma" & idt_Time_T2
sce_melanoma = scDataviz::processFCS(
  files = panel1[idt],
  metadata = metadata_panel1[idt, ],
  transformation = TRUE, # Transform data
  transFun = function (x) asinh(x), # Using asinh function (cofactor 5 is used by default)
  downsampleVar = NULL, # No downsampling by variance
  downsample = NULL, # Downsample to n cells per downsample_grouping
  colsRetain = inclusions,
  asinhFactor = 150,
  bgNoiseThreshold = -Inf,
  colsDiscard = exclusions,
  newColnames = markernames) # Discard columns not selected for downstream analysis
saveRDS(sce_melanoma, "/Volumes/TRANSCEND/IBIMA/Citometria/IMK_results/RDS/panel1_Melanoma_T2.RDS")


# Panel2 ------------------------------------------------------------------

panel2 <- metadata$FCS_final[metadata$Panel=="Panel2"]
metadata_panel2 <-  metadata[metadata$Panel=="Panel2", ]
rownames(metadata_panel2) <- metadata_panel2$FCS_final

fcs <- CytoTree::runExprsExtract(panel2[1])
exclusions <- c("Original_ID<Original_ID>", "Time")
inclusions <- colnames(fcs)[!colnames(fcs)%in%exclusions]
inclusions <- unlist(lapply(strsplit(inclusions, "<"), `[`, 1))
markernames <- colnames(fcs)[!colnames(fcs)%in%exclusions]
markernames <- gsub(">", "", markernames)
markernames <- unlist(lapply(strsplit(markernames, "<"), `[`, 2))
names(markernames) <- inclusions

sce = scDataviz::processFCS(
  files = panel2,
  metadata = metadata_panel2,
  transformation = TRUE, # Transform data
  transFun = function (x) asinh(x), # Using asinh function (cofactor 5 is used by default)
  downsampleVar = NULL, # No downsampling by variance
  downsample = NULL, # Downsample to n cells per downsample_grouping
  colsRetain = inclusions,
  asinhFactor = 150,
  colsDiscard = exclusions,
  newColnames = markernames) # Discard columns not selected for downstream analysis
saveRDS(sce, "/Volumes/TRANSCEND/IBIMA/Citometria/IMK_results/RDS/panel2.RDS")

idt <- metadata_panel2$Time=="T1"
sce_T1 = scDataviz::processFCS(
  files = panel2[idt],
  metadata = metadata_panel2[idt, ],
  transformation = TRUE, # Transform data
  transFun = function (x) asinh(x), # Using asinh function (cofactor 5 is used by default)
  downsampleVar = NULL, # No downsampling by variance
  downsample = NULL, # Downsample to n cells per downsample_grouping
  colsRetain = inclusions,
  asinhFactor = 150,
  colsDiscard = exclusions,
  newColnames = markernames) # Discard columns not selected for downstream analysis
saveRDS(sce_T1, "/Volumes/TRANSCEND/IBIMA/Citometria/IMK_results/RDS/panel2_T1.RDS")

idt <- metadata_panel2$Diagnosis=="Lung"
sce_lung = scDataviz::processFCS(
  files = panel2[idt],
  metadata = metadata_panel2[idt, ],
  transformation = TRUE, # Transform data
  transFun = function (x) asinh(x), # Using asinh function (cofactor 5 is used by default)
  downsampleVar = NULL, # No downsampling by variance
  downsample = NULL, # Downsample to n cells per downsample_grouping
  colsRetain = inclusions,
  asinhFactor = 150,
  colsDiscard = exclusions,
  newColnames = markernames) # Discard columns not selected for downstream analysis
saveRDS(sce_lung, "/Volumes/TRANSCEND/IBIMA/Citometria/IMK_results/RDS/panel2_Lung.RDS")

idt <- metadata_panel2$Diagnosis=="Melanoma"
sce_melanoma = scDataviz::processFCS(
  files = panel2[idt],
  metadata = metadata_panel2[idt, ],
  transformation = TRUE, # Transform data
  transFun = function (x) asinh(x), # Using asinh function (cofactor 5 is used by default)
  downsampleVar = NULL, # No downsampling by variance
  downsample = NULL, # Downsample to n cells per downsample_grouping
  colsRetain = inclusions,
  asinhFactor = 150,
  colsDiscard = exclusions,
  newColnames = markernames) # Discard columns not selected for downstream analysis
saveRDS(sce_melanoma, "/Volumes/TRANSCEND/IBIMA/Citometria/IMK_results/RDS/panel2_Melanoma.RDS")


# Panel4 ------------------------------------------------------------------

Panel4 <- metadata$FCS_final[metadata$Panel=="Panel4"& metadata$Location=="Civil"]
metadata_Panel4 <-  metadata[metadata$Panel=="Panel4" & metadata$Location=="Civil", ]

Panel4 <- gsub("/Volumes/TRANSCEND/IBIMA/Citometria/IMK_results/Panel4/FCS/", "", Panel4)
Panel4 <- paste0("/Volumes/TRANSCEND/IBIMA/Citometria/IMK_results/Panel4/FCS_raw/", Panel4)

metadata_Panel4$FCS_final <- gsub("/Volumes/TRANSCEND/IBIMA/Citometria/IMK_results/Panel4/FCS/", "", metadata_Panel4$FCS_final)
metadata_Panel4$FCS_final <- paste0("/Volumes/TRANSCEND/IBIMA/Citometria/IMK_results/Panel4/FCS_raw/", metadata_Panel4$FCS_final)

rownames(metadata_Panel4) <- metadata_Panel4$FCS_final

fcs <- CytoTree::runExprsExtract(Panel4[1])
exclusions <- c("Original_ID<Original_ID>", "Time")
exclusions <- c("Original_ID<Original_ID>", "Time", "FSC-A<FSC-A>" ,"FSC-H<FSC-H>", "SSC-A<SSC-A>", "SSC-H<SSC-H>", "FSC-W<FSC-W>", "SSC-W<SSC-W>")
inclusions <- colnames(fcs)[!colnames(fcs)%in%exclusions]
inclusions <- unlist(lapply(strsplit(inclusions, "<"), `[`, 1))
markernames <- colnames(fcs)[!colnames(fcs)%in%exclusions]
markernames <- gsub(">", "", markernames)
markernames <- unlist(lapply(strsplit(markernames, "<"), `[`, 2))
names(markernames) <- inclusions

sce = scDataviz::processFCS(
  files = Panel4,
  metadata = metadata_Panel4,
  transformation = TRUE, # Transform data
  transFun = function (x) asinh(x), # Using asinh function (cofactor 5 is used by default)
  downsampleVar = NULL, # No downsampling by variance
  downsample = NULL, # Downsample to n cells per downsample_grouping
  colsRetain = inclusions,
  asinhFactor = 150,
  bgNoiseThreshold = -Inf,
  colsDiscard = exclusions,
  newColnames = markernames) # Discard columns not selected for downstream analysis
saveRDS(sce, "/Volumes/TRANSCEND/IBIMA/Citometria/IMK_results/RDS/panel4.RDS")


