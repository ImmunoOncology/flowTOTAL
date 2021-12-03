
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

do_bg_pca <- function(ff, filename, output.dir, chnl = c("FSC-A", "SSC-A"), channel_bg, min.count = 5000){

  bg_filter <- openCyto::gate_flowclust_1d(fr=ff, channel_bg, K=2)
  if(bg_filter@min<min.count) bg_filter@min <- min.count
  chnl_filter <- openCyto::gate_flowclust_1d(fr=ff, chnl[2], K=2)
  idt <- ff@exprs[, channel_bg]>bg_filter@min & ff@exprs[, chnl[2]]<chnl_filter@min

  peaks_chnl1 <- openCyto:::.find_peaks(ff@exprs[idt, chnl[1]], num_peaks = 1)
  peaks_chnl2 <- openCyto:::.find_peaks(ff@exprs[idt, chnl[2]], num_peaks = 1)

  k <- nrow(openCyto:::.find_peaks(ff@exprs[idt, chnl[1]]))*nrow(openCyto:::.find_peaks(ff@exprs[idt, chnl[2]]))

  filter <- openCyto::gate_flowclust_2d(fr=ff[idt, ], quantile = 0.99, xChannel = chnl[1], yChannel = chnl[2], K=6,target = c(peaks_chnl1$x[1], peaks_chnl2$x[1]))
  p1 <- ggcyto::autoplot(ff, x = chnl[2], y = channel_bg, bins = 100)+geom_gate(bg_filter)+geom_gate(chnl_filter)
  p2 <- ggcyto::autoplot(ff[idt, ], x = chnl[1], y = chnl[2], bins = 100)+geom_gate(filter)+
    xlim(range(ff@exprs[, chnl[1]]))+ylim(range(ff@exprs[, chnl[2]]))
  p3 <- ggcyto::autoplot(ff, x = chnl[1], y = chnl[2], bins = 100)+geom_gate(filter)

  p <- ggpubr::ggarrange(plotlist = list(as.ggplot(p1), as.ggplot(p2), as.ggplot(p3)))
  ggsave(paste0(output.dir, "/", filename, ".pdf"), plot = p, device = "pdf", width = 18, height = 12)
}

do_double_bg_pca <- function(ff, filename, output.dir, chnl = c("FSC-A", "SSC-A"), channel_bg, min.count = 10000){

  #bg_filter <- openCyto::gate_flowclust_1d(fr=ff, channel_bg, K=2)
  #if(bg_filter@min<min.count) bg_filter@min <- min.count
  chnl_filter <- openCyto::gate_flowclust_1d(fr=ff, chnl[2], K=2)
  idt <- ff@exprs[, chnl[2]]<chnl_filter@min
  bg_filter <- openCyto::gate_flowclust_1d(fr=ff[idt, ], channel_bg, K=2)
  if(bg_filter@min<min.count) bg_filter@min <- min.count
  idt <- ff@exprs[, channel_bg]>bg_filter@min & ff@exprs[, chnl[2]]<chnl_filter@min

  peaks_chnl1 <- openCyto:::.find_peaks(ff@exprs[idt, chnl[1]], num_peaks = 1)
  peaks_chnl2 <- openCyto:::.find_peaks(ff@exprs[idt, chnl[2]], num_peaks = 1)

  k <- nrow(openCyto:::.find_peaks(ff@exprs[idt, chnl[1]]))*nrow(openCyto:::.find_peaks(ff@exprs[idt, chnl[2]]))

  filter <- openCyto::gate_flowclust_2d(fr=ff[idt, ], quantile = 0.99, xChannel = chnl[1], yChannel = chnl[2], K=6,target = c(peaks_chnl1$x[1], peaks_chnl2$x[1]))
  p1 <- ggcyto::autoplot(ff, x = chnl[2], y = channel_bg, bins = 100)+geom_gate(bg_filter)+geom_gate(chnl_filter)
  p2 <- ggcyto::autoplot(ff[idt, ], x = chnl[1], y = chnl[2], bins = 100)+geom_gate(filter)+
    xlim(range(ff@exprs[, chnl[1]]))+ylim(range(ff@exprs[, chnl[2]]))
  p3 <- ggcyto::autoplot(ff, x = chnl[1], y = chnl[2], bins = 100)+geom_gate(filter)

  p <- ggpubr::ggarrange(plotlist = list(as.ggplot(p1), as.ggplot(p2), as.ggplot(p3)))
  ggsave(paste0(output.dir, "/", filename, ".pdf"), plot = p, device = "pdf", width = 18, height = 12)
}


log_file <- "/Volumes/ExtremeSSD/IBIMA/Citometria/Log_file_error.txt"
log_file_error <- function(messages){
  sink(log_file, append = T)
  for(i in messages){
    cat(i, "\n")
  }
  sink()
}


metadata_panel1 <- metadata[metadata$panel%in%"Panel2", ]
rownames(metadata_panel1) <- 1:nrow(metadata_panel1)


cl <- parallel::makeCluster(6, setup_strategy = "sequential")
doParallel::registerDoParallel(cl)
library(doParallel)
system.time(clean <- foreach(i = 1:nrow(metadata_panel1)) %dopar% {
  tryCatch({
    library(flowCore)
    library(ggcyto)
    library(ggplot2)
    library(openCyto)
    filename_clean <- metadata_panel1$filename_clean[i]
    ff <- flowCore::read.FCS(filename_clean)
    channel_bg <- as.vector(ff@parameters@data$name[which(ff@parameters@data$desc=="CD4")])
    output.dir <- "/Volumes/TRANSCEND/IBIMA/Citometria/pruebas"
    filename <- unlist(strsplit(filename_clean, "/"))
    filename <- filename[length(filename)]
    filename <- gsub(".fcs$", "", filename)
    filename <- paste0(filename, "-", i)
    do_bg_pca(ff, filename = filename, channel_bg = channel_bg, output.dir = output.dir)
   },
  error=function(e) {
    log_file_error(paste(e, "Iter-->", i, "\n", "File: ", metadata$filename[i]))
    return(e)
  })

})

parallel::stopCluster(cl)


ggcyto::autoplot(ff[my_clusters$clusters==1, ], x = "FSC-A", y = "SSC-A", bins = 100)
ggcyto::autoplot(ff[my_clusters$clusters==2, ], x = "FSC-A", y = "SSC-A", bins = 100)
ggcyto::autoplot(ff[my_clusters$clusters==3, ], x = "FSC-A", y = "SSC-A", bins = 100)

ggcyto::autoplot(ff, x = "PE-Cy7-A", y = "SSC-A", bins = 100) + ggcyto::scale_x_logicle()
ggplot(as.data.frame(ff@exprs), aes( x = `PE-Cy7-A`, y = `SSC-A`)) + geom_point(aes(col=my_clusters$clusters))

ggcyto::autoplot(ff[my_clusters$clusters==1, ], x = "PE-Cy7-A", y = "SSC-A", bins = 100) + ggcyto::scale_x_logicle()
ggcyto::autoplot(ff[my_clusters$clusters==2, ], x = "PE-Cy7-A", y = "SSC-A", bins = 100) + ggcyto::scale_x_logicle()
ggcyto::autoplot(ff[my_clusters$clusters==3, ], x = "PE-Cy7-A", y = "SSC-A", bins = 100) + ggcyto::scale_x_logicle()



set.seed(123);idt <- sample(1:nrow(ff@exprs), nrow(ff@exprs)*0.3)
data <- data.frame(ff@exprs[idt, c("PE-Cy7-A", "FSC-A", "SSC-A")])
kms <- kmeans(x = data, centers = 3)

plot(data, col = kms$cluster, pch = 19)

cols <- viridis::viridis(kms$nclust)[kms$label]
rgl::plot3d(kms$x, col = cols)

res.pca <- prcomp(data, scale = TRUE)

km <- fpc::dbscan(data = res.pca$x, eps = 2)

factoextra::fviz_pca_ind(res.pca,
             col.ind = km$cluster, # Color by the quality of representation
             gradient.cols = c("#00AFBB", "#E7B800", "#FC4E07"),
             label = "none"     # Avoid text overlapping
)

rownames(data) <- paste0("Cell-", 1:nrow(data))
sade <- Seurat::CreateSeuratObject(t(data))

sade <- ScaleData(sade)
sade <- RunPCA(sade, features = rownames(sade))
sade <- FindNeighbors(object = sade, dims=1:2, reduction = "pca", force.recalc = F)
sade <- FindClusters(object = sade, resolution = 0.1, print.output = 0, random.seed = 123)
DimPlot(sade)
prono1<-knn(train = data, test = data, cl = tabla_alumnos$interes)
prono1

rgl::plot3d(x=data$FSC.A, y=data$SSC.A, z=data$PE.Cy7.A)
rgl::plot3d(x=H[, 1], y=H[, 2], z=H[, 3], col = c("red", "yellow", "orange"), size = 3, add = TRUE)

plot(x=data$FSC.A, y=data$SSC.A, size=1)
points(x=H[, 1], y=H[, 2], col="red", size = 10)

kms_iris <- ks::kms(x = ff@exprs[, c("PE-Cy7-A", "FSC-A", "SSC-A")], H = H)

kdr_oval <- ks::kdr(x = samp_oval)

emst <- emstreeR::ComputeMST(x = kdr_oval$end.points, verbose = FALSE)



plot(kms_iris,
     col = viridis::viridis(kms_iris$nclust))



set.seed(123)
idt <- sample(1:nrow(metadata_panel1), 50)
metadata_panel1 <- metadata_panel1[idt, ]

rownames(metadata_panel1) <- 1:nrow(metadata_panel1)
for(i in 1:nrow(metadata_panel1)){
  filename_clean <- metadata_panel1$filename_clean[i]
  ff <- flowCore::read.FCS(filename_clean)
  marker_bg <- c("CD3+")
  dir_plot <- paste0(output, gsub("/Volumes/ExtremeSSD/IBIMA/Citometria/fcs_clean/", "", filename_clean))
  dir_plot <- gsub(".fcs$", ".pdf", dir_plot)
  dummy <- do_backgating(ff, marker_bg = marker_bg, dir_plot=dir_plot)
}





dummy$plot
dummy$plot_back
dummy$plot_filter

