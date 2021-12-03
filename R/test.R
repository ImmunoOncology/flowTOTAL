
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

do_bg_pca <- function(ff, filename, output.dir, chnl = c("FSC-A", "SSC-A"), channel_bg){

  res.pca <- prcomp(ff@exprs[, c(channel_bg, chnl[1], chnl[2])], scale = TRUE)
  my_clusters <- ClusterR::KMeans_rcpp(res.pca$x, 3)
  my_clusters_bg <- ClusterR::KMeans_rcpp(ff@exprs[, c(channel_bg, chnl[2])], 3)
  bg_filter <- openCyto::gate_flowclust_1d(fr=ff, channel_bg, K=2)
  idt_cluster <- which((table(my_clusters$clusters[ff@exprs[, channel_bg]>bg_filter@min])/table(my_clusters$clusters))>0.9)
  list_centroid <- list()
  list_bg_filter <- list()

  # p1 <- ggcyto::autoplot(ff, x = chnl[1], y = chnl[2], bins = 100)
  # p2 <- ggcyto::autoplot(ff[my_clusters$clusters%in%idt_cluster, ], x = chnl[1], y = chnl[2], bins = 100)+
  #   xlim(range(ff@exprs[, chnl[1]]))+ylim(range(ff@exprs[, chnl[2]]))
  # for(cluster in idt_cluster){
  #   list_centroid[[cluster]] <- colMeans(ff@exprs[my_clusters$clusters==cluster, c(chnl[1], chnl[2])])
  #   list_bg_filter[[cluster]] <- openCyto::gate_flowclust_2d(fr=ff[my_clusters$clusters==cluster, ], quantile = 0.99, xChannel = chnl[1], yChannel = chnl[2], K=3,target = list_centroid[[cluster]])
  #   p1 <- p1 + ggcyto::geom_gate(list_bg_filter[[cluster]])
  #   p2 <- p2 + ggcyto::geom_gate(list_bg_filter[[cluster]])
  # }
  # p3 <- ggcyto::autoplot(ff, x = channel_bg, y = chnl[2], bins = 100)+ggcyto::geom_gate(bg_filter)
  # p4 <- ggcyto::autoplot(ff[my_clusters$clusters==idt_cluster, ], x = channel_bg, y = chnl[2], bins = 100)+ggcyto::geom_gate(bg_filter)+
  #   xlim(range(ff@exprs[, channel_bg]))+ylim(range(ff@exprs[, chnl[2]]))
  # p <- ggpubr::ggarrange(plotlist = list(as.ggplot(p1), as.ggplot(p2), as.ggplot(p3), as.ggplot(p4)))

  idt_bg <- ff@exprs[, channel_bg]>bg_filter@min
  idt <- which(table(my_clusters_bg$clusters[idt_bg], my_clusters$clusters[idt_bg])==max(table(my_clusters_bg$clusters[idt_bg], my_clusters$clusters[idt_bg])), arr.ind = TRUE)
  cells <- my_clusters$clusters==idt[, "row"] & my_clusters_bg$clusters==idt[, "col"]
  centroid <- colMeans(ff@exprs[cells & idt_bg, c(chnl[1], chnl[2])])
  bg_filter_v2 <- openCyto::gate_flowclust_2d(fr=ff[cells & idt_bg, ], quantile = 0.99, xChannel = chnl[1], yChannel = chnl[2], K=3,target = centroid)

  p1 <- ggcyto::autoplot(ff, x = chnl[1], y = chnl[2], bins = 100)+ggcyto::geom_gate(bg_filter_v2)
  p1 <- ggcyto::autoplot(ff[cells & idt_bg, ], x = chnl[1], y = chnl[2], bins = 100)+ggcyto::geom_gate(bg_filter_v2)
  p2 <- ggcyto::autoplot(ff, x = channel_bg, y = chnl[2], bins = 100)+ggcyto::geom_gate(bg_filter)
  idt <- sample(1:nrow(res.pca$x), 0.3*nrow(res.pca$x))
  p3 <- ggplot(as.data.frame(res.pca$x[idt, c(1,2)]), aes(x=PC1, y=PC2))+geom_point(aes(col=my_clusters$clusters[idt]))
  p4 <- ggplot(as.data.frame(ff@exprs)[idt, ], aes(x=`FSC-A`, y=`SSC-A`))+geom_point(aes(col=my_clusters$clusters[idt]))

  p <- ggpubr::ggarrange(plotlist = list(as.ggplot(p1), as.ggplot(p2), (p3), (p4)))
  ggsave(paste0(output.dir, "/", filename, ".pdf"), plot = p, device = "pdf", width = 18, height = 12)
}


rownames(metadata_panel1) <- 1:nrow(metadata_panel1)
for(i in 1:nrow(metadata_panel1)){
  filename_clean <- metadata_panel1$filename_clean[i]
  ff <- flowCore::read.FCS(filename_clean)
  channel_bg <- as.vector(ff@parameters@data$name[which(ff@parameters@data$desc=="CD3")])
  output.dir <- "/Volumes/TRANSCEND/IBIMA/Citometria/pruebas"
  filename <- unlist(strsplit(filename_clean, "/"))
  filename <- filename[length(filename)]
  filename <- gsub(".fcs$", "", filename)
  filename <- paste0(filename, "-", i)
  do_bg_pca(ff, filename = filename, channel_bg = channel_bg, output.dir = output.dir)
}




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

