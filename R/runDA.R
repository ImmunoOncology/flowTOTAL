
doDA <- function(sce_file, output, response, response_label){

  dir.ouput <- file.path(output, "DA")
  if(!dir.exists(dir.ouput)) dir.create(dir.ouput)

  sce <- readRDS(sce_file)
  dummy_df <- as.data.frame(table(sce@metadata$Cluster, sce@metadata$ID)/rowSums(table(sce@metadata$Cluster, sce@metadata$ID)))
  p1 <- ggpubr::ggviolin(dummy_df, x="Var1", y = "Freq", add = "jitter")
  ggpubr::ggexport(filename = file.path(dir.ouput, "dummy_freq.pdf"), plotlist = list(p1), ncol = 1, nrow = 1)

  dummy_df <- as.data.frame(table(sce@metadata$Cluster, sce@metadata$ID))

  p2 <- ggpubr::ggviolin(dummy_df, x="Var1", y = "Freq", add = "jitter")
  ggpubr::ggexport(filename = file.path(dir.ouput, response, "dummy_count.pdf"), plotlist = list(p2), ncol = 1, nrow = 1)
  p3 <- ggpubr::ggviolin(dummy_df[dummy_df$Freq>10, ], x="Var1", y = "Freq", add = "jitter")
  ggpubr::ggexport(filename = file.path(dir.ouput, response, "dummy_count_nonzero.pdf"), plotlist = list(p3), ncol = 1, nrow = 1)

  p_cluster <- plotClusters(sce,
                            clusterColname = 'Cluster',
                            labSize = 7.0,
                            subtitle = 'UMAP performed on expression values',
                            caption = paste0('Note: clusters / communities identified via',
                                             '\nLouvain algorithm with multilevel refinement'),
                            axisLabSize = 20,
                            titleLabSize = 20,
                            subtitleLabSize = 16,
                            captionLabSize = 16)

  ggsave(plot = p_cluster, filename = file.path(dir.ouput, "p_cluster.pdf"), device = "pdf", height = 7, width = 7)

  if("Marker"%in%colnames(sce@metadata)){
    p_celltype <- plotClusters(sce,
                               clusterColname = 'Marker',
                               legendPosition = "bottom",
                               labSize = 7.0,
                               subtitle = 'UMAP performed on expression values',
                               caption = paste0('Note: clusters / communities identified via',
                                                '\nLouvain algorithm with multilevel refinement'),
                               axisLabSize = 20,
                               titleLabSize = 20,
                               subtitleLabSize = 16,
                               captionLabSize = 16)
    ggsave(plot = p_celltype, filename = file.path(dir.ouput, "p_Marker.pdf"), device = "pdf", height = 7, width = 7)

  }


  if(response%in%colnames(sce@metadata)){

    idt_response <- sce@metadata[, response]%in%response_label
    sce_response <- sce[, idt_response]
    sce_response@metadata <- sce@metadata[idt_response, ]

    dummy <- reshape2::melt(table(sce@metadata$Cluster, sce@metadata$ID))
    colnames(dummy) <- c("Cluster", "ID", "Value")
    dummy[, "Response"] <- sce@metadata[match(dummy$ID, sce@metadata$ID), response]
    dummy$Response <- factor(dummy$Response, levels = response_label)
    idt <- reshape2::dcast(data = dummy, formula = Cluster~Response, fun.aggregate = sum, value.var = "Value")
    rownames(idt) <- idt$Cluster
    idt$Cluster <- NULL
    dummy <- dummy[!dummy$Cluster%in%names(which(rowSums(idt)==0)), ]

    if(any(table(dummy$Response)==0)){
      return(NULL)
    }

    stat.test <- dummy %>% group_by(Cluster) %>%
      rstatix::wilcox_test(Value ~ Response) %>% as.data.frame()
    stat.test.sig <- stat.test[stat.test$p<0.05, ]

    p_all <- ggpubr::ggviolin(dummy, x = "Response" , y = "Value", color = "Response", add = "jitter", palette = "jco", facet.by = "Cluster")+ggpubr::stat_compare_means()
    ggpubr::ggexport(filename = file.path(dir.ouput, "p_all.pdf"), plotlist = list(p_all), ncol = 1, nrow = 1, height = 12, width =15)

    if(nrow(stat.test.sig)>0){
      p_stat <- ggpubr::ggviolin(dummy[dummy$Cluster%in%stat.test.sig$Cluster, ], x = "Response" , y = "Value", color = "Response", add = "jitter", palette = "jco", facet.by = "Cluster")+ggpubr::stat_compare_means()
      ggpubr::ggexport(filename = file.path(dir.ouput, "p_stat_sig.pdf"), plotlist = list(p_stat), ncol = 1, nrow = 1, height = 8, width = 5*nrow(stat.test.sig))
    }

    p_marker <- markerExpressionPerCluster(sce,
                                           clusters = c(unique(sce@metadata$Cluster)),
                                           clusterAssign = sce@metadata$Cluster,
                                           markers = rownames(sce),
                                           nrow = 3, ncol = round(length(unique(sce@metadata$Cluster))/3+0.4),
                                           caption = 'Cluster assignments based on UMAP performed on PC eigenvectors',
                                           stripLabSize = 22,
                                           axisLabSize = 22,
                                           titleLabSize = 22,
                                           subtitleLabSize = 18,
                                           captionLabSize = 18)

    write.table(stat.test, file.path(dir.ouput, "stat_test.txt"), col.names = T, row.names = F, sep = "\t", quote = F)
    ggsave(plot = p_marker, filename = file.path(dir.ouput, "p_marker.pdf"), device = "pdf", height = 12, width = 15)

  }else{
    stat.test <- NA
  }

  return(stat.test)

}

runDA <- function(output, response, response_label){
  doDA(file.path(output, "SEDA", "sce_CytoTree.rds"), output, response, response_label)
}
