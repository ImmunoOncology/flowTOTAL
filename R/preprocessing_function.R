#' Filtering singlets
#'
#' Double checking for singlets by looking `FSC-A` and `FSC-H` channels.
#' @param fC flowCore with the preprocess FCS data.
#' @param chnl channels used to identify singlets. Default  `FSC-A` and `FSC-H`.
#' @keywords singlets
#' @export
#' @examples
#' filter_singlets()
filter_singlets <- function(fC, chnl = c("FSC-A", "FSC-H")){
  fC_PeacoQC <- PeacoQC::RemoveDoublets(fC, channel1 = chnl[1], channel2 = chnl[2])
  gate_singlet <- openCyto:::.singletGate(fC_PeacoQC, channels = chnl)
  idt_singlet <- flowCore::filter(fC_PeacoQC, gate_singlet)@subSet
  fC_singlet <- fC_PeacoQC
  fC_singlet@exprs <- fC_PeacoQC@exprs[idt_singlet, ]
  return(fC_singlet)
}


#' Preproccessing Function
#'
#' This function allows run preprocessing analysis for raw FCS file. It does the compensation and
#' the QC (remove doublets and anomalies).
#' @param file path to the file FCS to run.
#' @param filename filename for the cleaned FCS.
#' @param output path to the location for the cleaned FCS.
#' @param report should report number of anomalies and doublets. Default set to TRUE.
#' @keywords preprocessing
#' @export
#' @examples
#' run_Preprocessing()
run_Preprocessing <- function(file, filename, output, report=T){

  if(!dir.exists(output)){
    message("Creating directory -->", output)
    dir.create(output)
  }

  if(!grepl(".fcs$", filename)){
    message("Adding extension .fcs")
    filename <- paste0(filename, ".fcs")
  }

  ff <- flowCore::read.FCS(file)
  flowCore::identifier(ff) <- gsub(".fcs$", "", filename)
  ff_comp <- flowCore::compensate(ff, spillover = flowCore::spillover(ff)$SPILL)
  res_QC <- flow_auto_qc_custom(ff_comp, filename = filename, ChExcludeFS = NULL, ChExcludeFM=NULL, mini_report="Preprocessing", folder_results=paste0(output, "/resultsQC"))
  ff_QC <- res_QC$FCS
  ff_singlet <- filter_singlets(ff_QC)
  res_QC$minireport$RemoveDoublets <- nrow(ff_QC@exprs)-nrow(ff_singlet@exprs)
  res_QC$minireport$N.final.events <- nrow(ff_singlet@exprs)
  res_QC$minireport$File <- file

  flowCore::write.FCS(ff_singlet, filename = paste0(output, "/", filename))

  if(report){
    reporte_filename <- paste0(output, "/resultsQC/Preprocessing.txt")
    if(!dir.exists(paste0(output, "/resultsQC"))) dir.create(paste0(output, "/resultsQC"))
    write.table(
      x = res_QC$minireport,
      file = reporte_filename
      , col.names = !file.exists(reporte_filename)
      , sep = "\t"
      , row.names = F
      , append = file.exists(reporte_filename)
      , quote = F
    )
  }
}


#' Simplify flowCore object
#'
#' Remove empty channels
#' @param fC flowCore to be simplified.
#' @param shape_channel channels used to identify shape and are. Default  `FSC-A`, `FSC-H`, `SSC-A` and `SSC-H`.
#' @keywords flowCore
#' @export
#' @examples
#' simplify_flowCore()
simplify_flowCore <- function(fC, shape_channel = c("FSC-A", "FSC-H", "SSC-A", "SSC-H")){
  idt_channel <- which(!is.na(fC@parameters@data$desc))
  idt_shape <- which(fC@parameters@data$name%in%shape_channel)

  fC_simplify <- fC[, c(idt_channel, idt_shape)]
  return(fC_simplify)
}


apply_gate <- function(ff, channel, cutpoint_min, cutpoint_max, return_plot = TRUE, join=T){


  # Identify direction
  channel_sign <- paste0("[", sapply(channel, function(x) substr(x, nchar(x), nchar(x))))
  channel <-  sapply(channel, function(x) substr(x, 1, nchar(x)-1))
  names(channel) <- channel

  `[+` <- function(a, b){ a>b}
  `[-` <- function(a, b){ a<b}

  `[+]` <- function(a){ min(a)*1.05 }
  `[-]` <- function(a){ min(a)*1.05 }

  `[+[` <- function(a){ max(a)*0.95 }
  `[-[` <- function(a){ min(a)*1.05 }

  idt <- rep(TRUE, nrow(ff@exprs))

  # Filter bg channel iteratively
  gate_bg <- list()
  for(i in 1:length(channel)){
    gate <- openCyto::gate_flowclust_1d(fr=ff, channel[i], K=2, cutpoint_min=cutpoint_min, cutpoint_max = cutpoint_max)
    events_bg <-  do.call(channel_sign[i], args = list(ff@exprs[, names(gate@min)], gate@min))
    idt <- idt & events_bg
    gate_bg[[channel[i]]] <- gate
  }

  border <- list()
  target <- list()

  for(i in 1:length(channel)){
    border[[channel[i]]] <- eval(call(paste0(channel_sign[i], "]"), ff@exprs[idt, channel[i]]))
    target[[channel[i]]] <- eval(call(paste0(channel_sign[i], "["), ff@exprs[idt, channel[i]]))
  }

  border <- unlist(border)
  target <- unlist(target)

  if(sum(idt)>200 & join){
    gate_peak <- openCyto:::.find_peaks(ff@exprs[idt, channel[1]], adjust = 3)
    gate_bg[["join"]] <- openCyto::gate_flowclust_2d(fr=ff[idt, ], quantile = 0.9, xChannel = channel[1], yChannel = channel[2], K=2, target = target, min.count = 100,
                                                     min = border)
    idt <- idt & flowCore::filter(ff, gate_bg[["join"]])@subSet
  }

  result = list(idt=idt, channel=channel, gate_bg=gate_bg)

  if(return_plot){
    p1 <- ggcyto::autoplot(ff, channel, bins = 100)
    for(gate in gate_bg){
      p1 <- p1 + ggcyto::geom_gate(gate)
    }
    result[["plot"]] <- p1
  }

  return(result)
}

identify_cluster <- function(ff, chnl, idt){

  # Look for the maximun peak in chnl1 and the peak with smallest complex
  peaks_chnl1 <- openCyto:::.find_peaks(ff@exprs[, chnl[1]], adjust = 2)
  peaks_chnl1 <- peaks_chnl1[order(peaks_chnl1$y, decreasing = T), ]
  peaks_chnl2 <- openCyto:::.find_peaks(ff@exprs[, chnl[2]], adjust = 2)
  peaks_chnl2 <- peaks_chnl2[order(peaks_chnl2$x, decreasing = F), ]

  # First bg clusters
  k <- nrow(peaks_chnl1)*min(nrow(peaks_chnl2), 3)
  k <- max(k, 3)
  k <- min(k, 12)

  a=flowClust::flowClust(ff, varNames = chnl, K = k)
  a@label <- factor(a@label)
  target_cluster <- which.max(table(a@label[idt & !a@flagOutliers])/table(a@label[!a@flagOutliers]))
  target_cluster <- max(table(a@label[idt & !a@flagOutliers])/table(a@label[!a@flagOutliers]), na.rm = T)
  target_cluster <- names(table(a@label))[round(table(a@label[idt & !a@flagOutliers])/table(a@label[!a@flagOutliers]), 1)%in%round(target_cluster, 1)]

  my_cluster <- a@label%in%target_cluster
  my_cluster <- !a@flagOutliers & my_cluster

  #plot(a, data = ff,level = 0.8, z.cutoff = 0)

  return(my_cluster)
}

do_backgating_cluster <- function(ff, filename, output.dir, chnl = c("FSC-A", "SSC-A"), channel_bg, channel_filter = c("SSC-A-"), cutpoint_min = 1000, cutpoint_max = 50000, downsampling = FALSE, n_sample = NULL){

  idt_final <- c()

  # In case of 1 backgating channel based, the second channel to be used is SSC-A
  if(length(channel_bg)==1)
    channel_bg <- c(channel_bg, channel_filter)

  result_gate <- apply_gate(ff = ff, channel = channel_bg, cutpoint_min = cutpoint_min, cutpoint_max = cutpoint_max, return_plot = TRUE)
  idt_bg <- result_gate$idt
  gate_bg <- result_gate$gate_bg
  channel_bg <- result_gate$channel
  p1 <- do_ggcyto(ff, channel_bg, gate_bg, logicle_chnls = channel_bg[1])

  idt_downsample <- rep(FALSE, nrow(ff@exprs))
  idt_downsample[sample(1:length(idt_downsample), 0.3*length(idt_downsample))] <- TRUE
  idt_downsample <- idt_downsample | idt_bg

  peaks_chnl1 <- openCyto:::.find_peaks(ff@exprs[idt_bg, chnl[1]], adjust = 4)
  peaks_chnl2 <- openCyto:::.find_peaks(ff@exprs[idt_bg, chnl[2]], adjust = 4)

  cnt <- 1
  gate_filter <- list()

  for(num_peak in 1:nrow(peaks_chnl1)){

    # idt_downsample <- idt_downsample & (ff@exprs[, chnl[1]] > min(ff@exprs[idt_bg, chnl[1]]))
    # idt_downsample <- idt_downsample & (ff@exprs[, chnl[2]] > min(ff@exprs[idt_bg, chnl[2]]))

    k <- max(nrow(openCyto:::.find_peaks(ff@exprs[idt_bg, chnl[1]], adjust = 3)), nrow(openCyto:::.find_peaks(ff@exprs[idt_bg, chnl[2]], adjust = 3)))
    k <- max(k, 3)
    k <- min(k, 6)

    gate_filter[[cnt]] <- openCyto::gate_flowclust_2d(fr=ff[idt_downsample, ], quantile = 0.9, xChannel = chnl[1], yChannel = chnl[2], K=k, target = colMeans(ff@exprs[idt_bg, chnl]))

    if(sum(idt_bg)>200){

      border_min <- apply(ff@exprs[idt_bg, chnl], 2, min)
      border_max <- apply(ff@exprs[idt_bg, chnl], 2, max)

      ff_iter <- ff[flowCore::filter(ff, gate_filter[[cnt]])@subSet, ]
      ff_iter <- ff_iter[ff_iter@exprs[, chnl[1]]>border_min[1] & ff_iter@exprs[, chnl[1]]<border_max[1], ]
      ff_iter <- ff_iter[ff_iter@exprs[, chnl[2]]>border_min[2] & ff_iter@exprs[, chnl[2]]<border_max[2], ]

      gate_filter[[cnt]] <- openCyto::gate_flowclust_2d(fr=ff_iter, quantile = 0.8, xChannel = chnl[1], yChannel = chnl[2], K=1, target = colMeans(ff_iter@exprs[, chnl]), min = border_min, max = border_max, min.count = 100)
    }

    if(length(idt_final)==0) idt_final <- flowCore::filter(ff, gate_filter[[cnt]])@subSet

    idt_final <- idt_final | flowCore::filter(ff, gate_filter[[cnt]])@subSet

    cnt <- cnt+1
  }

  p2 <- do_ggcyto(ff[idt_final, ], channel_bg, gate_bg)
  p3 <- do_ggcyto(ff[idt_bg, ], chnl, gate_filter)
  p4 <- do_ggcyto(ff, chnl, gate_filter)

  p <- ggpubr::ggarrange(plotlist = list(p1, p2, p3, p4), nrow = 2, ncol = 2)
  #p <- as.ggplot(p1)
  ggsave(paste0(output.dir, "/PDF/", filename, ".pdf"), plot = p, device = "pdf", width = 15, height = 12)
  flowCore::write.FCS(ff[idt_final, ], filename = paste0(output.dir, "/FCS/", filename, ".fcs"))
}

do_backgating_cluster_old <- function(ff, filename, output.dir, chnl = c("FSC-A", "SSC-A"), channel_bg, channel_filter = c("SSC-A-"), cutpoint_min = 3000, cutpoint_max = 50000, downsampling = FALSE, n_sample = NULL){

  idt_final <- c()

  # In case of 1 backgating channel based, the second channel to be used is SSC-A
  if(length(channel_bg)==1)
    channel_bg <- c(channel_bg, channel_filter)

  result_gate <- apply_gate(ff = ff, channel = channel_bg, cutpoint_min = cutpoint_min, cutpoint_max = cutpoint_max, return_plot = TRUE)
  idt_bg <- result_gate$idt
  gate_bg <- result_gate$gate_bg
  channel_bg <- result_gate$channel
  p1 <- do_ggcyto(ff, channel_bg, gate_bg, logicle_chnls = channel_bg[1])

  if(sum(idt_bg)<200){

    idt_downsample <- rep(FALSE, nrow(ff@exprs))
    idt_downsample[sample(1:length(idt_downsample), 0.3*length(idt_downsample))] <- TRUE
    idt_downsample <- idt_downsample | idt_bg

    k <- nrow(openCyto:::.find_peaks(ff@exprs[idt_downsample, chnl[1]], adjust = 3))*nrow(openCyto:::.find_peaks(ff@exprs[idt_downsample, chnl[2]], adjust = 3))
    k <- max(k, 3)
    k <- min(k, 6)

    gate_filter <- openCyto::gate_flowclust_2d(fr=ff[idt_downsample, ], quantile = 0.9, xChannel = chnl[1], yChannel = chnl[2], K=k,target = c(50000, 12500))
    idt_final <-  flowCore::filter(ff, gate_filter)@subSet
    gate_filter <- list(gate_filter)

  }else{

    idt_cluster <- identify_cluster(ff = ff, chnl = chnl, idt = idt_bg)

    peaks_chnl1 <- openCyto:::.find_peaks(ff@exprs[idt_cluster, chnl[1]], adjust = 4)
    peaks_chnl2 <- openCyto:::.find_peaks(ff@exprs[idt_cluster, chnl[2]], adjust = 4)

    cnt <- 1
    gate_filter <- list()

    for(num_peak in 1:nrow(peaks_chnl1)){

      gate_iter <-  openCyto::gate_flowclust_2d(fr=ff[idt_bg & idt_cluster, ], quantile = 0.9, xChannel = chnl[1], yChannel = chnl[2], K=nrow(peaks_chnl1), target = c(peaks_chnl1$x[num_peak], peaks_chnl2$x[1]))
      idt_iter <- flowCore::filter(ff, gate_iter)@subSet

      idt_downsample <- rep(FALSE, nrow(ff@exprs))
      idt_downsample[sample(1:length(idt_downsample), 0.3*length(idt_downsample))] <- TRUE
      idt_downsample <- idt_downsample | idt_iter
      idt_downsample <- idt_downsample & (ff@exprs[, chnl[1]] > min(ff@exprs[idt_bg, chnl[1]]))
      idt_downsample <- idt_downsample & (ff@exprs[, chnl[2]] > min(ff@exprs[idt_bg, chnl[2]]))

      k <- nrow(openCyto:::.find_peaks(ff@exprs[idt_downsample, chnl[1]], adjust = 3))*nrow(openCyto:::.find_peaks(ff@exprs[idt_downsample, chnl[2]], adjust = 3))
      k <- max(k, 3)
      k <- min(k, 6)

      gate_filter[[cnt]] <- openCyto::gate_flowclust_2d(fr=ff[idt_downsample, ], quantile = 0.9, xChannel = chnl[1], yChannel = chnl[2], K=k, target = colMeans(ff@exprs[idt_iter, chnl]), min = apply(ff@exprs[idt_bg, chnl], 2, min))

      if(length(idt_final)==0) idt_final <- flowCore::filter(ff, gate_filter[[cnt]])@subSet

      idt_final <- idt_final | flowCore::filter(ff, gate_filter[[cnt]])@subSet

      cnt <- cnt+1
    }
    idt_bg <- idt_bg & idt_cluster
  }

  p2 <- do_ggcyto(ff[idt_final, ], channel_bg, gate_bg)
  p3 <- do_ggcyto(ff[idt_bg, ], chnl, gate_filter)
  p4 <- do_ggcyto(ff, chnl, gate_filter)

  p <- ggpubr::ggarrange(plotlist = list(p1, p2, p3, p4), nrow = 2, ncol = 2)
  #p <- as.ggplot(p1)
  ggsave(paste0(output.dir, "/PDF/", filename, ".pdf"), plot = p, device = "pdf", width = 15, height = 12)
  flowCore::write.FCS(ff[idt_final, ], filename = paste0(output.dir, "/FCS/", filename, ".fcs"))
}

#' Create ggcyto plot
#'
#' Create a ggcyto plot based on two given channels (or one and by default `SSC-A` will be used as the second one). The ggcyto
#' could include gates and logicle scale.
#' @param fC flowCore.
#' @param channels channels used to be plot. Vector with two elements. In case of one element, `SSC-A` will be used as the second channel.
#' @param gates gates to be included in the ggycyto.
#' @param logicle_chnls channels to be transfoormed in logicle scale.
#' @keywords plot
#' @export
#' @examples
#' do_ggcyto()

do_ggcyto <- function(fC, channels, gates = NULL, logicle_chnls = NULL, main = NULL, legend = FALSE){

  if(length(logicle_chnls)>0){
    lgcl <- flowCore::logicleTransform()
    trans <- flowCore::transformList(c(logicle_chnls), lgcl)
    fC_plot <- ggcyto::transform(fC, trans)

    for(gate_idt in logicle_chnls){
      logicle_chnls_tr <- logicle_chnls[logicle_chnls%in%names(gates[[gate_idt]]@min)]
      gates[[gate_idt]] <- ggcyto::rescale_gate(gates[[gate_idt]], lgcl, logicle_chnls_tr)
    }

    if("join" %in% names(gates)){
      if(any(logicle_chnls%in%names(gates[["join"]]@mean)))
        gates[["join"]] <- ggcyto::rescale_gate(gates[["join"]], lgcl, logicle_chnls)
    }

  }else{
    fC_plot <- fC
  }

  if(length(channels)==1){
    channels <- c(channels, "SSC-A")
  }

  p <- ggcyto::autoplot(fC_plot, x = channels[1], y = channels[2], bins = 100)

  if(!is.null(gates)){
    for(gate in gates){
      p <- p + ggcyto::geom_gate(gate)
    }
  }

  p <- p + ggplot2::theme_bw()+ggplot2::theme(strip.text.x = ggplot2::element_text(size = 0, colour = "orange", angle = 90))

  #Labelling
  p <- ggcyto::as.ggplot(p)
  labels <- fC@parameters@data$desc[match(channels, fC@parameters@data$name)]

  for(i in 1:length(labels)){
    if(labels[i]!=channels[i]){
      labels[i] <- paste0(labels[i], " (", channels[i], ")")
    }
  }

  labels[match(logicle_chnls, channels)] <- paste0(labels[match(logicle_chnls, channels)], " - logicleTransform")

  p <- p + xlab(labels[1]) + ylab(labels[2])

  if(!is.null(main))
    p <- p+ggtitle(label = main)+ggplot2::theme(plot.title = ggplot2::element_text(size = 15, face = "bold", hjust = 0.5))

  if(!legend){
    p <- p+ggplot2::theme(legend.position = "none")
  }

  return(p)
}



do_traditional <- function(filename, id = NULL, panel, batch=NULL, output.dir, info_panel, cutpoint_min=1, cutpoint_max=500000){

  if(is.null(id)) id <- filename

  ff.raw <- flowCore::read.FCS(filename)
  n.total <- nrow(ff.raw)

  if(!is.null(batch) & grepl("Batch", colnames(info_panel), ignore.case = T)){
    info_panel <- info_panel[info_panel$Batch==batch, ]
  }

  info_panel <- info_panel[info_panel$Panel==panel, ]

  df_traditional <- data.frame(File=filename, ID=id)
  df_traditional$N.total <- n.total
  df_traditional[, info_panel$Label] <- NA

  plot_traditional <- list()

  for(i in 1:nrow(info_panel)){

    pattern_list <- info_panel$Pattern[i]
    label <- info_panel$Label[i]

    pattern_list <- unlist(strsplit(pattern_list, "[.]"))
    pattern_list <- sapply(pattern_list, strsplit, "[;]")

    plot_list <- list()
    ff <- ff.raw

    for(pattern in pattern_list){

      channel_sign <- sapply(pattern, function(x) substring(x, nchar(x), nchar(x)))
      channel <- sapply(pattern, function(x) substring(x, 1, nchar(x)-1))
      names(channel) <- channel
      names(channel_sign) <- names(channel)
      all(channel_sign%in%c("+","-"))

      if(all(!names(channel)%in%ff@parameters@data$desc)){
        message("Error --> ", id)
        return(0)
      }

      channel <- ff@parameters@data$name[match(names(channel), ff@parameters@data$desc)]
      names(channel) <- names(channel_sign)
      channel <- paste0(channel, channel_sign)
      res_gate <- apply_gate(ff, channel, cutpoint_min = cutpoint_min, cutpoint_max = cutpoint_max, return_plot = F, join = F)

      idt_bg <- res_gate$idt
      gate_bg <- res_gate$gate_bg
      channel_bg <- res_gate$channel

      plot_list[[paste0(pattern, collapse = ":")]] <- do_ggcyto(ff, channel_bg, gate_bg, logicle_chnls = channel_bg, main = label)
      ff <- ff[idt_bg, ]

    }

    df_traditional[, label] <- nrow(ff)
    plot_traditional[[label]] <- plot_list[[length(plot_list)]]
  }

  nrow <- ifelse(length(plot_traditional)%%2==0, length(plot_traditional)/2, length(plot_traditional)/2+1)
  result_plot <- ggpubr::ggarrange(plotlist = plot_traditional, ncol = 2, nrow = nrow)
  result_filename <- paste0(output.dir, "/Traditional_counts.txt")

  if(!dir.exists(paste0(output.dir, "/plots"))) dir.create(paste0(output.dir, "/plots"))
  if(!dir.exists(paste0(output.dir, "/plots/", panel))) dir.create(paste0(output.dir, "/plots/", panel))
  ggsave(paste0(output.dir, "/plots/", panel, "/", gsub(".fcs$", "", id), ".pdf"), plot = result_plot, device = "pdf", width = 15, height = 6*nrow)

  df_traditional$Sys_Time <- Sys.time()

  write.table(
    x = df_traditional,
    file = result_filename
    , col.names = !file.exists(result_filename)
    , sep = "\t"
    , row.names = F
    , append = file.exists(result_filename)
    , quote = F
  )

}

#' Do back-gating
#'
#' Do back-gating for the first filter.
#' @param fC flowCore to be simplified.
#' @param shape_channel channels used to identify shape and are. Default  `FSC-A`, `FSC-H`, `SSC-A` and `SSC-H`.
#' @keywords flowCore
#' @export
#' @examples
#' do_backgating()
do_backgating_old <- function(fC, marker_bg, n.cutoff=200, chnl = c("FSC-A", "SSC-A"),  min.chnl1 = 25000, max.chnl1 = 150000, min.chnl2 = 100, max.chnl2 = 50000,
                          expected_population_chnl1=3, expected_population_chnl2=2, expected_population=3, logicle_x=TRUE, logicle_y=TRUE,
                          dir_plot=NULL, cutpoint_min=3000){

  marker_bg <- unlist(strsplit(marker_bg, ":"))
  marker_bg_sign <- paste0("[", substring(marker_bg, nchar(marker_bg)))
  marker_bg <- substring(marker_bg, 1, nchar(marker_bg)-1)

  if(!all(marker_bg%in%fC@parameters@data$desc)){
    message("Provide correct format for marker_bg!")
  }

  channel_bg <- as.vector(fC@parameters@data$name[match(marker_bg, fC@parameters@data$desc)])

  if(!all(marker_bg_sign%in%c("[+", "[-"))){
    message("Provide correct format for marker_bg!")
  }

  fC_bg <- fC

  # gate_bg <- lapply(channel_bg, openCyto::gate_flowclust_1d, fr=fC, K=2, cutpoint_min=cutpoint_min)
  # names(gate_bg) <- channel_bg
  #
  # `[+` <- function(a, b){ a>b}
  # `[-` <- function(a, b){ a<b}
  #
  # p_bg <- do_ggcyto(fC, channels = channel_bg, gates = gate_bg, logicle_chnls = channel_bg)
  #
  # for(i in 1:length(gate_bg)){
  #   gate <- gate_bg[[i]]
  #
  #   events_bg <-  do.call(marker_bg_sign[i], args = list(fC_bg@exprs[, names(gate@min)], gate@min))
  #   if(sum(events_bg)>n.cutoff)
  #     fC_bg@exprs <- fC_bg@exprs[events_bg, ]
  # }


  gate_bg <- get_2d_gate(fC = fC, chnl = channel_bg, min.chnl1 = cutpoint_min, max.chnl1 = NULL
                             , min.chnl2 = NULL, max.chnl2 = NULL, expected_population_chnl1 = expected_population_chnl1
                             , expected_population_chnl2 = expected_population_chnl2, expected_population = expected_population, conf=0.90)

  p_bg <- do_ggcyto(fC, channels = channel_bg, gates = gate_bg, logicle_chnls = channel_bg)

  to_filter <- filter(fC, gate_bg)@subSet
  fC_bg <- fC
  fC_bg@exprs <- fC@exprs[to_filter, ]

  gate_target <- get_2d_gate(fC = fC_bg, chnl = chnl, min.chnl1 = min.chnl1, max.chnl1 = max.chnl1
                            , min.chnl2 = min.chnl2, max.chnl2 = max.chnl2, expected_population_chnl1 = expected_population_chnl1
                            , expected_population_chnl2 = expected_population_chnl2, expected_population = expected_population, conf=0.95)

  p <- ggcyto::autoplot(fC, x = chnl[1], y = chnl[2], bins = 100)
  p <- p + ggcyto::geom_gate(gate_target)+
    ggplot2::xlab(chnl[1])+
    ggplot2::ylab(chnl[2])+
    ggplot2::theme_bw()+ggplot2::theme(strip.text.x = ggplot2::element_text(size = 0, colour = "orange", angle = 90))

  p_back <- ggcyto::autoplot(fC_bg, x = chnl[1], y = chnl[2], bins = 100)
  p_back <- p_back + ggcyto::geom_gate(gate_target)+
    ggplot2::xlab(chnl[1])+ggplot2::xlim(c(0, max(fC_bg@exprs[, chnl[1]])))+
    ggplot2::ylab(chnl[2])+ggplot2::ylim(c(0, max(fC_bg@exprs[, chnl[2]])))+
    ggplot2::theme_bw()+ggplot2::theme(strip.text.x = ggplot2::element_text(size = 0, colour = "orange", angle = 90))


  to_filter <- filter(fC, gate_target)@subSet
  fC_filter <- fC
  fC_filter@exprs <- fC@exprs[to_filter, ]

  if(!is.null(dir_plot)){
    plot_list <- ggpubr::ggarrange(plotlist = list(ggcyto::as.ggplot(p), ggcyto::as.ggplot(p_bg), ggcyto::as.ggplot(p_back)), nrow = 1)
    ggplot2::ggsave(filename = dir_plot, plot = plot_list, device = "pdf", width = 12, height = 5)
  }

  return(list(plot=p, plot_back=p_bg, plot_filter=p_back, fC_raw=fC, fC_filter=fC_filter))

}


get_2d_gate <- function(fC, chnl, min.chnl1, max.chnl1, min.chnl2, max.chnl2, expected_population_chnl1, expected_population_chnl2,
                        expected_population, conf){

  fC_bg <- fC
  gate_bg <- openCyto::gate_flowclust_1d(fC, chnl[1], K=2)
  fC_bg@exprs <- fC@exprs[fC@exprs[, names(gate_bg@min)]>=gate_bg@min, ]

  if(length(chnl)==1)
    chnl <- c(chnl, "SSC-A")

  x.axis <- fC_bg@exprs[, chnl[1]]
  if(is.null(min.chnl1)) min.chnl1 <- min(x.axis)
  if(is.null(max.chnl1)) max.chnl1 <- max(x.axis)

  th.x.axis <- openCyto:::.find_peaks(x.axis, num_peaks = expected_population_chnl1)
  th.x.axis <- th.x.axis[th.x.axis$x>min.chnl1 & th.x.axis$x<max.chnl1, ]
  th.chnl1 <- th.x.axis[which.max(th.x.axis$x), "x"]

  y.axis=fC_bg@exprs[, chnl[2]]
  if(is.null(min.chnl2)) min.chnl2 <- min(y.axis)
  if(is.null(max.chnl2)) max.chnl2 <- max(y.axis)

  th.y.axis <- openCyto:::.find_peaks(y.axis, num_peaks = expected_population_chnl2)
  th.y.axis <- th.y.axis[th.y.axis$x>min.chnl2 & th.y.axis$x<max.chnl2, ]
  th.chnl2 <- th.y.axis[which.min(th.y.axis$x), "x"]

  expected_population <- max(nrow(th.y.axis), nrow(th.x.axis))

  gate_target <- openCyto:::.flowClust.2d(fC, channels = chnl, K=3, target=c(th.chnl1,th.chnl2), quantile=conf)

  return(gate_target)

}


do_gating <- function(fC, marker_bg, n.cutoff=1000){

  marker_bg <- c("CD127+:CD25-")
  marker_bg <- unlist(strsplit(marker_bg, ":"))
  marker_bg_sign <- paste0("[", substring(marker_bg, nchar(marker_bg)))
  marker_bg <- substring(marker_bg, 1, nchar(marker_bg)-1)

  if(!all(marker_bg%in%fC@parameters@data$desc)){
    message("Provide correct format for marker_bg!")
  }

  channel_bg <- as.vector(fC@parameters@data$name[match(marker_bg, fC@parameters@data$desc)])

  if(!all(marker_bg_sign%in%c("+", "-"))){
    message("Provide correct format for marker_bg!")
  }

  gate_bg <- lapply(channel_bg, openCyto::gate_flowclust_1d, fr=fC, K=2)

  `[+` <- function(a, b){ a>b}
  `[-` <- function(a, b){ a<b}

  fC_bg <- fC

  for(i in 1:length(gate_bg)){
    gate <- gate_bg[[i]]
    events_bg <-  do.call(marker_bg_sign[i], args = list(fC_bg@exprs[, names(gate@min)], gate@min))
    if(sum(events_bg)>n.cutoff)
      fC_bg@exprs <- fC_bg@exprs[events_bg, ]
  }

  x.axis <- fC_bg@exprs[, chnl[1]]
  th.x.axis <- openCyto:::.find_peaks(x.axis, num_peaks = 3)
  th.x.axis <- th.x.axis[th.x.axis$x>min.fsca, ]
  max.fsca <- min(th.x.axis[which.max(th.x.axis$y), "x"], max.fsca)

  y.axis=fC_bg@exprs[, chnl[2]]
  th.y.axis <- openCyto:::.find_peaks(y.axis, num_peaks = 2)
  th.y.axis <- th.y.axis[th.y.axis$x<max.ssca, ]
  max.ssca <- th.y.axis[which.max(th.y.axis$y), "x"]

  gate_target <- openCyto:::.flowClust.2d(fC, channels = chnl, K=3, target=c(max.fsca,max.ssca), quantile=0.8)

  p <- ggcyto::autoplot(fC, x = chnl[1], y = chnl[2], bins = 100)
  p <- p + ggcyto::geom_gate(gate_target)+
    ggplot2::scale_y_continuous(labels=function(x)x/1000)+ggplot2::ylab("SSC-A  (x 1.000)")+
    ggplot2::scale_x_continuous(labels=function(x)x/1000)+ggplot2::xlab("FSC-A  (x 1.000)")+
    ggplot2::theme_bw()+ggplot2::theme(strip.text.x = ggplot2::element_text(size = 0, colour = "orange", angle = 90))


  to_filter <- filter(fC, gate_target)@subSet
  fC_filter <- fC
  fC_filter@exprs <- fC@exprs[to_filter, ]

  return(plot=p, fC_raw=fC, fC_filter=fC_filter)

}





