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

  fC_simplify <- fC[, c(idt_shape, idt_shape)]
  return(fC_simplify)
}


do_backgating <- function(ff, filename, output.dir, chnl = c("FSC-A", "SSC-A"), channel_bg, channel_filter = c("SSC-A-"), cutpoint_min = 5000){

  channel_bg.raw <- channel_bg
  if(length(channel_bg)==1)
    channel_bg <- c(channel_bg, channel_filter)

  channel_sign <- paste0("[", sapply(channel_bg, function(x) substr(x, nchar(x), nchar(x))))
  channel_bg <-  sapply(channel_bg, function(x) substr(x, 1, nchar(x)-1))
  names(channel_bg) <- channel_bg

  `[+` <- function(a, b){ a>b}
  `[-` <- function(a, b){ a<b}

  gate_bg <- lapply(channel_bg, openCyto::gate_flowclust_1d, fr=ff, K=2, cutpoint_min=cutpoint_min)

  idt <- rep(TRUE, nrow(ff@exprs))

  for(i in length(gate_bg):1){
    gate <- gate_bg[[i]]
    events_bg <-  do.call(channel_sign[i], args = list(ff@exprs[, names(gate@min)], gate@min))
    idt <- idt & events_bg
  }

  # if(sum(idt)<200 & length(channel_filter)>1){
  #   i = 1
  #   message("Nº of events of channel ", channel_bg[i], " below cutpoint ", 100. )
  #   ff@parameters@data[which(ff@parameters@data$name==channel_bg[i]), "desc"] <- NA
  #   #channel_bg[i] <- ff@parameters@data$name[which(!is.na(ff@parameters@data$desc))[1]]
  #   #message("New channel selected --> ", channel_bg[i])
  #   #channel_sign <- gsub("[[]", "", channel_sign)
  #   #channel_bg[i] <- paste0(channel_bg[i], channel_sign[i])
  #   #channel_bg[-i] <- paste0(channel_bg[-i], channel_sign[-i])
  #   do_backgating(ff, filename = filename, channel_bg = channel_bg.raw, channel_filter = channel_filter[2], output.dir = output.dir)
  #   return(0)
  # }

  peaks_chnl1 <- openCyto:::.find_peaks(ff@exprs[idt, chnl[1]], num_peaks = 1)
  peaks_chnl1 <- openCyto:::.find_peaks(ff@exprs[idt, chnl[1]])
  peaks_chnl2 <- openCyto:::.find_peaks(ff@exprs[idt, chnl[2]], num_peaks = 1)
  peaks_chnl2 <- openCyto:::.find_peaks(ff@exprs[idt, chnl[2]])

  peaks_chnl1 <- peaks_chnl1[order(peaks_chnl1$y, decreasing = T), ]
  peaks_chnl2 <- peaks_chnl2[order(peaks_chnl2$x, decreasing = F), ]

  k <- nrow(openCyto:::.find_peaks(ff@exprs[idt, chnl[1]]))*nrow(openCyto:::.find_peaks(ff@exprs[idt, chnl[2]]))
  k <- max(3, k)
  k2 <- nrow(openCyto:::.find_peaks(ff@exprs[, chnl[1]]))*nrow(openCyto:::.find_peaks(ff@exprs[, chnl[2]]))
  k2 <- max(3, k2)
  idt_v2 <- rep(FALSE, nrow(ff@exprs))
  idt_v2[sample(1:length(idt_v2), 0.2*length(idt_v2))] <- TRUE
  idt_v2 <- idt_v2 | idt

  filter <- openCyto::gate_flowclust_2d(fr=ff[idt, ], quantile = 0.99, xChannel = chnl[1], yChannel = chnl[2], K=k,target = c(peaks_chnl1$x[1], peaks_chnl2$x[1]))
  filter_v2 <- openCyto::gate_flowclust_2d(fr=ff[idt_v2, ], quantile = 0.9, xChannel = chnl[1], yChannel = chnl[2], K=k,target = c(peaks_chnl1$x[1], peaks_chnl2$x[1]))
  filter_v3 <- openCyto::gate_flowclust_2d(fr=ff, quantile = 0.9, xChannel = chnl[1], yChannel = chnl[2], K=k2,target = c(filter_v2@mean))

  if(sum(idt)<200){
    filter_v2 <- openCyto::gate_flowclust_2d(fr=ff[idt_v2, ], quantile = 0.9, xChannel = chnl[1], yChannel = chnl[2], K=5,target = c(50000, 1000))
  }

  p1 <- ggcyto::autoplot(ff, x = channel_bg[1], y = channel_bg[2], bins = 100)
  for(gate in gate_bg){
    p1 <- p1 + ggcyto::geom_gate(gate)
  }

  p2 <- ggcyto::autoplot(ff[idt, ], x = chnl[1], y = chnl[2], bins = 100)+geom_gate(filter)+geom_gate(filter_v2, col = "blue")+geom_gate(filter_v3, col = "black")

  p3 <- ggcyto::autoplot(ff, x = chnl[1], y = chnl[2], bins = 100)+geom_gate(filter)+geom_gate(filter_v2, col = "blue")+geom_gate(filter_v3, col = "black")

  if(sum(filter(ff, filter_v3)@subSet)>sum(filter(ff, filter_v2)@subSet)){
    p4 <- ggcyto::autoplot(ff[filter(ff, filter_v3)@subSet, ], x = chnl[1], y = chnl[2], bins = 100)+ggcyto_par_set(limits = list(x = c(0,max(ff@exprs[, chnl[1]])), y = c(0, max(ff@exprs[, chnl[2]]))))#+xlim(0, max(ff@exprs[, chnl[1]]))+ylim(0, max(ff@exprs[, chnl[2]]))
  }else{
    p4 <- ggcyto::autoplot(ff[filter(ff, filter_v2)@subSet, ], x = chnl[1], y = chnl[2], bins = 100)+ggcyto_par_set(limits = list(x = c(0,max(ff@exprs[, chnl[1]])), y = c(0, max(ff@exprs[, chnl[2]]))))#+xlim(0, max(ff@exprs[, chnl[1]]))+ylim(0, max(ff@exprs[, chnl[2]]))
  }

  p <- ggpubr::ggarrange(plotlist = list(as.ggplot(p1), as.ggplot(p2), as.ggplot(p3), as.ggplot(p4)))
  ggsave(paste0(output.dir, "/", filename, ".pdf"), plot = p, device = "pdf", width = 18, height = 12)
}

 #' Create ggcyto plot
#'
#' Create a ggcyto plot based on two given channels (or one and by default `SSC-A` will be used as the second one). The ggcyto
#' could include gates and logicle scale.
#' @param fC flowCore to be simplified.
#' @param channels channels used to be plot. Vector with two elements. In case of one element, `SSC-A` will be used as the second channel.
#' @param gates gates to be included in the ggycyto.
#' @param logicle_chnls channels to be transfoormed in logicle scale.
#' @keywords plot
#' @export
#' @examples
#' do_ggcyto()
do_ggcyto <- function(fC, channels, gates, logicle_chnls){

  if(length(logicle_chnls)>0){
    lgcl <- logicleTransform()
    trans <- transformList(c(logicle_chnls), lgcl)
    fC_plot <- transform(fC, trans)
    gates <- ggcyto::rescale_gate(gates, lgcl, logicle_chnls)

    # for(i in 1:length(gates)){
    #   gates[[i]] <- ggcyto::rescale_gate(gates[[i]], lgcl, names(gates)[i])
    # }

  }else{
    fC_plot <- fC
  }

  if(length(channels)==1){
    channels <- c(channels, "SSC-A")
  }

  p <- ggcyto::autoplot(fC_plot, x = channels[1], y = channels[2], bins = 100)

  # for(i in 1:length(gates)){
  #   p <- p + ggcyto::geom_gate(gates[[i]])
  # }

  if(!is.null(gates))
    p <- p + ggcyto::geom_gate(gates)

  p <- p + ggplot2::theme_bw()+ggplot2::theme(strip.text.x = ggplot2::element_text(size = 0, colour = "orange", angle = 90))


  return(p)
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
do_backgating <- function(fC, marker_bg, n.cutoff=200, chnl = c("FSC-A", "SSC-A"),  min.chnl1 = 25000, max.chnl1 = 150000, min.chnl2 = 100, max.chnl2 = 50000,
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





