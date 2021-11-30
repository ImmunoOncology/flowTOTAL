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

    for(i in 1:length(gates)){
      gates[[i]] <- ggcyto::rescale_gate(gates[[i]], lgcl, names(gates)[i])
    }
  }else{
    fC_plot <- fC
  }

  if(length(channels)==1){
    channels <- c(channels, "SSC-A")
  }

  p <- ggcyto::autoplot(fC_plot, x = channels[1], y = channels[2], bins = 100)

  for(i in 1:length(gates)){
    p <- p + ggcyto::geom_gate(gates[[i]])
  }

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

  gate_bg <- lapply(channel_bg, openCyto::gate_flowclust_1d, fr=fC, K=2, cutpoint_min=cutpoint_min)
  names(gate_bg) <- channel_bg

  `[+` <- function(a, b){ a>b}
  `[-` <- function(a, b){ a<b}

  fC_bg <- fC

  p_bg <- do_ggcyto(fC, channels = channel_bg, gates = gate_bg, logicle_chnls = channel_bg)

  for(i in 1:length(gate_bg)){
    gate <- gate_bg[[i]]

    events_bg <-  do.call(marker_bg_sign[i], args = list(fC_bg@exprs[, names(gate@min)], gate@min))
    if(sum(events_bg)>n.cutoff)
      fC_bg@exprs <- fC_bg@exprs[events_bg, ]
  }


  x.axis <- fC_bg@exprs[, chnl[1]]
  if(is.null(min.chnl1)) min.chnl1 <- min(x.axis)
  if(is.null(max.chnl1)) max.chnl1 <- max(x.axis)

  th.x.axis <- openCyto:::.find_peaks(x.axis, num_peaks = expected_population_chnl1)
  th.x.axis <- th.x.axis[th.x.axis$x>min.chnl1 & th.x.axis$x<max.chnl1, ]
  max.chnl1 <- min(th.x.axis[which.max(th.x.axis$y), "x"], max.chnl1)

  y.axis=fC_bg@exprs[, chnl[2]]
  if(is.null(min.chnl2)) min.chnl2 <- min(y.axis)
  if(is.null(max.chnl2)) max.chnl2 <- max(y.axis)

  th.y.axis <- openCyto:::.find_peaks(y.axis, num_peaks = expected_population_chnl2)
  th.y.axis <- th.y.axis[th.y.axis$x>min.chnl2 & th.y.axis$x<max.chnl2, ]
  max.chnl2 <- th.y.axis[which.max(th.y.axis$y), "x"]

  gate_target <- openCyto:::.flowClust.2d(fC, channels = chnl, K=expected_population, target=c(max.chnl1,max.chnl2), quantile=0.8)

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


get_2d_gate <- function(fC_bg, chnl, min.chnl1, max.chnl1, min.chnl2, max.chnl2, expected_population_chnl1, expected_population_chnl2,
                        expected_population, fn_1, fn_2){

  x.axis <- fC_bg@exprs[, chnl[1]]
  if(is.null(min.chnl1)) min.chnl1 <- min(x.axis)
  if(is.null(max.chnl1)) max.chnl1 <- max(x.axis)

  th.x.axis <- openCyto:::.find_peaks(x.axis, num_peaks = expected_population_chnl1)
  th.x.axis <- th.x.axis[th.x.axis$x>min.chnl1 & th.x.axis$x<max.chnl1, ]
  th.chnl1 <- min(th.x.axis[fn_1(th.x.axis$y), "x"], max.chnl1)

  y.axis=fC_bg@exprs[, chnl[2]]
  if(is.null(min.chnl2)) min.chnl2 <- min(y.axis)
  if(is.null(max.chnl2)) max.chnl2 <- max(y.axis)

  th.y.axis <- openCyto:::.find_peaks(y.axis, num_peaks = expected_population_chnl2)
  th.y.axis <- th.y.axis[th.y.axis$x>min.chnl2 & th.y.axis$x<max.chnl2, ]
  th.chnl2 <- min(th.y.axis[fn_2(th.y.axis$y), "x"], max.chnl2)

  gate_target <- openCyto:::.flowClust.2d(fC, channels = chnl, K=expected_population, target=c(th.chnl1,th.chnl2), quantile=0.8)

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





