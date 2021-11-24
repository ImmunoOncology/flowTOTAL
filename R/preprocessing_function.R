
#' Preproccessing Function
#'
#' This function allows run preprocessing analysis for raw FCS file. It does the compensation and
#' the QC (remove doublets and anomalies)
#' @param file path to the file FCS to run.
#' @param filename filename for the cleaned FCS.
#' @param output path to the location for the cleaned FCS.
#' @keywords preprocessing
#' @export
#' @examples
#' run_Preprocessing()
run_Preprocessing <- function(file, filename, output){

  ff <- flowCore::read.FCS(file)
  ff_comp <- flowCore::compensate(ff, spillover = flowCore::spillover(ff)$SPILL)
  ff_QC <- flowAI::flow_auto_qc(ff_comp, ChExcludeFS = NULL, ChExcludeFM=NULL, html_report=F, fcs_QC=F, mini_report=F)
  ff_singlet <- PeacoQC::RemoveDoublets(ff_QC)
  flowCore::write.FCS(ff_singlet, filename = paste0(output, filename))

}


#' Filtering singlets
#'
#' Double checking for singlets by looking `FSC-A` and `FSC-H` channels
#' @param fC flowCore with the preprocess FCS data.
#' @param chnl channels used to identify singlets. Default  `FSC-A` and `FSC-H`.
#' @keywords singlets
#' @export
#' @examples
#' filter_singlets()
filter_singlets <- function(fC, chnl = c("FSC-A", "FSC-H")){
  gate_singlet <- openCyto:::.singletGate(fC, channels = chnl)
  idt_singlet <- filter(fC, gate_singlet)@subSet
  fC_singlet <- fC
  fC_singlet@exprs <- fC@exprs[idt_singlet, ]
  return(fC_singlet)
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


#' Do back-gating
#'
#' Do back-gating for the first filter.
#' @param fC flowCore to be simplified.
#' @param shape_channel channels used to identify shape and are. Default  `FSC-A`, `FSC-H`, `SSC-A` and `SSC-H`.
#' @keywords flowCore
#' @export
#' @examples
#' simplify_flowCore()
do_backgating <- function(fC, marker_bg, n.cutoff=200, current_channels = c("FSC-A", "SSC-A"), min.chnl1 = 25000, max.chnl1 = 150000, min.chnl1=NULL, min.chnl2 = NULL, max.chnl2 = 50000){

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

  x.axis <- fC_bg@exprs[, current_channels[1]]
  th.x.axis <- openCyto:::.find_peaks(x.axis, num_peaks = 3)
  th.x.axis <- th.x.axis[th.x.axis$x>min.fsca, ]
  max.fsca <- min(th.x.axis[which.max(th.x.axis$y), "x"], max.fsca)

  y.axis=fC_bg@exprs[, current_channels[2]]
  th.y.axis <- openCyto:::.find_peaks(y.axis, num_peaks = 2)
  th.y.axis <- th.y.axis[th.y.axis$x<max.ssca, ]
  max.ssca <- th.y.axis[which.max(th.y.axis$y), "x"]

  gate_target <- openCyto:::.flowClust.2d(fC, channels = current_channels, K=3, target=c(max.fsca,max.ssca), quantile=0.8)

  p <- ggcyto::autoplot(fC, x = current_channels[1], y = current_channels[2], bins = 100)
  p <- p + ggcyto::geom_gate(gate_target)+
    ggplot2::scale_y_continuous(labels=function(x)x/1000)+ggplot2::ylab("SSC-A  (x 1.000)")+
    ggplot2::scale_x_continuous(labels=function(x)x/1000)+ggplot2::xlab("FSC-A  (x 1.000)")+
    ggplot2::theme_bw()+ggplot2::theme(strip.text.x = ggplot2::element_text(size = 0, colour = "orange", angle = 90))

  to_filter <- filter(fC, gate_target)@subSet
  fC_filter <- fC
  fC_filter@exprs <- fC@exprs[to_filter, ]

  return(plot=p, fC_raw=fC, fC_filter=fC_filter)

}



do_backgating <- function(fC, marker_bg, n.cutoff=200, chnl = c("FSC-A", "SSC-A"), min.fsca = 25000, max.fsca = 150000, max.ssca = 50000){

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





