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
#' @param shape_channel channels used to identify shape. Default  `FSC-A`, `FSC-H`, `SSC-A` and `SSC-H`.
#' @keywords flowCore
#' @export
#' @examples
#' simplify_flowCore()
simplify_flowCore <- function(filename, keep = NULL){
  fC <- read.FCS(filename)

  parameters_name <- names(fC@parameters@data$name)
  parameters_desc <- names(fC@parameters@data$desc)
  fC@parameters@data$name <- fC@parameters@data$desc
  names(fC@parameters@data$name) <- parameters_name
  names(fC@parameters@data$desc) <- parameters_desc

  colnames(fC@exprs) <- fC@parameters@data$name

  if(!is.null(keep)){
    if(all(keep%in%fC@parameters@data$name)){
      fC <- fC[, fC@parameters@data$name %in% keep]
      write.FCS(fC, filename)
    }else{
      return(FALSE)
    }
  }else{
    write.FCS(fC, filename)
  }

  return(TRUE)
}


#' Apply gate_flowclust_1d gate
#'
#' Apply gate_flowclust_1d gate to one ore two channels. If enougth informatio about the population is kept and join parameter is active
#' a target gate_flowclust_2d is applied to increase resolution.
#' @param fC flowCore to be gated.
#' @param channel channel to be used during gating.
#' @param cutpoint_min numeric value that sets a minimum thresold for the cutpoint. If a value is provided, any cutpoint below this value will be set to the given minimum value. If NULL (default), there is no minimum cutpoint value.
#' @param cutpoint_max numeric value that sets a maximum thresold for the cutpoint. If a value is provided, any cutpoint above this value will be set to the given maximum value. If NULL (default), there is no maximum cutpoint value.
#' @param return_plot whether the result should return the plot with the gates
#' @param join whether the result should be increased by determining a target populatin throught the given channels.
#' @keywords flowCore
#' @export
#' @examples
#' apply_gate()
apply_gate <- function(ff, channel, cutpoint_min, cutpoint_max, return_plot = TRUE, join=T, traditional=F){


  # Identify direction
  channel_sign <- paste0("[", sapply(channel, function(x) substr(x, nchar(x), nchar(x))))
  channel <-  sapply(channel, function(x) substr(x, 1, nchar(x)-1))

  beggining <- sapply(channel, function(x) substring(x, nchar(x), nchar(x)))=="^"
  ending <- sapply(channel, function(x) substring(x, nchar(x), nchar(x)))=="$"
  channel[beggining] <-  sapply(channel[beggining], function(x) substr(x, 1, nchar(x)-1))
  channel[ending] <-  sapply(channel[ending], function(x) substr(x, 1, nchar(x)-1))

  names(channel) <- channel
  channel <- unlist(channel)

  `[+` <- function(a, b){ a>b}
  `[-` <- function(a, b){ a<b}

  `[+]` <- function(a){ min(a)*1.05 }
  `[-]` <- function(a){ min(a)*1.05 }

  # `[+[` <- function(a){ max(a)*0.95 }
  # `[-[` <- function(a){ min(a)*1.05 }

  `[+[` <- function(a){
    peak_a <- openCyto:::.find_peaks(a, adjust = 1)
    peak_a <- peak_a[peak_a$y > max(peak_a$y)*0.01, ]
    return(peak_a$x[which.max(peak_a$x)])
  }  #{ quantile(a, 0.9) }#0.50) }
  `[-[` <- function(a){
    peak_a <- openCyto:::.find_peaks(a, adjust = 1)
    peak_a <- peak_a[peak_a$y > max(peak_a$y)*0.01, ]
    return(peak_a$x[which.max(peak_a$x)])
    }  #function(a){ quantile(a, 0.1) }#0.05) }

  idt <- rep(TRUE, nrow(ff@exprs))
  idt_list <- list()

  # Filter bg channel iteratively

  gate_bg <- list()


  num_peaks <- lapply(channel, function(x){
    k_peak <- openCyto:::.find_peaks(ff@exprs[idt, x], adjust = 2)
    k_peak <- k_peak[k_peak$y > max(k_peak$y)*0.05, ]
    return(nrow(k_peak))
  })



  first <- T
  max_iter <- 3
  iter <- 1

  while((any(unlist(num_peaks)>1) | first) & iter <= max_iter & (sum(idt)>200 | first)){

    for(i in length(channel):1){
      if(sum(idt)==0 | traditional){
        idt <- rep(TRUE, nrow(ff@exprs))
      }

      # gate <- tryCatch({
      #   openCyto::gate_flowclust_1d(fr=ff[idt, ], channel[i], K=2, cutpoint_min=cutpoint_min, cutpoint_max = cutpoint_max, min.count = 100)
      # },
      # error=function(e) {
      #   openCyto::gate_flowclust_1d(fr=ff[, ], channel[i], K=2, cutpoint_min=cutpoint_min, cutpoint_max = cutpoint_max, min.count = 100)
      # })


      k_peak <- openCyto:::.find_peaks(ff@exprs[idt, channel[i]], adjust = 2)
      k_peak <- k_peak[k_peak$y > max(k_peak$y)*0.05, ]
      k <- nrow(k_peak)
      if(k==1 & first){
        #gate <- openCyto::gate_flowclust_1d(fr=ff[idt, ], channel[i], K=k, cutpoint_min=cutpoint_min, cutpoint_max = cutpoint_max, min.count = 100)
        k <- nrow(openCyto:::.find_peaks(ff@exprs[idt, channel[i]], adjust = 1))
        gate <- openCyto::gate_tail(fr = ff[idt, ], channel = channel[i], num_peaks = k, ref_peak = ifelse(k>1, 2, 1))

        if(traditional){
          k_peak <- openCyto:::.find_peaks(ff@exprs[idt, channel[i]], adjust = 1)
          k <- nrow(k_peak)
          if(k==1){
            gate <- openCyto::gate_tail(fr = ff[idt, ], channel = channel[i], num_peaks = k, ref_peak = k)
          }else{
            #k_ref <- which(order(k_peak$x, decreasing = F)==2)
            #gate <- openCyto::gate_tail(fr = ff[idt, ], channel = channel[i], num_peaks = k, ref_peak = k_ref, side = "left")
            gate <- openCyto::gate_mindensity(fr=ff[idt, ], channel[i], cutpoint_min=cutpoint_min, cutpoint_max = cutpoint_max, min.count = 100, adjust = 1)
          }
        }

        gate@min <- ifelse(gate@min<cutpoint_min, cutpoint_min, gate@min)
        gate@min <- ifelse(gate@min>cutpoint_max, cutpoint_max, gate@min)

      }else if(k>2){
        #gate <- openCyto::gate_tail(fr = ff[idt, ], channel = channel[i], num_peaks = k, ref_peak = ifelse(channel_sign[i]=="[-", 1, k))


        #gate <- openCyto::gate_tail(fr = ff[idt, ], channel = channel[i], num_peaks = k, ref_peak = ifelse(channel_sign[i]=="[-", k-1, k), side = "left")

        gate <- openCyto::gate_tail(fr = ff[idt, ], channel = channel[i], num_peaks = k, ref_peak = ifelse(channel_sign[i]=="[-", k-1, k), side = ifelse(channel_sign[i]=="[-", "right", "left"))

        if(traditional){
          gate <- openCyto::gate_tail(fr = ff[idt, ], channel = channel[i], num_peaks = k, ref_peak = k)
        }

        gate@min <- ifelse(gate@min<cutpoint_min, cutpoint_min, gate@min)
        gate@min <- ifelse(gate@min>cutpoint_max, cutpoint_max, gate@min)
      }else if(k==2){
        gate <- openCyto::gate_mindensity(fr=ff[idt, ], channel[i], cutpoint_min=cutpoint_min, cutpoint_max = cutpoint_max, min.count = 100)
        gate@min <- ifelse(gate@min>cutpoint_max, cutpoint_max, gate@min)
      }else{
        gate <- gate_bg[[channel[i]]]
      }

      if(beggining[i]){
        gate <- openCyto::gate_tail(fr = ff[idt, ], channel = channel[i], side = "left")
      }

      if(ending[i]){
        gate <- openCyto::gate_tail(fr = ff[idt, ], channel = channel[i], side = "right")
      }

      events_bg <-  do.call(channel_sign[i], args = list(ff@exprs[, names(gate@min)], gate@min))
      idt <- idt & events_bg
      idt_list[[i]] <- idt
      gate_bg[[channel[i]]] <- gate
    }

    if(traditional){
      idt <- Reduce("&", idt_list)
    }

    if(sum(idt)>200){
      num_peaks <- lapply(channel, function(x){
        k_peak <- openCyto:::.find_peaks(ff@exprs[idt, x], adjust = 2)
        k_peak <- k_peak[k_peak$y > max(k_peak$y)*0.05, ]
        return(nrow(k_peak))
      })
    }

    first <- F
    iter <- iter+1
  }


  # if(sum(idt)<50  & join){
  #   idt <-idt_list[[1]]
  # }

  if(sum(idt)>200 & join){

    border <- list()
    target <- list()

    for(i in 1:length(channel)){
      border[[channel[i]]] <- eval(call(paste0(channel_sign[i], "]"), ff@exprs[idt, channel[i]]))
      target[[channel[i]]] <- eval(call(paste0(channel_sign[i], "["), ff@exprs[idt, channel[i]]))
    }

    border <- unlist(border)
    target <- unlist(target)

    sam <- data.frame(ff@exprs[idt, c(channel)], tms = Sys.time() + 1:sum(idt), id = gl(1, sum(idt)))
    sp::coordinates(sam) <- as.formula(paste0("~", paste0(colnames(sam)[1:2], collapse = "+")))
    tr <- trip(sam, c("tms", "id"))
    g <- tripGrid(tr)
    target <- coordinates(g)[which.max(g$z), ]

    gate_peak <- openCyto:::.find_peaks(ff@exprs[idt, channel[2]], adjust = 2)
    gate_peak <- gate_peak[gate_peak$y > max(gate_peak$y)*0.05, ]
    k_peak <- min(nrow(gate_peak), 2)
    gate_bg[["join"]] <- openCyto::gate_flowclust_2d(fr=ff[idt, ], quantile = 0.9, xChannel = channel[1], yChannel = channel[2], K=k_peak, target = target, min.count = 100,
                                                     min = border)
    idt <- idt & flowCore::filter(ff, gate_bg[["join"]])@subSet

    # idt <- tryCatch({
    #   gate_peak <- openCyto:::.find_peaks(ff@exprs[idt, channel[1]], adjust = 3)
    #   gate_bg[["join"]] <- openCyto::gate_flowclust_2d(fr=ff[idt, ], quantile = 0.9, xChannel = channel[1], yChannel = channel[2], K=2, target = target, min.count = 100,
    #                                                    min = border)
    #   idt <- idt & flowCore::filter(ff, gate_bg[["join"]])@subSet
    #   return(idt)
    #   },
    # error=function(e) {
    #   return(idt)
    # })

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


#' Backgating
#'
#' Apply backgating based on given channel and the shape and area to determining tune population.
#' @param fC flowCore to be gated.
#' @param filename filename to use for the new FCS file.
#' @param output.dir folder to saved the plots and the filtered FCS file.
#' @param chnl chnl with the shape and area.
#' @param channel_bg channel to be used during backgating.
#' @param channel_filter channel to be used in case channel_bg is onlye one. By default SSC-A with the - sign.
#' @param cutpoint_min numeric value that sets a minimum thresold for the cutpoint. If a value is provided, any cutpoint below this value will be set to the given minimum value. If NULL (default), there is no minimum cutpoint value.
#' @param cutpoint_max numeric value that sets a maximum thresold for the cutpoint. If a value is provided, any cutpoint above this value will be set to the given maximum value. If NULL (default), there is no maximum cutpoint value.
#' @param downsampling whether to downsample the flowCore object.
#' @param n_sample size for the downsample flowCore object.
#' @keywords flowCore
#' @export
#' @examples
#' do_backgating()
do_backgating <- function(ff, filename, output.dir, chnl = c("FSC-A", "SSC-A"), channel_bg, logicle_chnls=NULL, channel_filter = c("SSC-A-"), cutpoint_min = 1000, cutpoint_max = 50000, downsampling = FALSE, n_sample = NULL){

  idt_final <- c()

  # In case of 1 backgating channel based, the second channel to be used is SSC-A
  if(length(channel_bg)==1)
    channel_bg <- c(channel_bg, channel_filter)

  result_gate <- apply_gate(ff = ff, channel = channel_bg, cutpoint_min = cutpoint_min, cutpoint_max = cutpoint_max, return_plot = TRUE)

  idt_bg <- result_gate$idt
  gate_bg <- result_gate$gate_bg
  channel_bg <- result_gate$channel

  if(!is.null(logicle_chnls)){
    logicle_chnls <- channel_bg[logicle_chnls]
  }



  p1 <- do_ggcyto(ff, channel_bg, gate_bg, logicle_chnls = logicle_chnls)

  # idt_downsample <- rep(FALSE, nrow(ff@exprs))
  # idt_downsample[sample(1:length(idt_downsample), 0.3*length(idt_downsample))] <- TRUE
  # idt_downsample <- idt_downsample | idt_bg
  #
  # if(any(channel_bg %in% chnl)){
  #   chnl_dummy <- channel_bg[channel_bg %in% chnl]
  #   idt_downsample[ff@exprs[, chnl_dummy]<gate_bg[[chnl_dummy]]@min] <- FALSE
  # }
  #
  # peaks_chnl1 <- openCyto:::.find_peaks(ff@exprs[idt_bg, chnl[1]], adjust = 4)
  # peaks_chnl2 <- openCyto:::.find_peaks(ff@exprs[idt_bg, chnl[2]], adjust = 4)
  #
  # peaks_chnl1 <- peaks_chnl1[peaks_chnl1$y/max(peaks_chnl1$y)>0.1, ]
  #
  # cnt <- 1
  # gate_filter <- list()
  #
  # for(num_peak in 1:nrow(peaks_chnl1)){
  #
  #   # idt_downsample <- idt_downsample & (ff@exprs[, chnl[1]] > min(ff@exprs[idt_bg, chnl[1]]))
  #   # idt_downsample <- idt_downsample & (ff@exprs[, chnl[2]] > min(ff@exprs[idt_bg, chnl[2]]))
  #
  #   k <- max(nrow(openCyto:::.find_peaks(ff@exprs[idt_downsample, chnl[1]], adjust = 3)), nrow(openCyto:::.find_peaks(ff@exprs[idt_downsample, chnl[2]], adjust = 3)))
  #   k <- max(k, 3)
  #   k <- min(k, 6)
  #
  #   gate_filter[[cnt]] <- openCyto::gate_flowclust_2d(fr=ff[idt_downsample, ], quantile = 0.9, xChannel = chnl[1], yChannel = chnl[2], K=k, target = colMeans(ff@exprs[idt_bg, chnl]))
  #
  #   if(sum(idt_bg)>200){
  #
  #     border_min <- apply(ff@exprs[idt_bg, chnl], 2, min)
  #     border_max <- apply(ff@exprs[idt_bg, chnl], 2, max)
  #
  #     ff_iter <- ff[flowCore::filter(ff, gate_filter[[cnt]])@subSet | idt_bg, ]
  #     ff_iter <- ff_iter[ff_iter@exprs[, chnl[1]]>border_min[1] & ff_iter@exprs[, chnl[1]]<border_max[1], ]
  #     ff_iter <- ff_iter[ff_iter@exprs[, chnl[2]]>border_min[2] & ff_iter@exprs[, chnl[2]]<border_max[2], ]
  #
  #     gate_filter[[cnt]] <- openCyto::gate_flowclust_2d(fr=ff_iter, quantile = 0.8, xChannel = chnl[1], yChannel = chnl[2], K=1, target = colMeans(ff_iter@exprs[, chnl]), min = border_min, max = border_max, min.count = 100)
  #   }
  #
  #   if(length(idt_final)==0) idt_final <- flowCore::filter(ff, gate_filter[[cnt]])@subSet
  #
  #   idt_final <- idt_final | flowCore::filter(ff, gate_filter[[cnt]])@subSet
  #
  #   cnt <- cnt+1
  # }

  sam <- data.frame(ff@exprs[idt_bg, c(chnl, channel_bg[1])], tms = Sys.time() + 1:sum(idt_bg), id = gl(1, sum(idt_bg)))
  if(nrow(sam)<4){
    pt <- sam[which.max(sam$APC.A), c(1, 2)]
  }else{
    sp::coordinates(sam) <- as.formula(paste0("~", paste0(colnames(sam)[1:3], collapse = "+")))
    tr <- trip(sam, c("tms", "id"))
    g <- tripGrid(tr)
    pt <- coordinates(g)[which.max(g$z), ]
  }


  if(sum(idt_bg)<100){
    new_idt <- ff@exprs[, chnl[2]]<cutpoint_max
    pt <- apply(ff@exprs[new_idt, chnl], 2, median)

    num_peaks <- openCyto:::.find_peaks(ff@exprs[new_idt, chnl[2]], adjust = 1)
    num_peaks <- num_peaks[num_peaks$y > max(num_peaks$y)*0.05, ]
    k <- nrow(num_peaks)

    sam <- data.frame(ff@exprs[new_idt, c(chnl, channel_bg[1])], tms = Sys.time() + 1:sum(new_idt), id = gl(1, sum(new_idt)))
    sp::coordinates(sam) <- as.formula(paste0("~", paste0(colnames(sam)[1:3], collapse = "+")))
    tr <- trip(sam, c("tms", "id"))
    g <- tripGrid(tr)
    pt <- coordinates(g)[which.max(g$z), ]

    my_gate <- openCyto::gate_flowclust_2d(fr=ff[new_idt, ], quantile = 0.9, xChannel = chnl[1], yChannel = chnl[2], K=k, target = pt, min.count = 100)

    gate_filter <- list(my_gate)


  }else{


    # num_peaks_x <- openCyto:::.find_peaks(ff@exprs[idt_bg, chnl[1]], adjust = 1)
    # num_peaks_x <- num_peaks_x[num_peaks_x$y > max(num_peaks_x$y)*0.1, ]
    #
    # num_peaks_y <- openCyto:::.find_peaks(ff@exprs[idt_bg, chnl[2]], adjust = 1)
    # num_peaks_y <- num_peaks_y[num_peaks_y$y > max(num_peaks_y$y)*0.1, ]
    #
    # k <- max(nrow(num_peaks_y), nrow(num_peaks_x))
    # k <- min(3, k)

    my_gate <- openCyto::gate_flowclust_2d(fr=ff[idt_bg, ], quantile = 0.9, xChannel = chnl[1], yChannel = chnl[2], target = pt, min.count = 100)

    # # Paranoia (OJO AQUI)
    # num_peaks_x <- openCyto:::.find_peaks(ff@exprs[, chnl[1]], adjust = 1)
    # num_peaks_x <- num_peaks_x[num_peaks_x$y > max(num_peaks_x$y)*0.1, ]
    # num_peaks_x[nrow(num_peaks_x)+1, ] <- c(max(ff@exprs[, chnl[1]]), 1)
    # cnt <- 0
    # mid_down <- 0
    # for(i in 1:(nrow(num_peaks_x)-1)){
    #   mid <- num_peaks_x$x[i]+(num_peaks_x$x[i+1]-num_peaks_x$x[i])/2
    #   peak_idt <- ff@exprs[, chnl[1]]<mid & ff@exprs[, chnl[1]]>mid_down
    #   if(sum(peak_idt)>0){
    #     num_peaks_y <- openCyto:::.find_peaks(ff@exprs[peak_idt, chnl[2]], adjust = 1)
    #     num_peaks_y <- num_peaks_y[num_peaks_y$y > max(num_peaks_y$y)*0.1, ]
    #     cnt <- cnt + nrow(num_peaks_y)
    #     mid_down <- num_peaks_x$x[i]
    #   }
    # }
    # my_gate <- openCyto::gate_flowclust_2d(fr=ff, K=cnt, quantile = 0.9, xChannel = chnl[1], yChannel = chnl[2], target = pt, min.count = 100)

    #my_gate <- openCyto::gate_flowclust_2d(fr=ff[idt_bg, ], quantile = 0.9, xChannel = chnl[1], yChannel = chnl[2], target = pt, min.count = 100)

    # idt_downsample <- rep(FALSE, nrow(ff@exprs))
    # idt_downsample[sample(1:length(idt_downsample), 0.1*length(idt_downsample))] <- TRUE
    # idt_downsample <- idt_downsample | flowCore::filter(ff, gate)@subSet
    #
    # pt_target <- apply(ff@exprs[flowCore::filter(ff, my_gate)@subSet, chnl], 2, median)
    # max_target <- apply(ff@exprs[flowCore::filter(ff, my_gate)@subSet, chnl], 2, max)
    # min_target <- apply(ff@exprs[flowCore::filter(ff, my_gate)@subSet, chnl], 2, min)
    #
    # my_gate <- openCyto::gate_flowclust_2d(fr=ff[idt_downsample, ], quantile = 0.75, xChannel = chnl[1], yChannel = chnl[2], K=k, target = pt_target, min.count = 100
    #                                        #, min = min_target, max = max_target
    #                                        )
    gate_filter <- list(my_gate)

    num_peaks <- openCyto:::.find_peaks(ff@exprs[idt_bg, chnl[1]], adjust = 1)
    num_peaks <- num_peaks[num_peaks$y > max(num_peaks$y)*0.25, ]

    # if(nrow(num_peaks)>1){
    #   gate_filter <- list()
    #   cnt <- 1
    #   for(my_peak in num_peaks$x){
    #     pt_dummy <- pt
    #     pt_dummy[1] <- my_peak
    #     gate_filter[[cnt]] <- openCyto::gate_flowclust_2d(fr=ff[idt_bg, ], quantile = 0.9, xChannel = chnl[1], yChannel = chnl[2], target = pt_dummy, min.count = 100, K=nrow(num_peaks))
    #     cnt <- cnt+1
    #     }
    # }
  }


  if(length(gate_filter)>2){
    idt_final <- Reduce("|", lapply(gate_filter, function(x) flowCore::filter(ff, x)@subSet))
  }else{
    idt_final <- flowCore::filter(ff, gate_filter[[1]])@subSet
  }

  p2 <- do_ggcyto(ff[idt_final, ], channel_bg, gate_bg)
  p3 <- do_ggcyto(ff[idt_bg, ], chnl, gate_filter)
  p4 <- do_ggcyto(ff, chnl, gate_filter)

  do_cluster <- T
  if(do_cluster){

    idt_down <- sample(1:nrow(ff), nrow(ff)*0.3)
    idt_down <- idt_down | idt_final
    my.matrix <- ff@exprs[, unique(chnl, channel_bg)]
    rownames(my.matrix) <- paste0("cell", 1:nrow(my.matrix))
    my.matrix.raw <- my.matrix
    my.matrix <- my.matrix[idt_down, ]

    nn <- Seurat::FindNeighbors(my.matrix, distance.matrix = F,
                                k.param = 50, compute.SNN = T, prune.SNN = 1/15,
                                nn.method = "rann", annoy.metric = "euclidean",
                                nn.eps = 0, verbose = T, force.recalc = F)

    nnc <- Seurat::FindClusters(nn$snn, modularity.fxn = 1,
                                initial.membership = NULL, weights = NULL,
                                node.sizes = NULL, resolution = 1.2, method = "matrix",
                                algorithm = 1, n.start = 10, n.iter = 10,
                                random.seed = 0, group.singletons = T,
                                temp.file.location = NULL, edge.file.name = NULL,
                                verbose = T)

    dummy <- as.data.frame(cbind(my.matrix, label=nnc[, 1]))
    dummy$label <- factor(dummy$label)

    if(any(table(dummy$label[idt_final])/table(dummy$label)>0.05)){
      cluster <- names(which(table(dummy$label[idt_final])/table(dummy$label)>0.05))
      dummy <- dummy[dummy$label%in%cluster, ]
      idt_cluster <- rownames(my.matrix.raw)%in%rownames(dummy)
    }else{
      idt_cluster <- idt_final
    }

    dummy <- as.data.frame(cbind(my.matrix, label=nnc[, 1]))
    dummy <- dummy[sample(1:nrow(dummy), nrow(dummy)*0.1), ]
    p2_cluster <- ggplot(dummy)+geom_point(aes(x=`FSC-A`, y=`SSC-A`, col=label))
    dummy$label <- 0
    dummy$label[rownames(dummy)%in%rownames(my.matrix.raw)[idt_cluster]] <- 1
    p3_cluster <- ggplot(dummy)+geom_point(aes(x=`FSC-A`, y=`SSC-A`, col=label))
    my_gate <- openCyto::gate_flowclust_2d(fr=ff[idt_cluster, ], K=1, quantile = 0.9, xChannel = chnl[1], yChannel = chnl[2], target = pt, min.count = 100)
    p4_cluster <- do_ggcyto(ff, chnl, list(my_gate))
    idt_final <- idt_cluster
  }



  p <- ggpubr::ggarrange(plotlist = list(p1, p4, p2_cluster, p4_cluster), nrow = 2, ncol = 2)

  filename <- gsub(".fcs$", "", filename)

  ggsave(paste0(output.dir, "/PDF/", filename, ".pdf"), plot = p, device = "pdf", width = 15, height = 12)
  flowCore::write.FCS(ff[idt_final, ], filename = paste0(output.dir, "/FCS/", filename, ".fcs"))

  return(sum(idt_final))

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

    if(!is.null(gates)){
      for(gate_idt in logicle_chnls){
        logicle_chnls_tr <- logicle_chnls[logicle_chnls%in%names(gates[[gate_idt]]@min)]
        gates[[gate_idt]] <- ggcyto::rescale_gate(gates[[gate_idt]], lgcl, logicle_chnls_tr)
      }
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

      beggining <- sapply(channel, function(x) substring(x, nchar(x), nchar(x)))=="^"
      ending <- sapply(channel, function(x) substring(x, nchar(x), nchar(x)))=="$"
      channel[beggining] <-  sapply(channel[beggining], function(x) substr(x, 1, nchar(x)-1))
      channel[ending] <-  sapply(channel[ending], function(x) substr(x, 1, nchar(x)-1))

      channel_sign[beggining] <- paste0("^", channel_sign[beggining])
      channel_sign[ending] <- paste0("$", channel_sign[ending])

      names(channel) <- channel
      names(channel_sign) <- names(channel)

      if(all(!names(channel)%in%ff@parameters@data$desc)){
        message("Error --> ", id)
        return(0)
      }

      channel <- ff@parameters@data$name[match(names(channel), ff@parameters@data$desc)]
      names(channel) <- names(channel_sign)
      channel <- paste0(channel, channel_sign)
      res_gate <- apply_gate(ff, channel, cutpoint_min = cutpoint_min, cutpoint_max = cutpoint_max, return_plot = F, join = F, traditional = T)

      idt_bg <- res_gate$idt
      gate_bg <- res_gate$gate_bg
      channel_bg <- res_gate$channel

      if(length(gate_bg)==0) gate_bg <- NULL

      plot_list[[paste0(pattern, collapse = ":")]] <- do_ggcyto(ff, channel_bg, gate_bg, logicle_chnls = channel_bg, main = label)
      ff <- ff[idt_bg, ]

    }

    df_traditional[, label] <- nrow(ff)
    plot_traditional[[label]] <- plot_list[[length(plot_list)]]
  }

  nrow <- ifelse(length(plot_traditional)%%2==0, length(plot_traditional)/2, length(plot_traditional)/2+1)
  result_plot <- ggpubr::ggarrange(plotlist = plot_traditional, ncol = 2, nrow = nrow)
  result_filename <- paste0(output.dir, "/Traditional_counts_", panel, ".txt")

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




