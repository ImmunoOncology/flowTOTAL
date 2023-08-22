
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
  labels[is.na(labels)] <- channels[is.na(labels)]
  
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
      
      adj <- 2
      if(traditional) adj <- 4
      
      k_peak <- openCyto:::.find_peaks(ff@exprs[idt, channel[i]], adjust = adj)
      k_peak <- k_peak[k_peak$y > max(k_peak$y)*0.05, ]
      k <- nrow(k_peak)
      
      idt_k_peak <- idt
      if(traditional){
        idt_k_peak <- idt_k_peak & ff@exprs[, channel[i]]<quantile(ff@exprs[, channel], 0.95)
        k_peak <- openCyto:::.find_peaks(ff@exprs[idt_k_peak, channel[i]], adjust = adj)
        k_peak <- k_peak[k_peak$y > max(k_peak$y)*0.05, ]
        k <- nrow(k_peak)
      }
      
      if(k==1 & first){
        
        k <- nrow(openCyto:::.find_peaks(ff@exprs[idt, channel[i]], adjust = 1))
        gate <- openCyto::gate_tail(fr = ff[idt, ], channel = channel[i], num_peaks = k, ref_peak = ifelse(k>1, 2, 1))
        
        if(traditional){
          k_peak <- openCyto:::.find_peaks(ff@exprs[idt, channel[i]], adjust = 1)
          k <- nrow(k_peak)
          if(k==1){
            gate <- openCyto::gate_tail(fr = ff[idt, ], channel = channel[i], num_peaks = k, ref_peak = k)
          }else{
            gate <- openCyto::gate_mindensity(fr=ff[idt, ], channel[i], cutpoint_min=cutpoint_min, cutpoint_max = cutpoint_max, min.count = 100, adjust = 1)
            #gate <- openCyto::gate_mindensity(fr=ff[idt, ], channel[i], cutpoint_min=cutpoint_min, cutpoint_max = cutpoint_max, min.count = 100, adjust = 0.9)
          }
        }
        
        gate@min <- ifelse(gate@min<cutpoint_min, cutpoint_min, gate@min)
        gate@min <- ifelse(gate@min>cutpoint_max, cutpoint_max, gate@min)
        
      }else if(k>2){
        
        gate <- openCyto::gate_tail(fr = ff[idt, ], channel = channel[i], num_peaks = k, ref_peak = ifelse(channel_sign[i]=="[-", k-1, k), side = ifelse(channel_sign[i]=="[-", "right", "left"))
        if(traditional){
          gate <- openCyto::gate_tail(fr = ff[idt, ], channel = channel[i], num_peaks = k, ref_peak = k)
        }
        
        gate@min <- ifelse(gate@min<cutpoint_min, cutpoint_min, gate@min)
        gate@min <- ifelse(gate@min>cutpoint_max, cutpoint_max, gate@min)
      }else if(k==2){
        gate <- openCyto::gate_mindensity(fr=ff[idt_k_peak, ], channel[i], cutpoint_min=cutpoint_min, cutpoint_max = cutpoint_max, min.count = 100)
        gate@min <- ifelse(gate@min>cutpoint_max, cutpoint_max, gate@min)
      }else{
        gate <- gate_bg[[channel[i]]]
      }
      
      if(beggining[i]){
        k_peak <- openCyto:::.find_peaks(ff@exprs[idt, channel[i]], adjust = 0.5)
        k_peak <- k_peak[order(k_peak$x, decreasing = F), ]
        rownames(k_peak) <- 1:nrow(k_peak)
        max_peak <- which.max(k_peak$y)
        gate <- openCyto::gate_tail(fr = ff[idt, ], channel = channel[i], side = "left", ref_peak = max_peak, num_peaks = nrow(k_peak), adjust = 0.5)
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

estimateProportion <- function(filename, id = NULL, panel, batch=NULL, output.dir, info_panel, cutpoint_min=1, cutpoint_max=500000){
  
  if(is.null(id)) id <- filename
  if(!dir.exists(paste0(output.dir, "/", panel))) dir.create(paste0(output.dir, "/", panel))
  
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
    
    if(!dir.exists(paste0(output.dir, "/", panel, "/", label))) dir.create(paste0(output.dir, "/", panel,  "/", label))
    flowCore::write.FCS(ff, filename = paste0(output.dir, "/", panel,  "/", label, "/", id, ".fcs"))
    
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