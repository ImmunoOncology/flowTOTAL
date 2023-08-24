#' Check FCS
#'
#' Remove empty channels
#' @param ff file to FCS
#' @keywords FCS
#' @export
#' @examples
#' check_FCS()
check_FCS <- function(ff){
  idt <- grep("SC-", ff@parameters@data$name, fixed = T)
  
  if(any(is.na(ff@parameters@data$desc[idt])))
    ff@parameters@data$desc[idt] <- ff@parameters@data$name[idt]
  
  idt <- grep(" ", ff@parameters@data$desc, fixed = T)
  
  if(length(idt)>0)
    ff@parameters@data$desc[idt]-unlist(lapply(strsplit(ff@parameters@data$desc[idt], " "), `[`, 1))
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
  tryCatch({
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
  },
  error=function(e) {
    file.remove(filename)
    return(F)
  })
  
  
  
  return(TRUE)
}


get_contour <- function(ff){
  
  mkseq <- function(my_vector){
    my_vector_cp <- my_vector
    i <- 1
    j <- i
    my_seq <- F
    while(i<=(length(my_vector)-1)){
      if(my_vector[i]==(my_vector[i+1]-1)){
        if(!my_seq){
          j <- i
        }
        my_seq <- T
        my_vector_cp[i+1] <- my_vector[j]
      }else{
        my_seq <- F
        j <- i+1
      }
      i <- i+1
    }
    return(unique(my_vector_cp))
  }
  
  samp <- ff@exprs[, c("FSC-A", "SSC-A")]
  Hns <- ks::Hns(x = samp, deriv.order=0)
  kdde_0 <- ks::kde(x = samp, H = diag(diag((Hns))))
  
  dummy <- unlist(lapply(1:99, function(pct_i) length(with(kdde_0,contourLines(x=eval.points[[1]],y=eval.points[[2]],
                                                                               z=estimate,levels=cont[paste0(pct_i, "%")][[1]])))))
  
  
  pcts <- which(dummy==max(dummy))
  
  if(length(mkseq(pcts))==1){
    dummy_idt <- dummy
    dummy_idt[pcts] <- 0
    if(max(dummy_idt)!=1)
      pcts <- c(pcts, which(dummy_idt==max(dummy_idt)))
  }
  
  
  pcts <- mkseq(pcts)
  
  
  contour.dummy <- do.call("c", lapply(pcts, function(pct_i) with(kdde_0,contourLines(x=eval.points[[1]],y=eval.points[[2]],
                                                                                      z=estimate,levels=cont[paste0(pct_i, "%")]))))
  contour.dummy <- contour.dummy[order(unlist(lapply(contour.dummy, `[`, 1)), decreasing = T)]
  point.dummy <- lapply(1:length(contour.dummy), function(x) which(sp::point.in.polygon(point.x = ff@exprs[, "FSC-A"], point.y = ff@exprs[, "SSC-A"], pol.x = contour.dummy[[x]]$x, pol.y = contour.dummy[[x]]$y)==1))
  contour.dummy <- contour.dummy[unlist(lapply(point.dummy, length))>250]
  point.dummy <- lapply(1:length(contour.dummy), function(x) which(sp::point.in.polygon(point.x = ff@exprs[, "FSC-A"], point.y = ff@exprs[, "SSC-A"], pol.x = contour.dummy[[x]]$x, pol.y = contour.dummy[[x]]$y)==1))
  
  idt <- 1:length(point.dummy)
  res <- unlist(lapply(idt, function(idt2){
    sum(unlist(lapply(idt[!idt%in%idt2], function(idt3) ifelse(length(intersect(point.dummy[[idt2]], point.dummy[[idt3]]))>0, 1, 0))))
  }))
  
  cnt_dummy <- 2
  while(length(res)==1 & cnt_dummy<=length(which(dummy==max(dummy)))){
    pcts <- which(dummy==max(dummy))[cnt_dummy]
    cnt_dummy <- cnt_dummy+1
    
    contour.dummy <- do.call("c", lapply(pcts, function(pct_i) with(kdde_0,contourLines(x=eval.points[[1]],y=eval.points[[2]],
                                                                                        z=estimate,levels=cont[paste0(pct_i, "%")]))))
    contour.dummy <- contour.dummy[order(unlist(lapply(contour.dummy, `[`, 1)), decreasing = T)]
    point.dummy <- lapply(1:length(contour.dummy), function(x) which(sp::point.in.polygon(point.x = ff@exprs[, "FSC-A"], point.y = ff@exprs[, "SSC-A"], pol.x = contour.dummy[[x]]$x, pol.y = contour.dummy[[x]]$y)==1))
    contour.dummy <- contour.dummy[unlist(lapply(point.dummy, length))>250]
    point.dummy <- lapply(1:length(contour.dummy), function(x) which(sp::point.in.polygon(point.x = ff@exprs[, "FSC-A"], point.y = ff@exprs[, "SSC-A"], pol.x = contour.dummy[[x]]$x, pol.y = contour.dummy[[x]]$y)==1))
    
    idt <- 1:length(point.dummy)
    res <- unlist(lapply(idt, function(idt2){
      sum(unlist(lapply(idt[!idt%in%idt2], function(idt3) ifelse(length(intersect(point.dummy[[idt2]], point.dummy[[idt3]]))>0, 1, 0))))
    }))
    
  }
  
  while(any(res>0)){
    contour.dummy <- contour.dummy[-which.max(res)]
    point.dummy <- lapply(1:length(contour.dummy), function(x) which(sp::point.in.polygon(point.x = ff@exprs[, "FSC-A"], point.y = ff@exprs[, "SSC-A"], pol.x = contour.dummy[[x]]$x, pol.y = contour.dummy[[x]]$y)==1))
    
    
    idt <- 1:length(point.dummy)
    res <- unlist(lapply(idt, function(idt2){
      sum(unlist(lapply(idt[!idt%in%idt2], function(idt3) ifelse(length(intersect(point.dummy[[idt2]], point.dummy[[idt3]]))>0, 1, 0))))
    }))
    
  }
  
  
  contour.plot <-do.call('rbind', lapply(contour.dummy, function(df) data.frame(x=df[["x"]], y=df[["y"]])))
  ggplotify::as.ggplot(function(x) {
    plot(kdde_0)
    points(contour.plot$x, contour.plot$y, col="orange", cex=0.3)
  })
  
  return(contour.dummy)
}

doDensityBackgating <- function(ff, filename, output.dir, chnl = c("FSC-A", "SSC-A"), channel_bg, logicle_chnls=NULL, sd.max_it=0.75, min.pct_it=0.01, target.fsc=50000, target.ssc=17500, min.ff_subset=250){
  
  ff.raw <- ff
  
  `[+` <- function(a, b){ a>b}
  `[-` <- function(a, b){ a<b}
  
  
  panel_channel <- sapply(channel_bg, function(x) substr(x, 1, nchar(x)-1))
  channel_sign <- paste0("[", sapply(channel_bg, function(x) substr(x, nchar(x), nchar(x))))
  
  gate_chnl <- list()
  keep <- rep(TRUE, nrow(ff))
  
  for(j in 1:length(panel_channel)){
    n_peaks <- openCyto:::.find_peaks(ff@exprs[keep, panel_channel[j]], adjust = 0.5)
    n_peaks <- n_peaks[order(n_peaks$x, decreasing = F), ]
    
    n_peaks_dummy <- n_peaks[n_peaks$y>max(n_peaks$y)*0.01, ]
    n_peaks_dummy <- n_peaks_dummy[order(n_peaks_dummy$x, decreasing = F), ]
    
    dummy <- density(ff@exprs[keep, panel_channel[j]], adjust = 0.5)
    
    while(nrow(n_peaks_dummy)>1 & all(sapply(dummy$y[dummy$x>n_peaks_dummy$x[1] & dummy$x<n_peaks_dummy$x[2]][-1], function(x) dummy$y[dummy$x>n_peaks_dummy$x[1] & dummy$x<n_peaks_dummy$x[2]][1]-x)<0)){
      n_peaks_dummy <- n_peaks_dummy[2:nrow(n_peaks_dummy), ]
    }
    
    if(n_peaks$x[1]!=n_peaks_dummy$x[1]){
      ref <- which(n_peaks$x==n_peaks_dummy$x[1])
    }else{
      ref <- 1
    }
    
    #gate_chnl[[j]] <- openCyto::gate_tail(ff[keep, ], panel_channel[j], adjust = 0.5, ref_peak = ref, auto_tol=T, num_peaks = nrow(n_peaks))
    gate_chnl[[j]] <- gate_tail_custom(ff[keep, ], panel_channel[j], adjust = 0.5, ref_peak = ref, auto_tol=T, num_peaks = nrow(n_peaks))
    
    if(nrow(n_peaks)==1){
      n_peaks <- openCyto:::.find_peaks(ff@exprs[, panel_channel[j]], adjust = 1)
      #gate_chnl[[j]] <- openCyto::gate_tail(ff, panel_channel[j], adjust = 1, ref_peak = 1, auto_tol=T, num_peaks = nrow(n_peaks))
      gate_chnl[[j]] <- gate_tail_custom(ff, panel_channel[j], adjust = 1, ref_peak = 1, auto_tol=T, num_peaks = nrow(n_peaks))
      
    }
    
    keep <- keep & do.call(channel_sign[j], args = list(ff@exprs[, panel_channel[j]], gate_chnl[[j]]@min))
  }
  
  ff_subset <- ff[keep, ]
  
  p1 <- ggcyto::autoplot(ff, "FSC-A", "SSC-A", bins=100)
  p_gate <- ggcyto::autoplot(ff, panel_channel[1], panel_channel[2], bins=100)
  for(gate in gate_chnl){
    p_gate <- p_gate+ggcyto::geom_gate(gate)
  }
  
  ff@exprs <- rbind(ff@exprs, ff_subset@exprs)
  
  contour.pct <- get_contour(ff)
  
  my_gates <- list()
  p2 <- ggcyto::autoplot(ff.raw, "FSC-A", "SSC-A", bins=100)
  
  gate_keep <- list()
  big.cells <- list()
  ssc.cells <- list()
  fsc.cells <- list()
  pct_it_cnt <- 1
  
  for(pct_it in 1:length(contour.pct)){
    inner <- sp::point.in.polygon(point.x = ff@exprs[, "FSC-A"], point.y = ff@exprs[, "SSC-A"], pol.x = contour.pct[[pct_it]]$x, pol.y = contour.pct[[pct_it]]$y)
    
    K <- openCyto:::.find_peaks(ff@exprs[inner==1, "FSC-A"], adjust = 2)
    K <- K[K$y>K$y*0.1, ]
    
    for(j in 1:nrow(K)){
      
      gate_j <- tryCatch({
        suppressMessages(openCyto::gate_flowclust_2d(ff[inner==1, ], yChannel = "SSC-A", xChannel = "FSC-A", K=nrow(K), target = c(K$x[j], median(ff@exprs[inner==1, "SSC-A"]))))
      }, error=function(x){
        return(NA)
      })
      
      if(!is.na(gate_j)){
        my_gates[[pct_it_cnt]] <- gate_j
        p2 <- p2+ggcyto::geom_gate(my_gates[[pct_it_cnt]])
        
        keep <- flowCore::filter(ff.raw, my_gates[[pct_it_cnt]])@subSet
        big.cells[[pct_it_cnt]] <- mean(ff.raw@exprs[keep, c("SSC-A")])
        ssc.cells[[pct_it_cnt]] <- mean(ff.raw@exprs[keep, c("SSC-A")])
        fsc.cells[[pct_it_cnt]] <- mean(ff.raw@exprs[keep, c("FSC-A")])
        
        keep_gate <- rep(TRUE, nrow(ff.raw))
        
        for(k in 1:length(panel_channel)){
          keep_gate <- keep_gate & do.call(channel_sign[k], args = list(ff.raw@exprs[, panel_channel[k]], gate_chnl[[k]]@min))
        }
        
        
        keep_chl <- factor(ifelse(keep_gate, 1, 0), levels = c("0","1"))
        gate_keep[[pct_it_cnt]] <- keep_chl[keep]
        
        pct_it_cnt <- pct_it_cnt+1
      }
    }
  }
  
  
  too.small <- which(unlist(fsc.cells)<10000 & unlist(ssc.cells)<10000)
  if(length(too.small)>0 & length(too.small)!=length(fsc.cells)){
    gate_keep <- gate_keep[-too.small]
    fsc.cells <- fsc.cells[-too.small]
    ssc.cells <- ssc.cells[-too.small]
    my_gates <- my_gates[-too.small]
  }else if(length(too.small)>0 & length(too.small)==length(fsc.cells)){
    too.small <- which(unlist(fsc.cells)<250 & unlist(ssc.cells)<250)
    if(length(too.small)>0  & length(too.small)!=length(fsc.cells)){
      gate_keep <- gate_keep[-too.small]
      fsc.cells <- fsc.cells[-too.small]
      ssc.cells <- ssc.cells[-too.small]
      my_gates <- my_gates[-too.small]
    }
  }
  
  
  
  max_it <- max(unlist(lapply(gate_keep, function(x) (table(x)/length(x))[2])))
  pct_it <- unlist(lapply(gate_keep, function(x) (table(x)/length(x))[2]))
  
  if(all(pct_it<min.pct_it) | nrow(ff_subset)<min.ff_subset){
    pct_it <- which.min(abs(unlist(fsc.cells)-target.fsc)+abs(unlist(ssc.cells)-target.ssc))
  }else if(any(pct_it>0.5)){
    pct_it <- which(pct_it>0.4)
  }else{
    pct_it <- which(pct_it>(max_it*sd.max_it) & pct_it>min.pct_it)
    
    # if(length(pct_it)==2){
    #   pct_it <- which.min(unlist(fsc.cells))
    # }else{
    # }
  }
  
  idt_final <- rep(FALSE, nrow(ff.raw))
  for(pct_it_gate in pct_it){
    p2 <- p2+ggcyto::geom_gate(my_gates[[pct_it_gate]], col="blue")
    idt_final <- idt_final | flowCore::filter(ff.raw, my_gates[[pct_it_gate]])@subSet
  }
  
  # saveRDS(as.ggplot(p2), paste0(output.dir, "/PDF/", gsub(".fcs", "", my_sample, fixed = T), ".RDS"))
  
  
  p_final <- ggcyto::autoplot(ff.raw[idt_final, ], panel_channel[1], panel_channel[2], bins=100)
  p <- ggpubr::ggarrange(plotlist = list(ggcyto::as.ggplot(p1), ggcyto::as.ggplot(p_gate), ggcyto::as.ggplot(p2), ggcyto::as.ggplot(p_final)), nrow = 2, ncol = 2)
  
  ggplot2::ggsave(paste0(output.dir, "/PDF/", gsub(".fcs", "", filename, fixed = T), ".pdf"), plot = p, device = "pdf", width = 15, height = 12)
  flowCore::write.FCS(ff.raw[idt_final, ], filename = paste0(output.dir, "/FCS/", filename))
  
  return(sum(idt_final))
}


runDensityBackgating <- function(metadata, output, chnl = c("FSC-A", "SSC-A"), channel_bg, logicle_chnls=NULL, sd.max_it=0.75, 
                                min.pct_it=0.01, target.fsc=50000, target.ssc=17500, min.ff_subset=250, 
                                log_file="Log_file_error.txt", log_file_traditional="Log_file_error_traditional.txt", track_file="Log_file_track.txt", ncores=NULL){
  
  output.dir <- paste0(output, "/Backgating")
  if(!dir.exists(output.dir)) dir.create(output.dir)
  
  log_file <- paste0(output, "/", log_file)
  log_file_error <- function(messages){
    sink(log_file, append = T)
    for(i in messages){
      cat(i, "\n")
    }
    sink()
  }
  
  
  log_file_traditional <- paste0(output, "/", log_file_traditional)
  log_file_error_traditional <- function(messages){
    sink(log_file_traditional, append = T)
    for(i in messages){
      cat(i, "\n")
    }
    sink()
  }
  
  track_file <- paste0(output, "/", track_file)
  log_file_track <- function(messages){
    sink(track_file, append = T)
    for(i in messages){
      cat(i, "\n")
    }
    sink()
  }
  

  if(!dir.exists(paste0(output.dir, "/PDF"))) dir.create(paste0(output.dir, "/PDF"))
  if(!dir.exists(paste0(output.dir, "/FCS"))) dir.create(paste0(output.dir, "/FCS"))
  
  if(!is.null(ncores)){
    
    cl <- parallel::makeCluster(ncores)
    doParallel::registerDoParallel(cl)
    library(doParallel)
    
    foreach(i = 1:nrow(metadata)) %dopar% {
      library(flowCore)
      
      filename_clean <- metadata$filename_clean[i]
      filename <- unlist(lapply(strsplit(filename_clean, "/"), function(x) x[length(x)]))
      
      tryCatch({
        
        ff <- flowCore::read.FCS(filename_clean)
        res_bg <- doDensityBackgating(ff, filename = filename, output.dir = output.dir, channel_bg = channel_bg, logicle_chnls = logicle_chnls
                                        , sd.max_it=sd.max_it, min.pct_it=min.pct_it, target.fsc=target.fsc, target.ssc=target.ssc, min.ff_subset=min.ff_subset)
        log_file_track(paste(filename_clean, filename, paste0(output.dir, "/FCS/", filename), nrow(ff@exprs), res_bg, "Completed", i, sep = "\t"))
        
      },
      error=function(e) {
        log_file_error(paste(e, "Iter-->", i, "\n", "File: ", metadata$filename[i]))
        log_file_track(paste(filename_clean, filename, "\t", NA, NA, "Error", i, sep = "\t"))
        return(e)
      })
    }
  }else{
    for(i in 1:nrow(metadata)){
      filename_clean <- metadata$filename_clean[i]
      filename <- unlist(lapply(strsplit(filename_clean, "/"), function(x) x[length(x)]))
      
      tryCatch({
        
        ff <- flowCore::read.FCS(filename_clean)
        res_bg <- doDensityBackgating(ff, filename = filename, output.dir = output.dir, channel_bg = channel_bg, logicle_chnls = logicle_chnls
                                      , sd.max_it=sd.max_it, min.pct_it=min.pct_it, target.fsc=target.fsc, target.ssc=target.ssc, min.ff_subset=min.ff_subset)
        log_file_track(paste(filename_clean, filename, paste0(output.dir, "/FCS/", filename), nrow(ff@exprs), res_bg, "Completed", i, sep = "\t"))
        
      },
      error=function(e) {
        log_file_error(paste(e, "Iter-->", i, "\n", "File: ", metadata$filename[i]))
        log_file_track(paste(filename_clean, filename, "\t", NA, NA, "Error", i, sep = "\t"))
        return(e)
      })
    }
  }
}
  