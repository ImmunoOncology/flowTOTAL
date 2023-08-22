#' Filtering singlets
#'
#' Double checking for singlets by looking `FSC-A` and `FSC-H` channels.
#' @param fC flowCore with the preprocess FCS data.
#' @param chnl channels used to identify singlets. Default  `FSC-A` and `FSC-H`.
#' @keywords singlets
#' @export
#' @examples
#' filterSinglets()
filterSinglets <- function(fC, chnl = c("FSC-A", "FSC-H")){
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
#' doPreprocessing()
doPreprocessing <- function(file, filename, output, report=T){
  
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
  if("SPILL"%in%names(ff@description)){
    ff_comp <- flowCore::compensate(ff, spillover = flowCore::spillover(ff)$SPILL)
  }else{
    ff_comp <- ff
  }
  
  res_QC <- tryCatch({
    flow_auto_qc_custom(ff_comp, filename = filename, ChExcludeFS = NULL, ChExcludeFM=NULL, mini_report="Preprocessing", folder_results=paste0(output, "/resultsQC"))
  }, error=function(x){
    res_QC <- list(
      FCS=ff_comp,
      minireport=(data.frame("File"=NA, "N.initial.events"=NA, "FlowRateQC"=NA, "FlowSignalQC"=NA, "FlowMarginQC"=NA, "RemoveDoublets"=NA, "N.final.events"=NA))
    )
    return(res_QC)
  })
  ff_QC <- res_QC$FCS
  chnl <-  c("FSC-A", "FSC-H")[c("FSC-A", "FSC-H")%in%ff_QC@parameters@data$name]
  if(length(chnl)==2){
    ff_singlet <- filter_singlets(ff_QC, chnl = chnl)
  }else{
    ff_singlet <- ff_QC
  }
  
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



#' Run preproccessing Function
#'
#' This function allows run preprocessing analysis for raw FCS file. It does the compensation and
#' the QC (remove doublets and anomalies).
#' @param metadata data.frame  with filename column
#' @param filename filename for the cleaned FCS.
#' @param output path to the location for the cleaned FCS.
#' @param report should report number of anomalies and doublets. Default set to TRUE.
#' @param cluster run through parallel package
#' @param log_file Log file
#' @keywords preprocessing
#' @export
#' @examples
#' runPreprocessing()
runPreprocessing <- function(metadata, output, report=T, ncores=NULL, log_file="Log_file_error_raw2clean.txt"){
  
  log_file_error <- function(messages){
    sink(log_file, append = T)
    for(i in messages){
      cat(i, "\n")
    }
    sink()
  }
  
  if(!all(colnames(metadata)%in%c("filename")))
    stop("Error in metadata")
  
  if(is.null(ncores)){
    for(i in 1:nrow(metadata)){
      tryCatch({
        file <- metadata$filename[i]
        filename <- gsub("fcs_raw/", "", file)
        output <- "fcs_clean"
        doPreprocessing(file, filename, output, report=T)
      },
      error=function(e) {
        log_file_error(paste(e, "Iter-->", i, "\n", "File: ", metadata$filename[i]))
        return(e)
      })
    }
  }else{
    
    cl <- parallel::makeCluster(ncores)
    doParallel::registerDoParallel(cl)
    library(doParallel)
    
    foreach(i = 1:nrow(metadata)) %dopar% {
      tryCatch({
        file <- metadata$filename[i]
        filename <- gsub("fcs_raw/", "", file)
        output <- "fcs_clean"
        doPreprocessing(file, filename, output, report=T)
      },
      error=function(e) {
        log_file_error(paste(e, "Iter-->", i, "\n", "File: ", metadata$filename[i]))
        return(e)
      })
    }
    parallel::stopCluster(cl)
  }
  
}
