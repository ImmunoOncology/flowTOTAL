flow_auto_qc_custom <- function (fcsfiles, filename="V1", timeCh = NULL,
                                 second_fractionFR = 0.1, alphaFR = 0.01, decompFR = TRUE,
                                 ChExcludeFS = c("FSC", "SSC"), outlier_binsFS = FALSE, pen_valueFS = 500,
                                 max_cptFS = 3, ChExcludeFM = c("FSC", "SSC"), sideFM = "both",
                                 neg_valuesFM = 1, mini_report = "QC_report", folder_results = "resultsQC") {
  
  FileType <- "FCS"
  i <- 1
  set <- as(fcsfiles, "flowSet")
  names <- flowCore::identifier(fcsfiles)
  
  if (missing(timeCh) || is.null(timeCh)) {
    timeCh <- flowAI:::findTimeChannel(set[[1]])
  }
  if (is.null(timeCh)) {
    warning("Impossible to retrieve the time channel automatically. The quality control can only be performed on signal acquisition and dynamic range.",
            call. = FALSE)
  }
  
  word <- which(grepl("TIMESTEP", names(flowCore::keyword(set[[1]])),
                      ignore.case = TRUE))
  timestep <- as.numeric(flowCore::keyword(set[[1]])[[word[1]]])
  if (!length(timestep)) {
    if (FileType == "LMD") {
      timestep <- 0.0009765625
    }
    else {
      warning("The TIMESTEP keyword was not found and hence it was set to 0.01. Graphs labels indicating time might not be correct",
              call. = FALSE)
      timestep <- 0.01
    }
  }
  if (second_fractionFR == "timestep") {
    second_fractionFR <- timestep
  }else if (second_fractionFR < timestep) {
    stop("The argument second_fractionFR must be greater or equal to timestep.",
         call. = FALSE)
  }
  if (folder_results != FALSE) {
    folder_results <- flowAI:::strip.sep(folder_results)
    dir.create(folder_results, showWarnings = FALSE)
    folder_results <- paste0(folder_results, .Platform$file.sep)
  }else {
    folder_results <- ""
  }
  filename_ext <- flowCore::identifier(set[[i]])
  filename <- sub("^([^.]*).*", "\\1", filename_ext)
  
  cat(paste0("Quality control for the file: ", filename,
             "\n"))
  
  if (!is.null(timeCh)) {
    if (length(unique(flowCore::exprs(set[[i]])[, timeCh])) ==
        1) {
      cat("The time channel contains a single value. It cannot be used to recreate the flow rate. \n")
      warning(paste0("The time channel in ", filename_ext,
                     " contains a single value. It cannot be used to recreate the flow rate. \n"),
              call. = FALSE)
      TimeChCheck <- "single_value"
    }
    else {
      TimeChCheck <- NULL
    }
  }else {
    TimeChCheck <- "NoTime"
  }
  
  FSbinSize <- min(max(1, ceiling(nrow(set[[1]])/100)),
                   500)
  if (is.null(TimeChCheck)) {
    ordFCS <- flowAI:::ord_fcs_time(set[[i]], timeCh)
  }else {
    ordFCS <- set[[i]]
  }
  
  origin_cellIDs <- 1:nrow(ordFCS)
  FR_bin_arg <- list(second_fraction = second_fractionFR,
                     timeCh = timeCh, timestep = timestep)
  FR_QC_arg <- list(alpha = alphaFR, use_decomp = decompFR)
  FS_bin_arg <- list(binSize = FSbinSize, timeCh = timeCh,
                     timestep = timestep, TimeChCheck = TimeChCheck)
  FS_QC_arg <- list(ChannelExclude = ChExcludeFS, pen_valueFS,
                    max_cptFS, outlier_binsFS)
  FM_QC_arg <- list(ChannelExclude = ChExcludeFM, side = sideFM,
                    neg_values = neg_valuesFM)
  if (is.null(TimeChCheck)) {
    FlowRateData <- do.call(flowAI:::flow_rate_bin, c(ordFCS,
                                                      FR_bin_arg))
    FlowRateQC <- do.call(flowAI:::flow_rate_check, c(ordFCS,
                                                      list(FlowRateData), FR_QC_arg))
  }else {
    FlowRateQC <- list()
    FlowRateQC$goodCellIDs <- origin_cellIDs
    FlowRateQC$res_fr_QC$badPerc <- 0
  }
  
  FlowSignalData <- do.call(flowAI:::flow_signal_bin, c(ordFCS,
                                                        FS_bin_arg))
  FlowSignalQC <- do.call(flowAI:::flow_signal_check, c(ordFCS,
                                                        list(FlowSignalData), FS_QC_arg))
  FlowMarginQC <- do.call(flowAI:::flow_margin_check, c(ordFCS,
                                                        FM_QC_arg))
  
  goodCellIDs <- intersect(FlowRateQC$goodCellIDs,
                           intersect(FlowSignalQC$goodCellIDs, FlowMarginQC$goodCellIDs))
  analysis <- "Flow Rate, Flow Signal and Flow Margin"
  
  badCellIDs <- setdiff(origin_cellIDs, goodCellIDs)
  totalBadPerc <- round(length(badCellIDs)/length(origin_cellIDs),
                        4)
  sub_exprs <- flowCore::exprs(ordFCS)
  params <-  flowCore::parameters(ordFCS)
  keyval <-  flowCore::keyword(ordFCS)
  
  
  goodfcs <- flowCore::flowFrame(exprs = sub_exprs[goodCellIDs, ], parameters = params, description = keyval)
  
  n_total <- as.integer(dim(set[[i]])[1])
  df_minireport <- data.frame(File=filename, N.initial.events=n_total,
                              FlowRateQC=n_total-length(FlowRateQC$goodCellIDs),
                              FlowSignalQC=n_total-length(FlowSignalQC$goodCellIDs),
                              FlowMarginQC=n_total-length(FlowMarginQC$goodCellIDs))
  
  return(list(FCS=goodfcs, minireport=df_minireport))
  
}

