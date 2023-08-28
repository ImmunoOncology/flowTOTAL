#' Custom Flow Auto QC from flowAI
#'
#' This function performs quality control for flow cytometry data using custom criteria based on `flow_auto_qc` function from flowAI package.
#' It calculates flow rates, signal quality, and flow margins to assess the data quality.
#'
#' Monaco G, Chen H (2022). flowAI: automatic and interactive quality control for flow cytometry data. R package version 1.30.0, 10.18129/B9.bioc.flowAI.
#'
#' @param fcsfiles A flowSet containing raw FCS files to be analyzed.
#' @param filename A character string to identify the analysis.
#' @param timeCh Name of the time channel in the FCS files.
#' @param second_fractionFR The second fraction for flow rate calculation. Default is 0.1.
#' @param alphaFR The significance level for flow rate check. Default is 0.01.
#' @param decompFR Logical indicating whether to use decomposition for flow rate check. Default is TRUE.
#' @param ChExcludeFS Channels to be excluded from flow signal QC. Default is c("FSC", "SSC").
#' @param outlier_binsFS Logical indicating whether to identify outlier bins for flow signal QC. Default is FALSE.
#' @param pen_valueFS Penalty value for flow signal QC. Default is 500.
#' @param max_cptFS Maximum number of change points for flow signal QC. Default is 3.
#' @param ChExcludeFM Channels to be excluded from flow margin QC. Default is c("FSC", "SSC").
#' @param sideFM Side for flow margin QC. Default is "both".
#' @param neg_valuesFM Number of negative values allowed for flow margin QC. Default is 1.
#' @param mini_report Filename for the mini report.
#' @param folder_results Path to the folder for saving results.
#'
#' @return A list containing cleaned flow cytometry data (flowFrame) and a summary report.
#'
#' @import flowAI
#' @import flowCore
#' @export
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


#' Custom Gate with Tail from OpenCyto
#'
#' This function creates a custom gate with a tail based on specified parameters. It is based on `gate_tail` from OpenCyto pacakge.
#'
#' Finak, Greg, Frelinger, Jacob, Jiang, Wenxin, Newell, Evan W., Ramey, John, Davis, Mark M., Kalams, Spyros A., De Rosa, Stephen C., Gottardo, Raphael (2014). “OpenCyto: An Open Source Infrastructure for Scalable, Robust, Reproducible, and Automated, End-to-End Flow Cytometry Data Analysis.” PLoS Computational Biology, 10(8), e1003806.
#'
#' @param fr A flowFrame containing flow cytometry data.
#' @param channel Name of the channel to create the gate on.
#' @param filterId ID for the created gate.
#' @param num_peaks Number of peaks to consider for tail gate calculation.
#' @param ref_peak Reference peak for tail gate calculation.
#' @param strict Logical indicating whether to strictly follow the cutoff point or not.
#' @param tol Tolerance for calculating the cutoff point.
#' @param side Side for gate creation ("left" or "right").
#' @param min Minimum value for truncating data before gate creation.
#' @param max Maximum value for truncating data before gate creation.
#' @param bias Bias to adjust the cutoff point.
#' @param positive Logical indicating whether gate should be positive or not.
#' @param ... arguments to pass to cytokine_cutpoint
#'
#' @return A rectangleGate representing the custom gate.
#'
#' @import flowCore
#' @import openCyto
#' @export
gate_tail_custom <- function(fr, channel, filterId = "", num_peaks = 1,
                      ref_peak = 1, strict = TRUE, tol = 1e-2, side = "right", min = NULL, max = NULL, bias = 0, positive = TRUE, ...) {

  side <- match.arg(side, c("right", "left"))
  if (!(is.null(min) && is.null(max))) {
    fr <- openCyto:::.truncate_flowframe(fr, channels = channel, min = min,
                              max = max)
  }
  # cutpoint is calculated using the first derivative of the kernel density
  # estimate.
  x <- as.vector(fr@exprs[, channel])
  cutpoint <- cytokine_cutpoint(x = x, num_peaks = num_peaks,
                                 ref_peak = ref_peak, tol = tol, side = side, strict = strict, ...)

  cutpoint <- cutpoint + bias
  if(positive){
    gate_coordinates <- list(c(cutpoint, Inf))
  }else{
    gate_coordinates <- list(c(-Inf, cutpoint))
  }
  names(gate_coordinates) <- channel
  flowCore::rectangleGate(gate_coordinates, filterId = filterId)

}

#' Calculate Cytokine Cutpoint from OpenCyto
#'
#' This function calculates the cutpoint for a cytokine gate based on specified parameters. It is based on `.cytokine_cutpoint` from OpenCyto pacakge.
#'
#' Finak, Greg, Frelinger, Jacob, Jiang, Wenxin, Newell, Evan W., Ramey, John, Davis, Mark M., Kalams, Spyros A., De Rosa, Stephen C., Gottardo, Raphael (2014). “OpenCyto: An Open Source Infrastructure for Scalable, Robust, Reproducible, and Automated, End-to-End Flow Cytometry Data Analysis.” PLoS Computational Biology, 10(8), e1003806.
#'
#' @param x A numeric vector containing data values.
#' @param num_peaks Number of peaks to consider for cutpoint calculation.
#' @param ref_peak Reference peak for cutpoint calculation.
#' @param method Method for calculating cutpoint ("first_deriv" or "second_deriv").
#' @param tol Tolerance for calculating the cutoff point.
#' @param adjust Adjustment factor for peak detection.
#' @param side Side for cutpoint calculation ("left" or "right").
#' @param strict Logical indicating whether to strictly follow the cutoff point or not.
#' @param plot Logical indicating whether to plot the peak detection.
#' @param auto_tol Logical indicating whether to automatically set tolerance based on peak values.
#' @param ... parameter to pass to deriv_density
#'
#' @return The calculated cutpoint.
#'
#' @import openCyto
#' @export
cytokine_cutpoint <- function(x, num_peaks = 1, ref_peak = 1,
                               method = c("first_deriv", "second_deriv"),
                               tol = 1e-2, adjust = 1, side = "right", strict = TRUE, plot = FALSE, auto_tol = FALSE, ...) {

  method <- match.arg(method)
  peaks <- sort(openCyto:::.find_peaks(x, num_peaks = num_peaks, adjust = adjust, plot = plot)[, "x"])

  #update peak count since it can be less than num_peaks
  num_peaks <- length(peaks)

  if (ref_peak > num_peaks) {
    outFunc <- ifelse(strict, stop, warning)
    outFunc("The reference peak is larger than the number of peaks found.",
            "Setting the reference peak to 'num_peaks'...",
            call. = FALSE)
    ref_peak <- num_peaks
  }

  # Double-check that a cutpoint minimum found via 'first_deriv'
  # passes the second-derivative test.

  if (method == "first_deriv") {
    # Finds the deepest valleys from the kernel density and sorts them.
    # The number of valleys identified is determined by 'num_peaks'
    deriv_out <- deriv_density(x = x, adjust = adjust, deriv = 1, ...)
    if(auto_tol){
      #Try to set the tolerance automatigically.
      tol = 0.01*max(abs(deriv_out$y))
    }
    if (side == "right") {

      deriv_valleys <- with(deriv_out, openCyto:::.find_valleys(x = x, y = y, adjust = adjust))
      deriv_valleys <- deriv_valleys[deriv_valleys > peaks[ref_peak]]
      deriv_valleys <- sort(deriv_valleys)[1]
      cutpoint <- with(deriv_out, x[x > deriv_valleys & abs(y) < tol])
      cutpoint <- cutpoint[1]

    } else if (side == "left") {

      deriv_out$y <- -deriv_out$y
      deriv_valleys <- with(deriv_out, openCyto:::.find_valleys(x = x, y = y, adjust = adjust))
      deriv_valleys <- deriv_valleys[deriv_valleys < peaks[ref_peak]]
      deriv_valleys <- sort(deriv_valleys, decreasing=TRUE)[1]
      cutpoint <- with(deriv_out, x[x < deriv_valleys & abs(y) < tol])
      cutpoint <- cutpoint[ length(cutpoint) ]

    } else {
      stop("Unrecognized 'side' argument (was '", side, "'.")
    }

  } else {
    # The cutpoint is selected as the first peak from the second derivative
    # density which is to the right of the reference peak.
    deriv_out <- deriv_density(x = x, adjust = adjust, deriv = 2, ...)

    if (side == "right") {
      deriv_peaks <- with(deriv_out, openCyto:::.find_peaks(x, y, adjust = adjust)[, "x"])
      deriv_peaks <- deriv_peaks[deriv_peaks > peaks[ref_peak]]
      cutpoint <- sort(deriv_peaks)[1]
    } else if (side == "left") {
      deriv_out$y <- -deriv_out$y
      deriv_peaks <- with(deriv_out, openCyto:::.find_peaks(x, y, adjust = adjust)[, "x"])
      deriv_peaks <- deriv_peaks[deriv_peaks < peaks[ref_peak]]
      cutpoint <- sort(deriv_peaks, decreasing=TRUE)[length(deriv_peaks)]
    } else {
      stop("Unrecognized 'side' argument (was '", side, "'.")
    }

  }

  cutpoint
}

#' Calculate Derivative of Density from OpenCyto
#'
#' This function calculates the derivative of the density function. It is based on `.deriv_density` from OpenCyto pacakge.
#'
#' Finak, Greg, Frelinger, Jacob, Jiang, Wenxin, Newell, Evan W., Ramey, John, Davis, Mark M., Kalams, Spyros A., De Rosa, Stephen C., Gottardo, Raphael (2014). “OpenCyto: An Open Source Infrastructure for Scalable, Robust, Reproducible, and Automated, End-to-End Flow Cytometry Data Analysis.” PLoS Computational Biology, 10(8), e1003806.
#'
#' @param x A numeric vector containing data values.
#' @param deriv Order of the derivative.
#' @param bandwidth Bandwidth for kernel density estimation.
#' @param adjust Adjustment factor for bandwidth.
#' @param num_points Number of points for estimation.
#' @param ... parameter to pass to drvkde
#'
#' @return A list with x and y components representing the calculated derivative.
#'
#' @import flowStats
#' @export
deriv_density <- function(x, deriv = 1, bandwidth = NULL, adjust = 1,
                           num_points = 10000, ...) {


  if (is.null(bandwidth)) {
    bandwidth <- ks::hpi(x, deriv.order = deriv)
  }
  #we use the private version of drvkde in flowStats (copied from feature package) to avoid the tcltk dependency
  deriv_x <- flowStats:::drvkde(x = x, drv = deriv, bandwidth = adjust * bandwidth,
                                gridsize = num_points, ...)
  list(x = deriv_x$x.grid[[1]], y = deriv_x$est)
}


