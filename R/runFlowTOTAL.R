#' Run the Complete Flow Analysis Pipeline
#'
#' This function runs the complete flow analysis pipeline in sequential or parallel mode.
#' It includes preprocessing, density-based backgating, and proportion estimation steps.
#'
#' @param fcs_path The path to the folder containing FCS files.
#' @param output The output directory for results.
#' @param panel_backgating The channel for density-based backgating. Vector of channel names and signs specifying the gating operations.
#' @param panel_estimate The panel for proportion estimation. File path.
#' @param log_file_track Path to log_file_track.txt file in order to run
#' runEstimateProportion (default: NULL. It takes the one created during runDensityBackgating in output/log_file_track.txt)
#' @param cluster An optional cluster object for parallel processing (default: NULL).
#' @param steps A character vector specifying the steps to run (default: runPreprocessing,
#' runDensityBackgating, runEstimateProportion, runSEDA, runDA).
#' @param metadata Optional data frame containing additional metadata information.
#' @param marker Whether to perform marker-specific analysis (see runSEDA fn).
#' @param downsampling The type of downsampling to perform ("random", "minMaxSamplingSEACELLS", "minMaxSampling", or none) (see runSEDA fn).
#' @param k_downsampling Percentage of umber of cells to downsample for each ID (see runSEDA fn).
#' @param seed Seed for reproducibility if downsampling is used.
#' @param batch Whether to perform batch effect correction (see runSEDA fn).
#' @param response Column name in the metadata indicating the response variable (see runDA fn).
#' @param response_label Labels for the response variable levels (see runDA fn).
#' @param conda_env Conda environment to use for running SEACells.
#'
#' @keywords flowCore
#' @export
#'
#' @examples
#' library(flowCore)
#' library(flowTOTAL)
#'
#' runFlowTOTAL(
#'   fcs_path = "path/to/fcs_files",
#'   output = "output_directory",
#'   panel_backgating = "CD3+",
#'   panel_estimate = "panel_estimate.txt",
#'   cluster = my_cluster_object,
#'   steps = c("runPreprocessing", "runDensityBackgating", "runEstimateProportion")
#' )
#'
runFlowTOTAL <- function(fcs_path, output, panel_backgating, panel_estimate, log_file_track=NULL, cluster = NULL,
                         metadata=NULL, marker=FALSE, downsampling="random", k_downsampling=0.01, seed=123, batch=FALSE,
                         response=NULL, response_label=NULL, conda_env=NULL,
                         steps = c("runPreprocessing", "runDensityBackgating", "runEstimateProportion")) {

  # Create a metadata data frame to store file information
  metadata <- data.frame(filename = list.files(fcs_path, full.names = TRUE))
  metadata$file <- sapply(strsplit(metadata$filename, "/"), function(x) x[length(x)])
  metadata$filename_clean <- file.path(output, "fcs_clean", metadata$file)

  # Step 1: Preprocessing
  if ("runPreprocessing" %in% steps) {
    message("Step 1: runPreprocessing")
    runPreprocessing(metadata = metadata, output = output, report = TRUE, cluster = cluster)
  } else {
    message("Step 1: runPreprocessing - skipped")
  }

  # Step 2: Density-Based Backgating
  if ("runDensityBackgating" %in% steps) {
    message("Step 2: runDensityBackgating")
    runDensityBackgating(metadata = metadata, output = output, channel_bg = panel_backgating, cluster = cluster)
  } else {
    message("Step 2: runDensityBackgating - skipped")
  }

  # Step 3: Proportion Estimation
  if ("runEstimateProportion" %in% steps) {
    message("Step 3: runEstimateProportion")
    if(is.null(log_file_track)) log_file_track <- file.path(output, "Log_file_track.txt")
    runEstimateProportion(log_file_track, panel_estimate, output, cluster = cluster)
  } else {
    message("Step 3: runEstimateProportion - skipped")
  }

  # Step 4: Single Event Downstream Analyses
  if ("runSEDA" %in% steps) {
    message("Step 4: runSEDA")
    runSEDA(output, metadata=metadata, marker=marker, downsampling=downsampling, k_downsampling=k_downsampling, seed=seed, batch=batch, conda_env=conda_env)
  } else {
    message("Step 4: runSEDA - skipped")
  }

  # Step 5: Differential Abundance Analyses
  if ("runDA" %in% steps) {
    message("Step 5: runDA")
    runDA(output, response=response, response_label=response_label)
  } else {
    message("Step 5: runDA - skipped")
  }
}
