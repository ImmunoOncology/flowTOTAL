% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/runFlowTOTAL.R
\name{runFlowTOTAL}
\alias{runFlowTOTAL}
\title{Run the Complete Flow Analysis Pipeline}
\usage{
runFlowTOTAL(
  fcs_path,
  output,
  panel_backgating,
  panel_estimate,
  log_file_track = NULL,
  cluster = NULL,
  metadata = NULL,
  marker = FALSE,
  downsampling = "random",
  k_downsampling = 0.01,
  seed = 123,
  batch = FALSE,
  response = NULL,
  response_label = NULL,
  conda_env = NULL,
  steps = c("runPreprocessing", "runDensityBackgating", "runEstimateProportion")
)
}
\arguments{
\item{fcs_path}{The path to the folder containing FCS files.}

\item{output}{The output directory for results.}

\item{panel_backgating}{The channel for density-based backgating. Vector of channel names and signs specifying the gating operations.}

\item{panel_estimate}{The panel for proportion estimation. File path.}

\item{log_file_track}{Path to log_file_track.txt file in order to run
runEstimateProportion (default: NULL. It takes the one created during runDensityBackgating in output/log_file_track.txt)}

\item{cluster}{An optional cluster object for parallel processing (default: NULL).}

\item{metadata}{Optional data frame containing additional metadata information.}

\item{marker}{Whether to perform marker-specific analysis (see runSEDA fn).}

\item{downsampling}{The type of downsampling to perform ("random", "minMaxSamplingSEACELLS", "minMaxSampling", or none) (see runSEDA fn).}

\item{k_downsampling}{Percentage of umber of cells to downsample for each ID (see runSEDA fn).}

\item{seed}{Seed for reproducibility if downsampling is used.}

\item{batch}{Whether to perform batch effect correction (see runSEDA fn).}

\item{response}{Column name in the metadata indicating the response variable (see runDA fn).}

\item{response_label}{Labels for the response variable levels (see runDA fn).}

\item{conda_env}{Conda environment to use for running SEACells.}

\item{steps}{A character vector specifying the steps to run (default: runPreprocessing,
runDensityBackgating, runEstimateProportion, runSEDA, runDA).}
}
\description{
This function runs the complete flow analysis pipeline in sequential or parallel mode.
It includes preprocessing, density-based backgating, and proportion estimation steps.
}
\examples{
library(flowCore)
library(flowTOTAL)

runFlowTOTAL(
  fcs_path = "path/to/fcs_files",
  output = "output_directory",
  panel_backgating = "CD3+",
  panel_estimate = "panel_estimate.txt",
  cluster = my_cluster_object,
  steps = c("runPreprocessing", "runDensityBackgating", "runEstimateProportion")
)

}
\keyword{flowCore}
