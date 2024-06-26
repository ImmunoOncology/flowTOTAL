% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/runEstimateProportion.R
\name{runEstimateProportion}
\alias{runEstimateProportion}
\title{Run Proportion Estimation and Gating for Multiple Files}
\usage{
runEstimateProportion(log_file_track, panel_estimate, output, cluster = NULL)
}
\arguments{
\item{log_file_track}{The path to the log file containing file tracking information.}

\item{panel_estimate}{File path to the panel for proportion estimation.
It contains information about patterns and labels.}

\item{output}{The output directory for results.}

\item{cluster}{An optional cluster object for parallel processing (default: NULL).}
}
\description{
This function runs the proportion estimation and gating process for multiple files.
It reads a log file to identify completed files, applies the estimation and gating functions,
and outputs results to the specified directory.
}
\examples{
runEstimateProportion(
  log_file_track = "output_directory/log_file.txt",
  panel_estimate = "panel_estimate.csv",
  output = "output_directory"
)

}
\keyword{flowCore}
