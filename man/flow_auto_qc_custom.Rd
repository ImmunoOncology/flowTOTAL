% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/helperFunctions.R
\name{flow_auto_qc_custom}
\alias{flow_auto_qc_custom}
\title{Custom Flow Auto QC from flowAI}
\usage{
flow_auto_qc_custom(
  fcsfiles,
  filename = "V1",
  timeCh = NULL,
  second_fractionFR = 0.1,
  alphaFR = 0.01,
  decompFR = TRUE,
  ChExcludeFS = c("FSC", "SSC"),
  outlier_binsFS = FALSE,
  pen_valueFS = 500,
  max_cptFS = 3,
  ChExcludeFM = c("FSC", "SSC"),
  sideFM = "both",
  neg_valuesFM = 1,
  mini_report = "QC_report",
  folder_results = "resultsQC"
)
}
\arguments{
\item{fcsfiles}{A flowSet containing raw FCS files to be analyzed.}

\item{filename}{A character string to identify the analysis.}

\item{timeCh}{Name of the time channel in the FCS files.}

\item{second_fractionFR}{The second fraction for flow rate calculation. Default is 0.1.}

\item{alphaFR}{The significance level for flow rate check. Default is 0.01.}

\item{decompFR}{Logical indicating whether to use decomposition for flow rate check. Default is TRUE.}

\item{ChExcludeFS}{Channels to be excluded from flow signal QC. Default is c("FSC", "SSC").}

\item{outlier_binsFS}{Logical indicating whether to identify outlier bins for flow signal QC. Default is FALSE.}

\item{pen_valueFS}{Penalty value for flow signal QC. Default is 500.}

\item{max_cptFS}{Maximum number of change points for flow signal QC. Default is 3.}

\item{ChExcludeFM}{Channels to be excluded from flow margin QC. Default is c("FSC", "SSC").}

\item{sideFM}{Side for flow margin QC. Default is "both".}

\item{neg_valuesFM}{Number of negative values allowed for flow margin QC. Default is 1.}

\item{mini_report}{Filename for the mini report.}

\item{folder_results}{Path to the folder for saving results.}
}
\value{
A list containing cleaned flow cytometry data (flowFrame) and a summary report.
}
\description{
This function performs quality control for flow cytometry data using custom criteria based on \code{flow_auto_qc} function from flowAI package.
It calculates flow rates, signal quality, and flow margins to assess the data quality.
}
\details{
Monaco G, Chen H (2022). flowAI: automatic and interactive quality control for flow cytometry data. R package version 1.30.0, 10.18129/B9.bioc.flowAI.
}
