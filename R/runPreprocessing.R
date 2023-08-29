#' Filter Singlet Cells from Flow Cytometry Data
#'
#' This function filters singlet cells from flow cytometry data using the PeacoQC package
#' to remove doublets and the openCyto package for singlet gating. It is a double checking for singlets by looking `FSC-A` and `FSC-H` channels.
#'
#' @param fC A flow cytometry dataset (typically a flowFrame object) containing the raw data from flowCore package.
#' @param chnl A character vector of length 2 specifying the channels to use for singlet gating.
#'             The default is c("FSC-A", "FSC-H").
#'
#' @return A flow cytometry dataset with only singlet cells based on the gating results.
#'
#' @import PeacoQC
#' @import flowCore
#' @import openCyto
#'
#' @examples
#' \dontrun{
#' filtered_data <- filterSinglets(fC = my_flow_data, chnl = c("FSC-A", "FSC-H"))
#' }
#'
#' @export
filterSinglets <- function(fC, chnl = c("FSC-A", "FSC-H")) {
  # Remove doublets using PeacoQC package
  fC_PeacoQC <- PeacoQC::RemoveDoublets(fC, channel1 = chnl[1], channel2 = chnl[2])

  # Perform singlet gating using openCyto package
  gate_singlet <- openCyto:::.singletGate(fC_PeacoQC, channels = chnl)

  # Apply singlet gating to the data
  idt_singlet <- flowCore::filter(fC_PeacoQC, gate_singlet)@subSet

  # Create a new flow cytometry dataset with only singlet cells
  fC_singlet <- fC_PeacoQC
  fC_singlet@exprs <- fC_PeacoQC@exprs[idt_singlet, ]

  return(fC_singlet)
}


#' Simplify flowCore Object
#'
#' This function simplifies a flowCore object by removing empty channels and renaming
#' channels based on their descriptions. It also allows for keeping specific channels if needed.
#'
#' @param filename The path to the FCS file to be simplified.
#' @param keep A character vector specifying the names of channels to keep. Default is NULL.
#'
#' @return If `keep` is NULL, returns the path to the simplified FCS file. If `keep` is specified,
#'         returns the path to the simplified FCS file after keeping only specified channels.
#'
#' @import flowCore
#'
#' @examples
#' \dontrun{
#' simplified_file <- simplify_flowCore(filename = "path/to/input.fcs")
#' simplified_kept_file <- simplify_flowCore(filename = "path/to/input.fcs",
#' keep = c("FSC-A", "SSC-A"))
#' }
#'
#' @export
simplify_flowCore <- function(filename, keep = NULL) {

  # Read the FCS file using flowCore
  fC <- flowCore::read.FCS(filename)

  # Identify indices of shape and time channels
  parameters_name <- names(fC@parameters@data$name)
  parameters_desc <- names(fC@parameters@data$desc)
  idt <- unique(c(grep("[FS]SC-", fC@parameters@data$name), grep("Time", fC@parameters@data$name), which(is.na(fC@parameters@data$desc))))

  # Update parameter names and descriptions
  fC@parameters@data$name[-idt] <- fC@parameters@data$desc[-idt]
  names(fC@parameters@data$name)[-idt] <- parameters_name[-idt]
  names(fC@parameters@data$desc)[-idt] <- parameters_desc[-idt]

  # Update column names in the exprs matrix
  colnames(fC@exprs) <- fC@parameters@data$name

  # Keep specified parameters and write to file if needed
  if (!is.null(keep)) {
    if (all(keep %in% fC@parameters@data$name)) {
      fC <- fC[, fC@parameters@data$name %in% keep]
      flowCore::write.FCS(fC, filename)
    } else {
      return(FALSE)
    }
  } else {
    flowCore::write.FCS(fC, filename)
  }

  return(filename)
}


#' Perform Preprocessing on Flow Cytometry Data
#'
#' This function performs a series of preprocessing steps on flow cytometry data, including quality control, compensation,
#' singlet gating, and generating a preprocessing report.
#'
#' @param file Path to the input FCS file.
#' @param filename Desired name for the output FCS file (without extension).
#' @param output Directory where output files and reports will be saved.
#' @param report Logical indicating whether to generate a preprocessing report. Default is TRUE.
#'
#' @return If successful, the function performs preprocessing steps and generates output files.
#'
#' @import flowCore
#'
#' @examples
#' \dontrun{
#' doPreprocessing(file = "path/to/input.fcs", filename =
#'   "preprocessed_input.fcs", output = "output_directory")
#' }
#'
#' @export
doPreprocessing <- function(file, filename, output, report = TRUE) {

  # Create the output directory if it doesn't exist
  if (!dir.exists(output)) {
    message("Creating directory -->", output)
    dir.create(output)
  }

  # Add .fcs extension if missing
  if (!grepl(".fcs$", filename)) {
    message("Adding extension .fcs")
    filename <- paste0(filename, ".fcs")
  }

  # Read the input FCS file
  ff <- flowCore::read.FCS(file)

  # Set identifier and compensate if necessary
  flowCore::identifier(ff) <- gsub(".fcs$", "", filename)
  if ("SPILL" %in% names(ff@description)) {
    ff_comp <- flowCore::compensate(ff, spillover = flowCore::spillover(ff)$SPILL)
  } else {
    ff_comp <- ff
  }

  cat(paste0("Quality control for the file: ", filename, "\n"))

  # Perform quality control and generate mini-report
  res_QC <- tryCatch({
    flow_auto_qc_custom(ff_comp, filename = filename, ChExcludeFS = NULL, ChExcludeFM = NULL, mini_report = "Preprocessing", folder_results = paste0(output, "/resultsQC"))
  }, error = function(x) {
    res_QC <- list(
      FCS = ff_comp,
      minireport = data.frame("File" = NA, "N.initial.events" = NA, "FlowRateQC" = NA, "FlowSignalQC" = NA, "FlowMarginQC" = NA, "RemoveDoublets" = NA, "N.final.events" = NA)
    )
    return(res_QC)
  })
  ff_QC <- res_QC$FCS

  # Identify appropriate channels for singlet gating
  chnl <- c("FSC-A", "FSC-H")[c("FSC-A", "FSC-H") %in% ff_QC@parameters@data$name]

  # Apply singlet gating if necessary
  if (length(chnl) == 2) {
    ff_singlet <- filterSinglets(ff_QC, chnl = chnl)
  } else {
    ff_singlet <- ff_QC
  }

  # Update the preprocessing mini-report
  res_QC$minireport$RemoveDoublets <- nrow(ff_QC@exprs) - nrow(ff_singlet@exprs)
  res_QC$minireport$N.final.events <- nrow(ff_singlet@exprs)
  res_QC$minireport$File <- file

  # Write and simplify the output FCS file
  flowCore::write.FCS(ff_singlet, filename = file.path(output, filename))
  simplify_flowCore(file.path(output, filename))

  # Generate and save the preprocessing report
  if (report) {
    reporte_filename <- file.path(output, "/resultsQC/Preprocessing.txt")
    if (!dir.exists(file.path(output, "/resultsQC"))) dir.create(file.path(output, "/resultsQC"))
    write.table(
      x = res_QC$minireport,
      file = reporte_filename,
      col.names = !file.exists(reporte_filename),
      sep = "\t",
      row.names = FALSE,
      append = file.exists(reporte_filename),
      quote = FALSE
    )
  }
}

#' Run Preprocessing Function
#'
#' This function runs the preprocessing analysis for raw FCS files. It performs compensation and quality control,
#' including removing doublets and anomalies.
#'
#' @param metadata A data frame with a column named "filename" containing paths to the raw FCS files.
#' @param output Path to the location where the cleaned FCS files will be saved.
#' @param report Logical indicating whether to report the number of anomalies and doublets. Default is TRUE.
#' @param cluster Number of parallel workers to be used for processing. Default is NULL.
#' @param log_file Path to the log file where error messages will be written.
#'
#' @return This function runs the preprocessing analysis and generates cleaned FCS files and reports.
#'
#' @import flowCore
#' @import doParallel
#' @import foreach
#' @export
runPreprocessing <- function(metadata, output, report = TRUE, cluster = NULL, log_file = "Log_file_error_raw2clean.txt") {

  log_file_error <- function(messages) {
    sink(log_file, append = TRUE)
    cat(paste(messages, "\n"), sep = "\n")
    sink()
  }

  # Set the output directory
  output <- file.path(output, "fcs_clean")

  # Check if "filename" column exists in metadata
  if (!"filename" %in% colnames(metadata))
    stop("Error in metadata")

  # Define the preprocessing function
  preprocess_function <- function(file) {
    filename <- basename(file)
    tryCatch({
      doPreprocessing(file, filename, output, report = report)
    }, error = function(e) {
      log_file_error(paste(e, "File:", filename))
    })
  }

  # Perform preprocessing sequentially or in parallel
  if (is.null(cluster)) {
    for (i in 1:nrow(metadata)) {
      preprocess_function(metadata$filename[i])
    }
  } else {
    requireNamespace('foreach')
    foreach(i = 1:nrow(metadata), .packages = c("flowCore")) %dopar% {
      preprocess_function(metadata$filename[i])
    }
    stopCluster(cl)
  }
}
