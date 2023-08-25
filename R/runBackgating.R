#' Check FCS
#'
#' Remove empty channels and update channel descriptions.
#'
#' This function takes an FCS (Flow Cytometry Standard) file and performs checks and updates
#' on its channel descriptions. It removes empty channels and ensures that channel descriptions
#' are appropriately set.
#'
#' @param ff An object of class 'flowFrame' representing the FCS file.
#'
#' @keywords FCS
#'
#' @return An updated flowFrame object with corrected channel descriptions.
#'
#' @examples
#' # Load the required libraries
#' library(flowCore)
#'
#' # Load an FCS file using flowCore
#' fcs_file <- read.FCS("path/to/your/fcsfile.fcs")
#'
#' # Perform FCS channel check and update
#' updated_fcs <- check_FCS(fcs_file)
#'
#' # Now you can use the updated_fcs with corrected channel descriptions
#'
#' @import flowCore
#' @export
check_FCS <- function(ff) {
  # Find channels with names starting with "SC-"
  idt <- grep("SC-", ff@parameters@data$name, fixed = TRUE)

  # Update missing channel descriptions with channel names
  if (any(is.na(ff@parameters@data$desc[idt]))) {
    ff@parameters@data$desc[idt] <- ff@parameters@data$name[idt]
  }

  # Find channels with spaces in their descriptions
  idt <- grep(" ", ff@parameters@data$desc, fixed = TRUE)

  # Update channel descriptions by removing text after space
  if (length(idt) > 0) {
    ff@parameters@data$desc[idt] <- unlist(lapply(strsplit(ff@parameters@data$desc[idt], " "), `[`, 1))
  }

  # Return the updated flowFrame object
  return(ff)
}


##' Get Contour
#'
#' Generate and plot contour lines of the density estimate for given FSC-A and SSC-A channels.
#'
#' This function calculates the kernel density estimate of the given FSC-A and SSC-A channels
#' and generates contour lines from the density estimate. It iteratively refines the contour lines
#' to identify dense regions and plots them using ggplot2.
#'
#' @param ff An object of class 'flowFrame' representing the flow cytometry data.
#'
#' @importFrom ks kde Hns
#' @importFrom sp point.in.polygon
#' @import ggplotify
#' @import ggplot2
#'
#' @return A list of contour lines containing x and y coordinates.
#'
#' @examples
#' # Load the required libraries
#' library(flowCore)
#' library(ggplot2)
#' library(ggplotify)
#'
#' # Load an FCS file using flowCore
#' fcs_file <- read.FCS("path/to/your/fcsfile.fcs")
#'
#' # Generate and plot contour lines
#' contour_lines <- get_contour(fcs_file)
#'
#' @export
get_contour <- function(ff) {
  # Helper function to create a sequence of consecutive numbers
  mkseq <- function(my_vector) {
    diffs <- diff(my_vector)
    unique_vals <- c(my_vector[1], my_vector[diffs != 1] + 1)
    return(unique_vals)
  }

  # Extract FSC-A and SSC-A channels from flowFrame object
  samp <- ff@exprs[, c("FSC-A", "SSC-A")]

  # Calculate kernel density estimate
  Hns <- ks::Hns(x = samp, deriv.order = 0)
  kdde_0 <- ks::kde(x = samp, H = diag(diag(Hns)))

  # Determine suitable percentile levels for contour lines
  dummy <- unlist(lapply(1:99, function(pct_i) length(with(kdde_0, contourLines(x = eval.points[[1]], y = eval.points[[2]],
                                                                                z = estimate, levels = cont[paste0(pct_i, "%")][[1]])))))
  pcts <- which(dummy == max(dummy))

  # Identify single percentile sequence
  if (length(mkseq(pcts)) == 1) {
    dummy_idt <- dummy
    dummy_idt[pcts] <- 0
    if (max(dummy_idt) != 1)
      pcts <- c(pcts, which(dummy_idt == max(dummy_idt)))
  }

  pcts <- mkseq(pcts)

  # Generate initial contour lines
  contour.dummy <- do.call("c",
                           lapply(pcts,
                                  function(pct_i){
                                    with(kdde_0,contourLines(x=eval.points[[1]],y=eval.points[[2]], z=estimate,levels=cont[paste0(pct_i, "%")]))
                                    }
                                  ))
  # Filter out sparse contour lines
  point.dummy <- lapply(1:length(contour.dummy), function(x) which(sp::point.in.polygon(point.x = ff@exprs[, "FSC-A"], point.y = ff@exprs[, "SSC-A"], pol.x = contour.dummy[[x]]$x, pol.y = contour.dummy[[x]]$y) == 1))
  contour.dummy <- contour.dummy[lengths(point.dummy) > 250]
  point.dummy <- lapply(1:length(contour.dummy), function(x) which(sp::point.in.polygon(point.x = ff@exprs[, "FSC-A"], point.y = ff@exprs[, "SSC-A"], pol.x = contour.dummy[[x]]$x, pol.y = contour.dummy[[x]]$y) == 1))

  # Initialize variables
  idt <- 1:length(point.dummy)
  res <- unlist(lapply(idt, function(idt2) {
    sum(unlist(lapply(idt[!idt %in% idt2], function(idt3) ifelse(length(intersect(point.dummy[[idt2]], point.dummy[[idt3]])) > 0, 1, 0))))
  }))

  # Iterate to refine contour lines
  cnt_dummy <- 2
  while (length(res) == 1 & cnt_dummy <= length(which(dummy == max(dummy)))) {
    pcts <- which(dummy == max(dummy))[cnt_dummy]
    cnt_dummy <- cnt_dummy + 1

    contour.dummy <- do.call("c", lapply(pcts, function(pct_i) with(kdde_0, contourLines(x = eval.points[[1]], y = eval.points[[2]],
                                                                                         z = estimate, levels = cont[paste0(pct_i, "%")]))))
  contour.dummy <- contour.dummy[lengths(point.dummy) > 250]
  point.dummy <- lapply(1:length(contour.dummy), function(x) which(sp::point.in.polygon(point.x = ff@exprs[, "FSC-A"], point.y = ff@exprs[, "SSC-A"], pol.x = contour.dummy[[x]]$x, pol.y = contour.dummy[[x]]$y) == 1))

  idt <- 1:length(point.dummy)
  res <- unlist(lapply(idt, function(idt2) {
    sum(unlist(lapply(idt[!idt %in% idt2], function(idt3) ifelse(length(intersect(point.dummy[[idt2]], point.dummy[[idt3]])) > 0, 1, 0))))
  }))
  }

  # Refine contour lines iteratively
  while (any(res > 0)) {
    contour.dummy <- contour.dummy[-which.max(res)]
    point.dummy <- lapply(1:length(contour.dummy), function(x) which(sp::point.in.polygon(point.x = ff@exprs[, "FSC-A"], point.y = ff@exprs[, "SSC-A"], pol.x = contour.dummy[[x]]$x, pol.y = contour.dummy[[x]]$y) == 1))

    idt <- 1:length(point.dummy)
    res <- unlist(lapply(idt, function(idt2) {
      sum(unlist(lapply(idt[!idt %in% idt2], function(idt3) ifelse(length(intersect(point.dummy[[idt2]], point.dummy[[idt3]])) > 0, 1, 0))))
    }))
  }

  # Create data frame for contour lines and plot using ggplot2
  contour.plot <- do.call('rbind', lapply(contour.dummy, function(df) data.frame(x = df[["x"]], y = df[["y"]])))
  ggplotify::as.ggplot(function(x) {
    plot(kdde_0)
    points(contour.plot$x, contour.plot$y, col = "orange", cex = 0.3)
  })

  # Return the refined contour lines
  return(contour.dummy)
}

#' Perform density-based gating on flow cytometry data
#'
#' This function performs density-based gating on flow cytometry data using various parameters.
#'
#' @param ff FlowFrame object containing the flow cytometry data.
#' @param filename Name of the input file.
#' @param output.dir Directory where output files will be saved.
#' @param chnl Vector of channel names for scatter plots (default: c("FSC-A", "SSC-A")).
#' @param channel_bg Vector of channel names for background gating.
#' @param logicle_chnls Vector of channel names for logicle transformation (if applicable).
#' @param sd.max_it Standard deviation factor for adjusting gating percentages.
#' @param min.pct_it Minimum gating percentage for selection.
#' @param target.fsc Target FSC value for gating adjustment.
#' @param target.ssc Target SSC value for gating adjustment.
#' @param min.ff_subset Minimum number of events to keep in the subset.
#'
#' @return Number of selected events after gating.
#'
#' @import ggcyto
#' @import ggplot2
#' @importFrom flowCore filter write.FCS
#' @importFrom stats density
#' @importFrom sp point.in.polygon
#' @importFrom openCyto .find_peaks gate_tail_custom gate_flowclust_2d
#'
#' @export
#'
#' @examples
#' # Load required packages
#' library(ggcyto)
#' library(flowCore)
#'
#' # Load flow cytometry data (replace with actual data)
#' data_file <- system.file("extdata", "sample_data.fcs", package="flowCore")
#' ff <- read.FCS(data_file)
#'
#' # Define channel names and background gating channels
#' channel_bg <- c("CD3-FITC", "CD4-PE", "CD8-APC")
#'
#' # Perform density-based gating
#' doDensityBackgating(ff, "sample_data.fcs", output.dir = "output_directory", channel_bg = channel_bg)
#'
#' @seealso
#' \code{\link{ggcyto::autoplot}}
#' \code{\link{ggplot2::ggsave}}
#' \code{\link{flowCore::filter}}
#' \code{\link{flowCore::write.FCS}}
#' \code{\link{openCyto::gate_tail_custom}}
#' \code{\link{openCyto::gate_flowclust_2d}}
#'
#' @keywords flow cytometry gating density-based
doDensityBackgating <- function(ff, filename, output.dir, chnl = c("FSC-A", "SSC-A"), channel_bg, logicle_chnls=NULL, sd.max_it=0.75, min.pct_it=0.01, target.fsc=50000, target.ssc=17500, min.ff_subset=250) {

  # Make a copy of the original data
  ff.raw <- ff

  # Define comparison functions for gating operations
  `[+]` <- function(a, b) { a > b }
  `[-]` <- function(a, b) { a < b }

  # Extract panel channel names and signs from channel_bg
  panel_channel <- sapply(channel_bg, function(x) substr(x, 1, nchar(x) - 1))
  channel_sign <- paste0("[", sapply(channel_bg, function(x) substr(x, nchar(x), nchar(x))))

  # Initialize lists for storing gating results and subsets
  gate_chnl <- list()
  keep <- rep(TRUE, nrow(ff))

  # Iterate through each panel channel
  for (j in 1:length(panel_channel)) {
    # Find peaks and calculate densities
    n_peaks <- openCyto:::.find_peaks(ff@exprs[keep, panel_channel[j]], adjust = 0.5)
    n_peaks <- n_peaks[order(n_peaks$x, decreasing = FALSE), ]

    n_peaks_dummy <- n_peaks[n_peaks$y > max(n_peaks$y) * 0.01, ]
    n_peaks_dummy <- n_peaks_dummy[order(n_peaks_dummy$x, decreasing = FALSE), ]

    dummy <- density(ff@exprs[keep, panel_channel[j]], adjust = 0.5)

    # Refine peak positions based on density
    while (nrow(n_peaks_dummy) > 1 && all(sapply(dummy$y[dummy$x > n_peaks_dummy$x[1] & dummy$x < n_peaks_dummy$x[2]][-1], function(x) dummy$y[dummy$x > n_peaks_dummy$x[1] & dummy$x < n_peaks_dummy$x[2]][1] - x) < 0)) {
      n_peaks_dummy <- n_peaks_dummy[2:nrow(n_peaks_dummy), ]
    }

    # Determine the reference peak
    if (n_peaks$x[1] != n_peaks_dummy$x[1]) {
      ref <- which(n_peaks$x == n_peaks_dummy$x[1])
    } else {
      ref <- 1
    }

    # Apply custom gating to the channel
    gate_chnl[[j]] <- gate_tail_custom(ff[keep, ], panel_channel[j], adjust = 0.5, ref_peak = ref, auto_tol = TRUE, num_peaks = nrow(n_peaks))

    # Handle the case of a single peak
    if (nrow(n_peaks) == 1) {
      n_peaks <- openCyto:::.find_peaks(ff@exprs[, panel_channel[j]], adjust = 1)
      gate_chnl[[j]] <- gate_tail_custom(ff, panel_channel[j], adjust = 1, ref_peak = 1, auto_tol = TRUE, num_peaks = nrow(n_peaks))
    }

    # Update the subset using the gating results
    keep <- keep & do.call(channel_sign[j], args = list(ff@exprs[, panel_channel[j]], gate_chnl[[j]]@min))
  }

  # Subset the original data based on the gating results
  ff_subset <- ff[keep, ]

  # Create density plots and gating visualization plots
  p1 <- ggcyto::autoplot(ff, "FSC-A", "SSC-A", bins = 100)
  p_gate <- ggcyto::autoplot(ff, panel_channel[1], panel_channel[2], bins = 100)
  for (gate in gate_chnl) {
    p_gate <- p_gate + ggcyto::geom_gate(gate)
  }

  # Append the subset data to the original data
  ff@exprs <- rbind(ff@exprs, ff_subset@exprs)

  # Calculate and store contours
  contour.pct <- get_contour(ff)

  # Initialize lists and plots for further analysis
  my_gates <- list()
  p2 <- ggcyto::autoplot(ff.raw, "FSC-A", "SSC-A", bins = 100)
  gate_keep <- list()
  big.cells <- list()
  ssc.cells <- list()
  fsc.cells <- list()
  pct_it_cnt <- 1

  # Loop through contour percentages
  for (pct_it in 1:length(contour.pct)) {
    inner <- sp::point.in.polygon(point.x = ff@exprs[, "FSC-A"], point.y = ff@exprs[, "SSC-A"], pol.x = contour.pct[[pct_it]]$x, pol.y = contour.pct[[pct_it]]$y)

    # Find peaks within the inner polygon
    K <- openCyto:::.find_peaks(ff@exprs[inner == 1, "FSC-A"], adjust = 2)
    K <- K[K$y > K$y * 0.1, ]

    # Loop through the peaks
    for (j in 1:nrow(K)) {
      # Attempt to gate using FlowClust
      gate_j <- tryCatch({
        suppressMessages(openCyto::gate_flowclust_2d(ff[inner == 1, ], yChannel = "SSC-A", xChannel = "FSC-A", K = nrow(K), target = c(K$x[j], median(ff@exprs[inner == 1, "SSC-A"]))))
      }, error = function(x) {
        return(NA)
      })

      # If gating is successful, update plots and store data
      if (!is.na(gate_j)) {
        my_gates[[pct_it_cnt]] <- gate_j
        p2 <- p2 + ggcyto::geom_gate(my_gates[[pct_it_cnt]])

        keep <- flowCore::filter(ff.raw, my_gates[[pct_it_cnt]])@subSet
        big.cells[[pct_it_cnt]] <- mean(ff.raw@exprs[keep, c("SSC-A")])
        ssc.cells[[pct_it_cnt]] <- mean(ff.raw@exprs[keep, c("SSC-A")])
        fsc.cells[[pct_it_cnt]] <- mean(ff.raw@exprs[keep, c("FSC-A")])

        # Apply gating to all channels
        keep_gate <- rep(TRUE, nrow(ff.raw))
        for (k in 1:length(panel_channel)) {
          keep_gate <- keep_gate & do.call(channel_sign[k], args = list(ff.raw@exprs[, panel_channel[k]], gate_chnl[[k]]@min))
        }

        # Convert gating results to a factor for storage
        keep_chl <- factor(ifelse(keep_gate, 1, 0), levels = c("0", "1"))
        gate_keep[[pct_it_cnt]] <- keep_chl[keep]

        pct_it_cnt <- pct_it_cnt + 1
      }
    }
  }

  # Remove small cells based on FSC and SSC values
  too.small <- which(unlist(fsc.cells) < 10000 & unlist(ssc.cells) < 10000)
  if (length(too.small) > 0 & length(too.small) != length(fsc.cells)) {
    gate_keep <- gate_keep[-too.small]
    fsc.cells <- fsc.cells[-too.small]
    ssc.cells <- ssc.cells[-too.small]
    my_gates <- my_gates[-too.small]
  } else if (length(too.small) > 0 & length(too.small) == length(fsc.cells)) {
    too.small <- which(unlist(fsc.cells) < 250 & unlist(ssc.cells) < 250)
    if (length(too.small) > 0 & length(too.small) != length(fsc.cells)) {
      gate_keep <- gate_keep[-too.small]
      fsc.cells <- fsc.cells[-too.small]
      ssc.cells <- ssc.cells[-too.small]
      my_gates <- my_gates[-too.small]
    }
  }

  # Calculate maximum percentage and gating percentages
  max_it <- max(unlist(lapply(gate_keep, function(x) (table(x) / length(x))[2])))
  pct_it <- unlist(lapply(gate_keep, function(x) (table(x) / length(x))[2]))

  # Adjust gating based on percentages and thresholds
  if (all(pct_it < min.pct_it) || nrow(ff_subset) < min.ff_subset) {
    pct_it <- which.min(abs(unlist(fsc.cells) - target.fsc) + abs(unlist(ssc.cells) - target.ssc))
  } else if (any(pct_it > 0.5)) {
    pct_it <- which(pct_it > 0.4)
  } else {
    pct_it <- which(pct_it > (max_it * sd.max_it) & pct_it > min.pct_it)
  }

  # Create a final subset of events based on gating results
  idt_final <- rep(FALSE, nrow(ff.raw))
  for (pct_it_gate in pct_it) {
    p2 <- p2 + ggcyto::geom_gate(my_gates[[pct_it_gate]], col = "blue")
    idt_final <- idt_final | flowCore::filter(ff.raw, my_gates[[pct_it_gate]])@subSet
  }

  # Create a final density plot of the selected events
  p_final <- ggcyto::autoplot(ff.raw[idt_final, ], panel_channel[1], panel_channel[2], bins = 100)

  # Arrange plots and save as PDF
  p <- ggpubr::ggarrange(plotlist = list(ggcyto::as.ggplot(p1), ggcyto::as.ggplot(p_gate), ggcyto::as.ggplot(p2), ggcyto::as.ggplot(p_final)), nrow = 2, ncol = 2)
  ggplot2::ggsave(paste0(output.dir, "/PDF/", gsub(".fcs", "", filename, fixed = TRUE), ".pdf"), plot = p, device = "pdf", width = 15, height = 12)

  # Write the selected events to an FCS file
  flowCore::write.FCS(ff.raw[idt_final, ], filename = paste0(output.dir, "/FCS/", filename))

  # Return the count of selected events
  return(sum(idt_final))
}

#' Run Density-Based Gating on Flow Cytometry Data
#'
#' This function performs density-based gating on a batch of flow cytometry data
#' files specified in the metadata. It applies density-based gating to each file,
#' saves the results, and logs in the output directory.
#'
#' @param metadata A data frame containing information about the files to be processed.
#'   It should have at least one column named 'filename_clean' containing the paths to the FCS files.
#'   Example: data.frame(filename_clean = c("file1.fcs", "file2.fcs"))
#' @param output Output directory path where the results and logs will be stored.
#' @param chnl Vector of channel names for which density-based gating is performed.
#'   Default: c("FSC-A", "SSC-A")
#' @param channel_bg Vector of channel names and signs specifying the gating operations.
#' @param logicle_chnls Vector of channel names that require logicle transformation.
#' @param sd.max_it Standard deviation threshold for peak detection.
#' @param min.pct_it Minimum percentage of gated events required for iteration.
#' @param target.fsc Target value for FSC-A gating.
#' @param target.ssc Target value for SSC-A gating.
#' @param min.ff_subset Minimum number of events required for the final subset.
#' @param log_file File name for the error log.
#' @param log_file_traditional File name for the traditional error log.
#' @param track_file File name for the tracking log.
#' @param cluster Number of cluster nodes for parallel processing. If NULL, no parallel processing is used.
#'
#' @export
#' @examples
#' runDensityBackgating(metadata = my_metadata, output = "results", channel_bg = c("FL1-H", "FL2-H"),
#'                      sd.max_it = 0.5, target.fsc = 50000, target.ssc = 20000)
#'
#' @import flowCore
#' @import foreach
#' @import doParallel
runDensityBackgating <- function(metadata, output, chnl = c("FSC-A", "SSC-A"), channel_bg, logicle_chnls=NULL, sd.max_it=0.75,
                                 min.pct_it=0.01, target.fsc=50000, target.ssc=17500, min.ff_subset=250,
                                 log_file="Log_file_error.txt", log_file_traditional="Log_file_error_traditional.txt", track_file="Log_file_track.txt", cluster=NULL){

  # Create output directories if they do not exist
  output.dir <- file.path(output, "Backgating")
  if(!dir.exists(output.dir)) dir.create(output.dir, recursive = TRUE)
  if(!dir.exists(file.path(output.dir, "PDF"))) dir.create(file.path(output.dir, "PDF"))
  if(!dir.exists(file.path(output.dir, "FCS"))) dir.create(file.path(output.dir, "FCS"))

  # Define error logging functions
  log_file <- file.path(output, log_file)
  log_file_error <- function(messages){
    sink(log_file, append = TRUE)
    cat(paste(messages, collapse = "\n"), "\n")
    sink()
  }

  log_file_traditional <- file.path(output, log_file_traditional)
  log_file_error_traditional <- function(messages){
    sink(log_file_traditional, append = TRUE)
    cat(paste(messages, collapse = "\n"), "\n")
    sink()
  }

  # Define tracking log function
  track_file <- file.path(output, track_file)
  log_file_track <- function(messages){
    sink(track_file, append = TRUE)
    cat(paste(messages, collapse = "\n"), "\n")
    sink()
  }


  # Define a function to process a single file
  process_file <- function(metadata, i, output.dir, channel_bg, logicle_chnls, sd.max_it, min.pct_it, target.fsc, target.ssc, min.ff_subset) {
    filename_clean <- metadata$filename_clean[i]
    filename <- unlist(strsplit(filename_clean, "/"))[length(strsplit(filename_clean, "/"))]

    tryCatch({
      ff <- flowCore::read.FCS(filename_clean)
      res_bg <- doDensityBackgating(ff, filename = filename, output.dir = output.dir, channel_bg = channel_bg, logicle_chnls = logicle_chnls,
                                    sd.max_it=sd.max_it, min.pct_it=min.pct_it, target.fsc=target.fsc, target.ssc=target.ssc, min.ff_subset=min.ff_subset)
      log_file_track(paste(filename_clean, filename, file.path(output.dir, "FCS", filename), nrow(ff@exprs), res_bg, "Completed", i, sep = "\t"))
    },
    error=function(e) {
      log_file_error(paste(e, "Iter-->", i, "\n", "File: ", metadata$filename[i]))
      log_file_track(paste(filename_clean, filename, "\t", NA, NA, "Error", i, sep = "\t"))
    })
  }

  # Set up parallel processing
  if(!is.null(cluster)){
    library(doParallel)
    foreach(i = 1:nrow(metadata), .packages = c("flowCore")) %dopar% {
      process_file(metadata, i, output.dir, channel_bg, logicle_chnls, sd.max_it, min.pct_it, target.fsc, target.ssc, min.ff_subset)
    }
    stopCluster(cl)
  } else {
    for(i in 1:nrow(metadata)){
      process_file(metadata, i, output.dir, channel_bg, logicle_chnls, sd.max_it, min.pct_it, target.fsc, target.ssc, min.ff_subset)
    }
  }
}
