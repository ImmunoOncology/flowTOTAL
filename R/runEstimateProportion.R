#' Create ggcyto plot
#'
#' Create a ggcyto plot based on two given channels (or one and by default \code{SSC-A} will be used as the second one). The ggcyto
#' could include gates and logicle scale.
#'
#' @param fC flowCore object.
#' @param channels Vector of length 2 specifying the channels to be plotted.
#'                 If only one channel is provided, \code{SSC-A} will be used as the second channel.
#' @param gates List of gates to be included in the ggcyto plot.
#' @param logicle_chnls Vector of channel names to be transformed using logicle scale.
#' @param main Main title for the plot.
#' @param legend Logical, whether to show the legend in the plot.
#'
#' @return A ggplot object.
#'
#' @keywords plot
#' @export
#'
#' @examples
#' do_ggcyto(fC = your_flowCore_object, channels = c("FL1-A", "FL2-A"))
#'
do_ggcyto <- function(fC, channels, gates = NULL, logicle_chnls = NULL, main = NULL, legend = FALSE) {

  # Apply logicle transformation to specified channels
  if (length(logicle_chnls) > 0) {
    lgcl <- flowCore::logicleTransform()
    trans <- flowCore::transformList(logicle_chnls, lgcl)
    fC_plot <- ggcyto::transform(fC, trans)

    # Rescale gates if provided
    if (!is.null(gates)) {
      for (gate_idt in logicle_chnls) {
        logicle_chnls_tr <- logicle_chnls[logicle_chnls %in% names(gates[[gate_idt]]@min)]
        gates[[gate_idt]] <- ggcyto::rescale_gate(gates[[gate_idt]], lgcl, logicle_chnls_tr)
      }
    }

    # Rescale join gate if necessary
    if ("join" %in% names(gates)) {
      if (any(logicle_chnls %in% names(gates[["join"]]@mean)))
        gates[["join"]] <- ggcyto::rescale_gate(gates[["join"]], lgcl, logicle_chnls)
    }
  } else {
    fC_plot <- fC
  }

  # If only one channel provided, add SSC-A as the second channel
  if (length(channels) == 1) {
    channels <- c(channels, "SSC-A")
  }

  # Create ggcyto plot
  p <- ggcyto::autoplot(fC_plot, x = channels[1], y = channels[2], bins = 100)

  # Add gates to the plot
  if (!is.null(gates)) {
    for (gate in gates) {
      p <- p + ggcyto::geom_gate(gate)
    }
  }

  # Customize plot appearance
  p <- p + ggplot2::theme_bw() +
    ggplot2::theme(strip.text.x = ggplot2::element_text(size = 0, colour = "orange", angle = 90))

  # Labelling
  p <- ggcyto::as.ggplot(p)
  labels <- fC@parameters@data$desc[match(channels, fC@parameters@data$name)]
  labels[is.na(labels)] <- channels[is.na(labels)]

  for (i in 1:length(labels)) {
    if (labels[i] != channels[i]) {
      labels[i] <- paste0(labels[i], " (", channels[i], ")")
    }
  }

  labels[match(logicle_chnls, channels)] <- paste0(labels[match(logicle_chnls, channels)], " - logicleTransform")

  p <- p + ggplot2::xlab(labels[1]) + ggplot2::ylab(labels[2])

  # Add main title and adjust legend
  if (!is.null(main))
    p <- p + ggplot2::ggtitle(label = main) +
    ggplot2::theme(plot.title = ggplot2::element_text(size = 15, face = "bold", hjust = 0.5))

  if (!legend) {
    p <- p + ggplot2::theme(legend.position = "none")
  }

  return(p)
}


#' Apply openCyto gate
#'
#' Apply openCyto gate to one or two channels. If enough information about the population is kept and the join parameter is active,
#' a target gate_flowclust_2d is applied to increase resolution.
#'
#' @param ff A flowFrame object to be gated.
#' @param channel A character vector specifying the channel(s) to be used during gating.
#' @param cutpoint_min A numeric value that sets a minimum threshold for the cutpoint.
#'                    If a value is provided, any cutpoint below this value will be set to the given minimum value.
#'                    If NULL (default), there is no minimum cutpoint value.
#' @param cutpoint_max A numeric value that sets a maximum threshold for the cutpoint.
#'                    If a value is provided, any cutpoint above this value will be set to the given maximum value.
#'                    If NULL (default), there is no maximum cutpoint value.
#' @param return_plot A logical value indicating whether the result should include the plot with the gates.
#' @param join A logical value indicating whether the result should be improved by determining a target population through the given channels.
#' @param traditional A logical value indicating whether to use a traditional gating approach.
#'
#' @return A list containing the following components:
#'   \item{plot}{A ggplot object showing the plot with applied gates (if return_plot is TRUE).}
#'   \item{idt}{A logical vector indicating gated events.}
#'   \item{channel}{The channels used for gating.}
#'   \item{gate_bg}{A list of gates applied to each channel.}
#'
#' @keywords flowCore
#' @export
#'
#' @examples
#' apply_gate(ff = your_flowFrame_object, channel = c("FL1-A", "FL2-A"))
#'
apply_gate <- function(ff, channel, cutpoint_min, cutpoint_max, return_plot = TRUE, join=T, traditional=F){


  # Identify the direction of gating for each channel
  channel_sign <- paste0("[", sapply(channel, function(x) substr(x, nchar(x), nchar(x))))
  channel <-  sapply(channel, function(x) substr(x, 1, nchar(x)-1))

  # Handle channel names with '^' or '$'
  beggining <- sapply(channel, function(x) substring(x, nchar(x), nchar(x)))=="^"
  ending <- sapply(channel, function(x) substring(x, nchar(x), nchar(x)))=="$"
  channel[beggining] <-  sapply(channel[beggining], function(x) substr(x, 1, nchar(x)-1))
  channel[ending] <-  sapply(channel[ending], function(x) substr(x, 1, nchar(x)-1))

  names(channel) <- channel
  channel <- unlist(channel)

  # Define comparison functions for gating
  `[+` <- function(a, b){ a>b}
  `[-` <- function(a, b){ a<b}

  # Define adjustment functions for cutpoints
  `[+]` <- function(a){ min(a)*1.05 }
  `[-]` <- function(a){ min(a)*1.05 }


  # Define functions to estimate cutpoints
  `[+[` <- function(a){
    peak_a <- openCyto:::.find_peaks(a, adjust = 1)
    peak_a <- peak_a[peak_a$y > max(peak_a$y)*0.01, ]
    return(peak_a$x[which.max(peak_a$x)])
  }

  `[-[` <- function(a){
    peak_a <- openCyto:::.find_peaks(a, adjust = 1)
    peak_a <- peak_a[peak_a$y > max(peak_a$y)*0.01, ]
    return(peak_a$x[which.max(peak_a$x)])
  }


  idt <- rep(TRUE, nrow(ff@exprs))
  idt_list <- list()
  gate_bg <- list()

  # Estimate the number of peaks in each channel
  num_peaks <- lapply(channel, function(x){
    k_peak <- openCyto:::.find_peaks(ff@exprs[idt, x], adjust = 2)
    k_peak <- k_peak[k_peak$y > max(k_peak$y)*0.05, ]
    return(nrow(k_peak))
  })

  first <- T
  max_iter <- 3
  iter <- 1

  # Iterate through channels and apply gating
  while((any(unlist(num_peaks)>1) | first) & iter <= max_iter & (sum(idt)>200 | first)){

    for(i in length(channel):1){
      if(sum(idt)==0 | traditional){
        idt <- rep(TRUE, nrow(ff@exprs))
      }

      adj <- 2
      if(traditional) adj <- 4

      k_peak <- openCyto:::.find_peaks(ff@exprs[idt, channel[i]], adjust = adj)
      k_peak <- k_peak[k_peak$y > max(k_peak$y)*0.05, ]
      k <- nrow(k_peak)

      idt_k_peak <- idt
      if(traditional){
        idt_k_peak <- idt_k_peak & ff@exprs[, channel[i]]<stats::quantile(ff@exprs[, channel], 0.95)
        k_peak <- openCyto:::.find_peaks(ff@exprs[idt_k_peak, channel[i]], adjust = adj)
        k_peak <- k_peak[k_peak$y > max(k_peak$y)*0.05, ]
        k <- nrow(k_peak)
      }

      # Apply gating based on the number of peaks and direction

      if(k==1 & first){
        # ... Gate single peak case
        k <- nrow(openCyto:::.find_peaks(ff@exprs[idt, channel[i]], adjust = 1))
        gate <- gate_tail_custom(fr = ff[idt, ], channel = channel[i], num_peaks = k, ref_peak = ifelse(k>1, 2, 1))

        if(traditional){
          k_peak <- openCyto:::.find_peaks(ff@exprs[idt, channel[i]], adjust = 1)
          k <- nrow(k_peak)
          if(k==1){
            gate <- gate_tail_custom(fr = ff[idt, ], channel = channel[i], num_peaks = k, ref_peak = k)
          }else{
            gate <- openCyto::gate_mindensity(fr=ff[idt, ], channel[i], cutpoint_min=cutpoint_min, cutpoint_max = cutpoint_max, min.count = 100, adjust = 1)
          }
        }

        gate@min <- ifelse(gate@min<cutpoint_min, cutpoint_min, gate@min)
        gate@min <- ifelse(gate@min>cutpoint_max, cutpoint_max, gate@min)

      }else if(k>2){
        # ... Gate multiple peak case
        gate <- gate_tail_custom(fr = ff[idt, ], channel = channel[i], num_peaks = k, ref_peak = ifelse(channel_sign[i]=="[-", k-1, k), side = ifelse(channel_sign[i]=="[-", "right", "left"))
        if(traditional){
          gate <- gate_tail_custom(fr = ff[idt, ], channel = channel[i], num_peaks = k, ref_peak = k)
        }
        gate@min <- ifelse(gate@min<cutpoint_min, cutpoint_min, gate@min)
        gate@min <- ifelse(gate@min>cutpoint_max, cutpoint_max, gate@min)
      }else if(k==2){
        # ... Gate two peak case
        gate <- openCyto::gate_mindensity(fr=ff[idt_k_peak, ], channel[i], cutpoint_min=cutpoint_min, cutpoint_max = cutpoint_max, min.count = 100)
        gate@min <- ifelse(gate@min>cutpoint_max, cutpoint_max, gate@min)
      }else{
        gate <- gate_bg[[channel[i]]]
      }

      # Apply gating for specific channel properties
      if(beggining[i]){
        # ... Handle channels starting with '^'
        k_peak <- openCyto:::.find_peaks(ff@exprs[idt, channel[i]], adjust = 0.5)
        k_peak <- k_peak[order(k_peak$x, decreasing = F), ]
        rownames(k_peak) <- 1:nrow(k_peak)
        max_peak <- which.max(k_peak$y)
        gate <- gate_tail_custom(fr = ff[idt, ], channel = channel[i], side = "left", ref_peak = max_peak, num_peaks = nrow(k_peak), adjust = 0.5)
      }
      if(ending[i]){
        # ... Handle channels ending with '$'
        gate <- gate_tail_custom(fr = ff[idt, ], channel = channel[i], side = "right")
      }

      # Update idt based on gated events
      events_bg <-  do.call(channel_sign[i], args = list(ff@exprs[, names(gate@min)], gate@min))
      idt <- idt & events_bg
      idt_list[[i]] <- idt
      gate_bg[[channel[i]]] <- gate
    }

    if(traditional){
      idt <- Reduce("&", idt_list)
    }

    if(sum(idt)>200){
      num_peaks <- lapply(channel, function(x){
        k_peak <- openCyto:::.find_peaks(ff@exprs[idt, x], adjust = 2)
        k_peak <- k_peak[k_peak$y > max(k_peak$y)*0.05, ]
        return(nrow(k_peak))
      })
    }

    first <- F
    iter <- iter+1
  }

  if(sum(idt)>200 & join){

    border <- list()
    target <- list()

    for(i in 1:length(channel)){
      border[[channel[i]]] <- eval(call(paste0(channel_sign[i], "]"), ff@exprs[idt, channel[i]]))
      target[[channel[i]]] <- eval(call(paste0(channel_sign[i], "["), ff@exprs[idt, channel[i]]))
    }

    border <- unlist(border)
    target <- unlist(target)

    sam <- data.frame(ff@exprs[idt, c(channel)], tms = Sys.time() + 1:sum(idt), id = gl(1, sum(idt)))
    sp::coordinates(sam) <- stats::as.formula(paste0("~", paste0(colnames(sam)[1:2], collapse = "+")))
    tr <- trip(sam, c("tms", "id"))
    g <- tripGrid(tr)
    target <- coordinates(g)[which.max(g$z), ]

    gate_peak <- openCyto:::.find_peaks(ff@exprs[idt, channel[2]], adjust = 2)
    gate_peak <- gate_peak[gate_peak$y > max(gate_peak$y)*0.05, ]
    k_peak <- min(nrow(gate_peak), 2)
    gate_bg[["join"]] <- openCyto::gate_flowclust_2d(fr=ff[idt, ], quantile = 0.9, xChannel = channel[1], yChannel = channel[2], K=k_peak, target = target, min.count = 100,
                                                     min = border)
    idt <- idt & flowCore::filter(ff, gate_bg[["join"]])@subSet

  }

  result = list(idt=idt, channel=channel, gate_bg=gate_bg)

  if(return_plot){
    p1 <- ggcyto::autoplot(ff, channel, bins = 100)
    for(gate in gate_bg){
      p1 <- p1 + ggcyto::geom_gate(gate)
    }
    result[["plot"]] <- p1
  }

  return(result)
}

#' Estimate Proportions and Create Plots
#'
#' This function reads FCS files, applies gating based on provided patterns, and estimates proportions for each pattern.
#' It also generates plots and outputs a summary file.
#'
#' @param filename The path to the FCS file.
#' @param id A character identifier for the file (default: NULL, uses filename).
#' @param output.dir The output directory for plots and summary files.
#' @param info_panel A data frame containing information about patterns and labels.
#' @param cutpoint_min A numeric value setting a minimum threshold for gating.
#' @param cutpoint_max A numeric value setting a maximum threshold for gating.
#'
#' @keywords flowCore
#' @export
#'
#' @examples
#' doEstimateProportion(
#'   filename = "path/to/your/fcs_file.fcs",
#'   output.dir = "output_directory",
#'   info_panel = info_panel_data_frame
#' )
#'
doEstimateProportion <- function(filename, id = NULL, output.dir, info_panel, cutpoint_min = 1, cutpoint_max = 500000) {

  # Set default id if not provided
  if (is.null(id)) id <- filename

  # Read FCS file and calculate total events
  ff.raw <- flowCore::read.FCS(filename)
  n.total <- nrow(ff.raw)

  # Create a data frame to store results
  df_traditional <- data.frame(File = filename, ID = id)
  df_traditional$N.total <- n.total
  df_traditional[, info_panel$Label] <- NA

  # Create a list to store plots
  plot_traditional <- list()

  # Iterate through each pattern in the info_panel
  for (i in 1:nrow(info_panel)) {

    # Parse pattern and label from info_panel
    pattern_list <- info_panel$Pattern[i]
    label <- info_panel$Label[i]

    # Split pattern into channel components
    pattern_list <- unlist(strsplit(pattern_list, "[.]"))
    pattern_list <- sapply(pattern_list, strsplit, "[;]")

    # Create a list to store plots for each subpattern
    plot_list <- list()
    ff <- ff.raw

    # Iterate through each subpattern in the pattern list
    for (pattern in pattern_list) {

      # Split channel names and signs from pattern
      channel_sign <- sapply(pattern, function(x) substring(x, nchar(x), nchar(x)))
      channel <- sapply(pattern, function(x) substring(x, 1, nchar(x) - 1))

      # Handle channels starting with '^' or ending with '$'
      beggining <- sapply(channel, function(x) substring(x, nchar(x), nchar(x)) == "^")
      ending <- sapply(channel, function(x) substring(x, nchar(x), nchar(x)) == "$")
      channel[beggining] <- sapply(channel[beggining], function(x) substr(x, 1, nchar(x) - 1))
      channel[ending] <- sapply(channel[ending], function(x) substr(x, 1, nchar(x) - 1))

      # Add '^' or '$' to channel_sign accordingly
      channel_sign[beggining] <- paste0("^", channel_sign[beggining])
      channel_sign[ending] <- paste0("$", channel_sign[ending])

      # Map channel names to corresponding FCS parameter names
      names(channel) <- channel
      names(channel_sign) <- names(channel)
      if (all(!names(channel) %in% ff@parameters@data$desc)) {
        message("Error --> ", id)
        return(0)
      }

      # Construct channel names with signs
      channel <- ff@parameters@data$name[match(names(channel), ff@parameters@data$desc)]
      names(channel) <- names(channel_sign)
      channel <- paste0(channel, channel_sign)

      # Apply gating and get results
      res_gate <- apply_gate(ff, channel, cutpoint_min = cutpoint_min, cutpoint_max = cutpoint_max, return_plot = FALSE, join = FALSE, traditional = TRUE)

      idt_bg <- res_gate$idt
      gate_bg <- res_gate$gate_bg
      channel_bg <- res_gate$channel

      # Prepare plot and update ff for the next iteration
      if (length(gate_bg) == 0) gate_bg <- NULL
      plot_list[[paste0(pattern, collapse = ":")]] <- do_ggcyto(ff, channel_bg, gate_bg, logicle_chnls = channel_bg, main = label)
      ff <- ff[idt_bg, ]
    }

    # Update results data frame and plot list
    df_traditional[, label] <- nrow(ff)
    plot_traditional[[label]] <- plot_list[[length(plot_list)]]

    # Create output directory if it doesn't exist
    if (!dir.exists(paste0(output.dir, "/", label))) dir.create(paste0(output.dir, "/", label))

    # Write gated FCS file
    flowCore::write.FCS(ff, filename = paste0(output.dir, "/", label, "/", id, ".fcs"))
  }

  # Arrange and save plots
  nrow <- ifelse(length(plot_traditional) %% 2 == 0, length(plot_traditional) / 2, length(plot_traditional) / 2 + 1)
  result_plot <- ggpubr::ggarrange(plotlist = plot_traditional, ncol = 2, nrow = nrow)
  result_filename <- paste0(output.dir, "/Traditional_counts.txt")

  # Create plots subdirectory if it doesn't exist
  if (!dir.exists(paste0(output.dir, "/plots"))) dir.create(paste0(output.dir, "/plots"))

  # Save the arranged plot as a PDF file
  ggplot2::ggsave(paste0(output.dir, "/plots/", gsub(".fcs$", "", id), ".pdf"), plot = result_plot, device = "pdf", width = 15, height = 6 * nrow)

  # Add timestamp to the results data frame and write to a summary file
  df_traditional$Sys_Time <- Sys.time()
  write.table(
    x = df_traditional,
    file = result_filename,
    col.names = !file.exists(result_filename),
    sep = "\t",
    row.names = FALSE,
    append = file.exists(result_filename),
    quote = FALSE
  )
}


doEstimateProportion2 <- function (filename, id = NULL, output.dir, info_panel, cutpoint_min = 1,
                                   cutpoint_max = 5e+05) {
  if (is.null(id))
    id <- filename
  ff.raw <- flowCore::read.FCS(filename)
  n.total <- nrow(ff.raw)
  df_traditional <- data.frame(File = filename, ID = id)
  df_traditional$N.total <- n.total
  df_traditional[, info_panel$Label] <- NA
  plot_traditional <- list()
  for (i in 1:nrow(info_panel)) {
    pattern_list <- info_panel$Pattern[i]
    label <- info_panel$Label[i]
    pattern_list <- unlist(strsplit(pattern_list, "[.]"))
    pattern_list <- sapply(pattern_list, strsplit, "[;:]")
    prev_gate <- NULL
    if(any(grepl(":", names(pattern_list)))){
      prev_gate <- names(pattern_list)[which(grepl(":", names(pattern_list)))-1]
    }
    plot_list <- list()
    ff <- ff.raw
    for (pattern in pattern_list) {
      pattern_name <- names(which(sapply(pattern_list, function(x) all(pattern%in%x))))
      channel_sign <- sapply(pattern, function(x) substring(x,
                                                            nchar(x), nchar(x)))
      channel <- sapply(pattern, function(x) substring(x,
                                                       1, nchar(x) - 1))
      beggining <- sapply(channel, function(x) substring(x,
                                                         nchar(x), nchar(x)) == "^")
      ending <- sapply(channel, function(x) substring(x,
                                                      nchar(x), nchar(x)) == "$")
      channel[beggining] <- sapply(channel[beggining],
                                   function(x) substr(x, 1, nchar(x) - 1))
      channel[ending] <- sapply(channel[ending], function(x) substr(x,
                                                                    1, nchar(x) - 1))
      channel_sign[beggining] <- paste0("^", channel_sign[beggining])
      channel_sign[ending] <- paste0("$", channel_sign[ending])
      names(channel) <- channel
      names(channel_sign) <- names(channel)
      if (all(!names(channel) %in% ff@parameters@data$desc)) {
        message("Error --> ", id)
        return(0)
      }
      channel <- ff@parameters@data$name[match(names(channel),
                                               ff@parameters@data$desc)]
      names(channel) <- names(channel_sign)
      channel <- paste0(channel, channel_sign)
      if(grepl(":", pattern_name)){
        res_gate <- apply_gate_cuadrant(ff = ff, channel = channel, channel_sign = channel_sign, main=label, prev_gate=prev_gate, output.dir=output.dir)
        plot_gate <- TRUE
      }else{
        res_gate <- apply_gate(ff, channel, cutpoint_min = cutpoint_min,
                               cutpoint_max = cutpoint_max, return_plot = FALSE,
                               join = FALSE, traditional = TRUE)
        plot_gate <- FALSE
      }

      idt_bg <- res_gate$idt
      gate_bg <- res_gate$gate_bg
      channel_bg <- res_gate$channel
      if (length(gate_bg) == 0)
        gate_bg <- NULL
      if(!plot_gate){
        plot_gate <- do_ggcyto(ff,
                               channel_bg, gate_bg, logicle_chnls = channel_bg,
                               main = label)
      }else{
        plot_gate <- res_gate$plot
      }
      plot_list[[paste0(pattern, collapse = ":")]] <- plot_gate
      ff <- ff[idt_bg, ]
    }
    df_traditional[, label] <- nrow(ff)
    plot_traditional[[label]] <- plot_list[[length(plot_list)]]
    if (!dir.exists(paste0(output.dir, "/", label)))
      dir.create(paste0(output.dir, "/", label))
    flowCore::write.FCS(ff, filename = paste0(output.dir,
                                              "/", label, "/", id, ".fcs"))
  }
  nrow <- ifelse(length(plot_traditional)%%2 == 0, length(plot_traditional)/2,
                 length(plot_traditional)/2 + 1)
  result_plot <- ggpubr::ggarrange(plotlist = plot_traditional,
                                   ncol = 2, nrow = nrow)
  result_filename <- paste0(output.dir, "/Traditional_counts.txt")
  if (!dir.exists(paste0(output.dir, "/plots")))
    dir.create(paste0(output.dir, "/plots"))
  ggplot2::ggsave(paste0(output.dir, "/plots/", gsub(".fcs$",
                                                     "", id), ".pdf"), plot = result_plot, device = "pdf",
                  width = 15, height = 6 * nrow)
  df_traditional$Sys_Time <- Sys.time()
  write.table(x = df_traditional, file = result_filename,
              col.names = !file.exists(result_filename), sep = "\t",
              row.names = FALSE, append = file.exists(result_filename),
              quote = FALSE)
}

gate_flowclust_2d_custom <- function (fr, xChannel, yChannel, filterId = "", K = 2, usePrior = "no",
                                      prior = list(NA), trans = 0, min.count = -1, max.count = -1,
                                      nstart = 1, plot = FALSE, target = NULL, transitional = FALSE,
                                      quantile = 0.9, translation = 0.25, transitional_angle = NULL,
                                      min = NULL, max = NULL, prev_gate="", output.dir, ...)
{
  options(cores = 1L)
  if (!is.null(target)) {
    target <- as.numeric(target)
    if (length(target) != 2) {
      warning("The 'target' location must be a numeric vector of length 2.\n               Using largest cluster instead...")
      target <- NULL
    }
  }
  if (!(is.null(min) && is.null(max))) {
    fr <- .truncate_flowframe(fr, channels = c(xChannel,
                                               yChannel), min = min, max = max)
  }
  if (usePrior == "yes" && identical(prior, list(NA))) {
    prior <- prior_flowclust(fr = fr, channels = c(xChannel,
                                                   yChannel), K = K)
  }

  if(!dir.exists(paste0(output.dir, "/tmp"))) dir.create(paste0(output.dir, "/tmp"))
  rds_obj <- paste0(output.dir, "/tmp/", basename(fr@description$FILENAME), "-", gsub(";", "", prev_gate), "-", xChannel, yChannel, "-", K, ".rds")
  if(file.exists(rds_obj)){
    message("Reading: ", rds_obj)
    tmix_results <- readRDS(rds_obj)
  }else{
    tmix_results <- flowClust(fr, varNames = c(xChannel, yChannel),
                              K = K, trans = trans, usePrior = usePrior, prior = prior,
                              min.count = min.count, max.count = max.count, nstart = nstart,
                              ...)
    message("Saving: ", rds_obj)
    saveRDS(tmix_results, rds_obj)
  }
  if (class(tmix_results) == "try-error") {
    tmix_results <- new("flowClust", varNames = c(xChannel,
                                                  yChannel), K = K, w = prior$w0, mu = prior$Mu0,
                        sigma = prior$Lambda0, nu = 4, prior = prior, ruleOutliers = c(0,
                                                                                       quantile, quantile))
  }
  fitted_means <- getEstimates(tmix_results)$locations
  if (is.null(target)) {
    cluster_selected <- which.max(tmix_results@w)
  }
  else {
    target_dist <- as.matrix(dist(rbind(fitted_means, target)))
    target_dist <- tail(target_dist, n = 1)[seq_len(K)]
    cluster_selected <- which.min(target_dist)
  }
  if (!transitional) {
    flowClust_gate <- openCyto:::.getEllipseGate(filter = tmix_results,
                                                 include = cluster_selected, quantile = quantile,
                                                 trans = trans)
  }
  else {
    chisq_quantile <- qchisq(quantile, df = 2)
    tol <- sqrt(.Machine$double.eps)
    xbar <- tmix_results@mu[cluster_selected, ]
    Sigma <- tmix_results@sigma[cluster_selected, , ]
    Sigma_eigen <- eigen(Sigma, symmetric = TRUE)
    u1 <- Sigma_eigen$vectors[, 1]
    u2 <- Sigma_eigen$vectors[, 2]
    lambda1 <- Sigma_eigen$values[1]
    lambda2 <- Sigma_eigen$values[2]
    u1_angle <- atan2(u1[2], u1[1])
    if (u1_angle < 0) {
      u1_angle <- u1_angle + 2 * pi
    }
    u2_angle <- atan2(u2[2], u2[1])
    if (u2_angle < 0) {
      u2_angle <- u2_angle + 2 * pi
    }
    if (u1_angle > pi) {
      R <- .rotation_matrix(pi)
      u1 <- as.vector(R %*% u1)
    }
    if (u2_angle > pi) {
      R <- .rotation_matrix(pi)
      u2 <- as.vector(R %*% u2)
    }
    u1_angle <- atan2(u1[2], u1[1])
    if (u1_angle < 0) {
      u1_angle <- u1_angle + 2 * pi
    }
    u2_angle <- atan2(u2[2], u2[1])
    if (u2_angle < 0) {
      u2_angle <- u2_angle + 2 * pi
    }
    eigen_angles <- c(u1_angle, u2_angle)
    if (is.null(transitional_angle)) {
      which_pos_quadrant <- which(0 < eigen_angles & eigen_angles <
                                    pi/2)
      transitional_angle <- eigen_angles[which_pos_quadrant]
      if (which_pos_quadrant == 1) {
        axis <- sqrt(lambda1 * chisq_quantile) * u1
        axis_perp <- sqrt(lambda2 * chisq_quantile) *
          u2
      }
      else {
        axis <- sqrt(lambda2 * chisq_quantile) * u2
        axis_perp <- sqrt(lambda1 * chisq_quantile) *
          u1
      }
    }
    else {
      theta_u1 <- transitional_angle - eigen_angles[1]
      theta_u2 <- transitional_angle - (pi/2) - eigen_angles[2]
      R1 <- .rotation_matrix(theta_u1)
      R2 <- .rotation_matrix(theta_u2)
      u1 <- as.vector(R1 %*% u1)
      u2 <- as.vector(R2 %*% u2)
      axis <- sqrt(lambda1 * chisq_quantile) * u1
      axis_perp <- -sqrt(lambda2 * chisq_quantile) * u2
    }
    gate_location <- xbar + translation * axis
    x <- exprs(fr)[, xChannel]
    y <- exprs(fr)[, yChannel]
    x_min <- min(x) - sd(x)
    y_min <- min(y) - sd(y)
    x_max <- max(x) + sd(x)
    y_max <- max(y) + sd(y)
    first_vertex <- gate_location + axis_perp
    fifth_vertex <- gate_location - axis_perp
    if (0 <= transitional_angle && transitional_angle <=
        pi/2) {
      second_vertex <- c(first_vertex[1], y_max)
      third_vertex <- c(x_max, y_max)
      fourth_vertex <- c(x_max, fifth_vertex[2])
    }
    else if (pi/2 < transitional_angle && transitional_angle <=
             pi) {
      second_vertex <- c(x_min, first_vertex[2])
      third_vertex <- c(x_min, y_max)
      fourth_vertex <- c(fifth_vertex[1], y_max)
    }
    else if (pi < transitional_angle && transitional_angle <=
             3 * pi/2) {
      second_vertex <- c(first_vertex[1], y_min)
      third_vertex <- c(x_min, y_min)
      fourth_vertex <- c(x_min, fifth_vertex[2])
    }
    else {
      second_vertex <- c(x_max, first_vertex[2])
      third_vertex <- c(x_max, y_min)
      fourth_vertex <- c(fifth_vertex[1], y_min)
    }
    polygon_gate <- rbind(first_vertex, second_vertex, third_vertex,
                          fourth_vertex, fifth_vertex, first_vertex)
    colnames(polygon_gate) <- c(xChannel, yChannel)
    flowClust_gate <- polygonGate(filterId = filterId, .gate = polygon_gate)
  }
  posteriors <- list(mu = tmix_results@mu, lambda = tmix_results@lambda,
                     sigma = tmix_results@sigma, nu = tmix_results@nu)
  if (plot) {
    plot(data = fr, tmix_results, main = filterId)
    if (transitional) {
      lines(rbind(xbar - axis, xbar + axis), col = "darkgreen")
      lines(rbind(xbar - axis_perp, xbar + axis_perp),
            col = "darkgreen")
      lines(polygon_gate, col = "red")
      lines(rbind(gate_location - axis_perp, gate_location +
                    axis_perp), col = "red")
      points(polygon_gate, col = "red", pch = 16)
    }
  }
  if (class(flowClust_gate) == "polygonGate")
    openCyto:::fcPolygonGate(flowClust_gate, prior, posteriors)
  else openCyto:::fcEllipsoidGate(flowClust_gate, prior, posteriors)
}


apply_gate_cuadrant <- function(ff, channel, channel_sign, main="", prev_gate="", output.dir){
  library(flowClust)
  channel <- sapply(channel, function(x) substr(x, 1, nchar(x) - 1))
  logicle_chnls <- channel
  lgcl <- flowCore::logicleTransform()
  trans <- flowCore::transformList(logicle_chnls, lgcl)
  fC_plot <- ggcyto::transform(ff, trans)

  range_x <- quantile(fC_plot@exprs[, channel[1]], c(0.1, 0.9))#range(fC_plot@exprs[, channel[1]])
  range_y <- quantile(fC_plot@exprs[, channel[2]], c(0.1, 0.9))#range(fC_plot@exprs[, channel[2]])

  cl1 <- c(range_x[1], range_y[1])
  cl2 <- c(range_x[1], range_y[2])
  cl3 <- c(range_x[2], range_y[1])
  cl4 <- c(range_x[2], range_y[2])
  K_cluster <- 12

  getGate <- function(fC_plot, channel, K_cluster, target, quantile=0.8, output.dir){

    gate_bg_v1 <- gate_flowclust_2d_custom(fC_plot, yChannel  = channel[2], xChannel = channel[1], K=K_cluster, target = target, quantile = quantile, prev_gate=prev_gate, output.dir=output.dir)
    gate_bg_v2 <- gate_flowclust_2d_custom(fC_plot, yChannel  = channel[2], xChannel = channel[1], K=4, target = target, quantile = quantile, prev_gate=prev_gate, output.dir=output.dir)

    if(sqrt(sum(gate_bg_v1@mean-cl2)^2)>0.5){
      K_cluster <- K_cluster*2
      gate_bg_v1 <- gate_flowclust_2d_custom(fC_plot, yChannel  = channel[2], xChannel = channel[1], K=K_cluster, target = target, quantile = quantile, prev_gate=prev_gate, output.dir=output.dir)
    }

    if(sum(flowCore::filter(fC_plot, gate_bg_v1)@subSet & flowCore::filter(fC_plot, gate_bg_v2)@subSet)/sum(flowCore::filter(fC_plot, gate_bg_v1)@subSet)>0.25){
      gate_bg <- gate_bg_v2
    }else{
      gate_bg <- gate_bg_v1
    }

    return(gate_bg)
  }

  idt <- rep(TRUE, nrow(fC_plot))
  if(channel_sign[1]=="-" & channel_sign[2]=="-"){
    gate_bg <- getGate(fC_plot, channel, K_cluster, cl1, output.dir=output.dir)
    idt <- idt & flowCore::filter(fC_plot, gate_bg)@subSet
  }else if(channel_sign[1]=="-" & channel_sign[2]=="+"){
    gate_bg <- getGate(fC_plot, channel, K_cluster, cl2, output.dir=output.dir)
    idt <- idt & flowCore::filter(fC_plot, gate_bg)@subSet
  }else if(channel_sign[1]=="+" & channel_sign[2]=="-"){
    gate_bg <- getGate(fC_plot, channel, K_cluster, cl3, output.dir=output.dir)
    idt <- idt & flowCore::filter(fC_plot, gate_bg)@subSet
  }else{
    gate_bg <- getGate(fC_plot, channel, K_cluster, cl4, output.dir=output.dir)
    idt <- idt & flowCore::filter(fC_plot, gate_bg)@subSet
  }

  p <- ggcyto::autoplot(fC_plot, x = channel[1], y = channel[2], bins = 100)+ggcyto::geom_gate(gate_bg)
  p <- p + ggplot2::theme_bw() + ggplot2::theme(strip.text.x = ggplot2::element_text(size = 0,
                                                                                     colour = "orange", angle = 90))
  p <- ggcyto::as.ggplot(p)
  labels <- fC_plot@parameters@data$desc[match(channel, fC_plot@parameters@data$name)]
  labels[is.na(labels)] <- channel[is.na(labels)]
  for (i in 1:length(labels)) {
    if (labels[i] != channel[i]) {
      labels[i] <- paste0(labels[i], " (", channel[i],
                          ")")
    }
  }
  labels[match(logicle_chnls, channel)] <- paste0(labels[match(logicle_chnls,
                                                               channel)], " - logicleTransform")
  p <- p + ggplot2::xlab(labels[1]) + ggplot2::ylab(labels[2])
  p <- p + ggplot2::ggtitle(label = main) + ggplot2::theme(plot.title = ggplot2::element_text(size = 15,
                                                                                              face = "bold", hjust = 0.5))
  p <- p + ggplot2::theme(legend.position = "none")
  result <- list(idt = idt, channel = channel, gate_bg = list("join"=gate_bg), plot=p)

  return(result)
}

#' Run Proportion Estimation and Gating for Multiple Files
#'
#' This function runs the proportion estimation and gating process for multiple files.
#' It reads a log file to identify completed files, applies the estimation and gating functions,
#' and outputs results to the specified directory.
#'
#' @param log_file_track The path to the log file containing file tracking information.
#' @param panel_estimate File path to the panel for proportion estimation.
#' It contains information about patterns and labels.
#' @param output The output directory for results.
#' @param cluster An optional cluster object for parallel processing (default: NULL).
#'
#' @keywords flowCore
#' @export
#'
#' @examples
#' runEstimateProportion(
#'   log_file_track = "output_directory/log_file.txt",
#'   panel_estimate = "panel_estimate.csv",
#'   output = "output_directory"
#' )
#'
runEstimateProportion <- function(log_file_track, panel_estimate, output, cluster = NULL) {

  # Read info_panel
  info_panel <- read.delim(panel_estimate)

  # Read the log file to identify completed files
  counts <- read.delim(log_file_track, header = FALSE)
  completed_files <- counts$V3[counts$V6 %in% "Completed"]

  # Create an output directory for results
  output.dir <- file.path(output, "Results")
  if (!dir.exists(output.dir)) dir.create(output.dir)

  # Set up a log file for errors during processing
  log_file_traditional <- file.path(output, "Log_file_error_estimate_proportion.txt")
  log_file_error_traditional <- function(messages) {
    cat(paste(messages, collapse = "\n"), "\n", file = log_file_traditional, append = TRUE)
  }

  # Helper function to process a single file
  process_file <- function(filename, output.dir, info_panel) {
    tryCatch({
      id <- gsub(".fcs$", "", basename(filename))
      doEstimateProportion2(filename = filename, id = id, info_panel = info_panel, output.dir = output.dir, cutpoint_min = 0, cutpoint_max = 30000)
    },
    error=function(e) {
      log_file_error_traditional(paste("File: ", filename, "\n", e))
    })
  }

  # Process files in parallel if a cluster is provided
  if (!is.null(cluster)) {
    requireNamespace('foreach')
    # Parallel processing using foreach
    foreach(i = completed_files, .packages = c("flowCore")) %dopar% {
      process_file(i, output.dir, info_panel)
    }
  } else {
    # Process files sequentially
    for (i in completed_files) {
      process_file(i, output.dir, info_panel)
    }
  }
}




