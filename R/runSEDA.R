

#' Extracts and processes FCS files with optional metadata and marker name.
#'
#' This function extracts data from FCS files, processes it, and returns a processed
#' single-cell experiment (SCE) object. It supports filtering and transformation of
#' the data, with the option to provide metadata and a marker name for the experiment.
#'
#' @param fcs_path A character vector of file paths to FCS files.
#' @param metadata Optional data frame containing metadata for FCS files.
#' @param maker_name Optional name for the marker to be added in metadata.
#' @return A processed single-cell experiment (SCE) object.
#' @export
#'
#' @examples
#' doSEDA(fcs_path = c("file1.fcs", "file2.fcs"), metadata = meta_data, maker_name = "Marker1")
#'
#'
doSEDA <- function(fcs_path, metadata = NULL, maker_name = NULL) {
  # Filter FCS files based on the number of rows
  keep <- sapply(fcs_path, function(x) {
    tryCatch({
      nrow(CytoTree::runExprsExtract(x)) > 2
    }, error = function(e) {
      return(FALSE)
    })
  })

  fcs_path <- fcs_path[keep]

  # Handle metadata if not provided
  if (is.null(metadata)) metadata <- data.frame(ID = basename(fcs_path), file_path = fcs_path)

  fcs_path_file <- basename(fcs_path)
  metadata <- metadata[match(fcs_path_file, metadata$ID), ]
  rownames(metadata) <- fcs_path

  # Extract marker names and inclusions
  fcs <- CytoTree::runExprsExtract(fcs_path[1])
  exclusions <- colnames(fcs)[grep("NA", colnames(fcs))]
  inclusions <- colnames(fcs)[!colnames(fcs) %in% exclusions]
  inclusions <- unlist(lapply(strsplit(inclusions, "<"), `[`, 1))
  markernames <- colnames(fcs)[!colnames(fcs) %in% exclusions]
  markernames <- gsub(">", "", markernames)
  markernames <- unlist(lapply(strsplit(markernames, "<"), `[`, 2))
  names(markernames) <- inclusions

  # Process FCS data
  sce <- scDataviz::processFCS(
    files = fcs_path,
    metadata = metadata,
    transformation = TRUE,
    transFun = function(x) asinh(x),
    downsampleVar = NULL,
    downsample = NULL,
    colsRetain = inclusions,
    asinhFactor = 150,
    filter = TRUE,
    bgNoiseThreshold = -Inf,
    colsDiscard = exclusions,
    newColnames = markernames
  )

  # Add marker name to metadata
  if (!is.null(maker_name)) sce@metadata$marker <- maker_name

  return(sce)
}

#' Perform Min-Max Gene Sampling on Expression Matrix
#'
#' This function performs Min-Max gene sampling on an expression matrix to select a subset of cells
#' that represent a broad expression range across genes. The algorithm aims to maximize the diversity
#' of gene expression profiles among the selected cells.
#'
#' @param expression_matrix A numeric matrix where rows represent genes and columns represent cells.
#' @param target_cells The desired number of cells to be selected from the expression matrix.
#' @return A vector of indices representing the selected cells.
#'
#' @examples
#' expression_matrix <- matrix(runif(1000), nrow = 100, ncol = 10)
#' selected_cells <- min_max_sampling(expression_matrix, target_cells = 5)
#'
#' @export
#'
min_max_sampling <- function(expression_matrix, target_cells) {
  # Calculate the expression range for each gene
  expression_range <- apply(expression_matrix, 1, function(x) max(x) - min(x))

  # Sort genes by expression range in descending order
  sorted_genes <- order(expression_range, decreasing = TRUE)

  # Initialize sets of selected and remaining cells
  selected_cells <- integer(0)
  remaining_cells <- 1:nrow(expression_matrix)

  while (length(selected_cells) < target_cells && length(remaining_cells) > 0) {
    selected_gene <- sorted_genes[1]
    remaining_genes <- sorted_genes[-1]

    gene_expression <- expression_matrix[selected_gene, ]
    gene_range <- max(gene_expression) - min(gene_expression)

    # Find the cell furthest from the current range
    distance_to_range <- abs(gene_expression - min(gene_expression)) / gene_range +
      abs(gene_expression - max(gene_expression)) / gene_range
    selected_cell <- which.max(distance_to_range)

    # Add the selected cell to the set and update remaining cells
    selected_cells <- c(selected_cells, selected_cell)
    remaining_cells <- setdiff(remaining_cells, selected_cell)
  }

  # Return the reduced expression matrix
  return(selected_cells)
}


#' Downsampling with Proportional Outliers
#'
#' This function performs downsampling on a single-cell gene expression matrix while maintaining the proportion of outliers.
#'
#' @param expression_matrix A matrix of gene expression data where rows represent cells and columns represent genes.
#' @param target_cells The desired number of cells after downsampling.
#'
#' @return A downsampled gene expression matrix.
#'
#' @examples
#' # Suppose 'expression_matrix' is your gene expression matrix and 'target_cells' is the desired number of cells after downsampling
#' sampled_expression_matrix <- downsampling_with_outliers(expression_matrix, target_cells)
#'
#'
#' @export
downsampling_with_outliers <- function(expression_matrix, target_cells) {
  # Calculate the dispersion measure, e.g., standard deviation
  sd_per_cell <- apply(expression_matrix, 1, sd)

  # Define a threshold to consider cells as outliers
  outlier_threshold <- quantile(sd_per_cell, 0.95)

  # Identify outlier cells
  outlier_cells <- which(sd_per_cell > outlier_threshold)

  # Calculate the proportion of outliers in the original matrix
  proportion_outliers <- length(outlier_cells) / nrow(expression_matrix)

  # Calculate how many non-outlier cells are needed in downsampling
  non_outlier_target <- target_cells * (1 - proportion_outliers)

  # Order non-outlier cells based on some measure (e.g., mean)
  non_outlier_means <- rowMeans(expression_matrix[setdiff(1:nrow(expression_matrix), outlier_cells), ])
  non_outlier_means <- non_outlier_means[!names(non_outlier_means)%in%names(outlier_cells)]
  sorted_non_outlier_indices <- order(non_outlier_means)

  # Select sampled non-outlier cell indices
  sampled_non_outlier_indices <- names(sorted_non_outlier_indices)[1:non_outlier_target]
  sampled_indices <- c(sampled_non_outlier_indices, names(outlier_cells))


  return(which(sampled_indices%in%rownames(expression_matrix)))
}




#' Run SEDA Analysis
#'
#' This function performs the SEDA analysis on FCS data, including optional marker-specific analysis and downsampling.
#'
#' @param output Output directory for saving the analysis results.
#' @param metadata Optional data frame containing additional metadata information.
#' @param marker Whether to perform marker-specific analysis.
#' @param downsampling The type of downsampling to perform ("random", "minMaxSamplingSEACELLS", "minMaxSampling", or none).
#' @param k_downsampling Percentage of umber of cells to downsample for each ID.
#' @param seed Seed for reproducibility if downsampling is used.
#' @param batch Whether to perform batch effect correction.
#' @param conda_env Conda environment to use for running SEACells.
#'
#' @import scDataviz
#'
#' @export
runSEDA <- function(output, metadata=NULL, marker=FALSE, downsampling="random", k_downsampling=0.01, seed=NULL, batch=FALSE, conda_env=NULL){

  if(!is.null(seed)) set.seed(seed)

  # Get a list of cleaned FCS file paths
  fcs_clean <- list.files(file.path(output, "fcs_clean"), full.names = T)

  # Check if exists
  if(file.exists(file.path(output, "SEDA", "sce_clean.rds"))){
    sce_clean <- readRDS(file.path(output, "SEDA", "sce_clean.rds"))
  }else{
    # Perform SEDA analysis on cleaned FCS files
    sce_clean <- doSEDA(fcs_path = fcs_clean)
    saveRDS(sce_clean, file.path(output, "SEDA", "sce_clean.rds"))
  }

  if(file.exists(file.path(output, "SEDA", "sce_final.rds"))){
    sce_final <- readRDS(file.path(output, "SEDA", "sce_final.rds"))
  }else if(marker){
    sce_clean@metadata$Marker <- "ND"

    # Extract marker-specific analysis results
    results_path <- list.dirs(file.path(output, "Results"))[-1]
    results_path <- results_path[-grep("plots", results_path)]
    fcs_path <- lapply(results_path, list.files, full.names=T)
    maker_name <- lapply(results_path, basename)

    # Perform marker-specific SEDA analysis
    sce_list <- lapply(1:length(maker_name), function(i) doSEDA(fcs_path = fcs_path[[i]], maker_name = maker_name[[i]]))

    id.sce_clean <- apply(sce_clean@assays@data$scaled, 2, paste0, collapse="-")

    for(i in 1:length(sce_list)){
      sce <- sce_list[[i]]
      id.sce <- apply(sce@assays@data$scaled, 2, paste0, collapse="-")
      sce_clean@metadata$Marker[id.sce_clean%in%id.sce] <- maker_name[[i]]
    }

    # Clean up temporary variables
    rm(sce_list)
    sce_final <- sce_clean
    rm(sce_clean)

  }else{
    sce_final <- sce_clean
    rm(sce_clean)
  }

  saveRDS(sce_final, file.path(output, "SEDA", "sce_final.rds"))

  # Remove non-informative genes and cells
  idt <- !apply(is.na(sce_final@assays@data$scaled), 2, any)
  sce_final <- sce_final[, idt]
  sce_final@metadata <- sce_final@metadata[idt, ]

  if(is.null(conda_env)){
    message("NULL in conda env. Then change to minMaxSampling.")
    downsampling <- "minMaxSampling"
  }

  # Perform specified downsampling method
  if(downsampling%in%"random"){
    idt <- unlist(lapply(unique(sce_final@metadata$ID), function(x) sample(which(sce_final@metadata$ID%in%x), length(which(sce_final@metadata$ID%in%x))*k_downsampling)))
  }else if(downsampling%in%"minMaxSamplingSEACELLS"){
    if(!dir.exists(file.path(output, "tmp"))) dir.create(file.path(output, "tmp"))
    tmp_file <- file.path(output, "tmp", "scaled-counts.txt")
    write.table(t(sce_final@assays@data$scaled), tmp_file, row.names = F, col.names = F, sep = "\t", quote = F)
    idt_file <- run_min_max_sampling(file_counts = tmp_file, output = file.path(output, "tmp"), n_SEACells = floor(ncol(sce_final)*k_downsampling), conda_env=conda_env)
    idt <- as.numeric(read.delim(idt_file, header = F)[, 1])
  }else if(downsampling%in%"minMaxSampling"){
    idt <- unlist(lapply(unique(sce_final@metadata$ID), function(x) downsampling_with_outliers(t(sce_final@assays@data$scaled)[which(sce_final@metadata$ID%in%x), ], floor(length(which(sce_final@metadata$ID%in%x))*k_downsampling))))
  }else{
    idt <- 1:nrow(ncol(sce_final))
  }

  # Apply downsampling to the SCE object
  sce_final <- sce_final[, idt]
  sce_final@metadata <- sce_final@metadata[idt, ]

  # Remove columns with zero variance
  idt <- matrixStats::colSds(sce_final@assays@data$scaled)==0
  sce_final <- sce_final[, !idt]
  sce_final@metadata <- sce_final@metadata[!idt, ]

  # Perform batch effect correction if specified
  if(batch & "batch" %in% colnames(sce_final@metadata)){
    # Apply batch effect correction using Harmony
    sce_final@assays@data$normexprs <- t(harmony::HarmonyMatrix(data_mat = t(sce_final@assays@data$scaled),
                                                             meta_data = sce_final@metadata,
                                                             vars_use = c("batch"),
                                                             do_pca = F, plot_convergence = T,
                                                             verbose = T))


    # Perform UMAP and clustering
    sce_final <- scDataviz::performUMAP(sce_final, assay = "normexprs")
    sce_final <- scDataviz::clusKNN(sce_final,
                        k.param = 20,
                        prune.SNN = 1/15,
                        resolution = 0.01,
                        algorithm = 2,
                        verbose = FALSE)

  }else{
    # Perform UMAP and clustering without batch effect correction
    sce_final <- scDataviz::performUMAP(sce_final)

    sce_final <- scDataviz::clusKNN(sce_final,
                                    k.param = 20,
                                    prune.SNN = 1/15,
                                    resolution = 0.01,
                                    algorithm = 2,
                                    verbose = FALSE)
  }

  if(!is.null(metadata) & "ID"%in%colnames(metadata)){
    keep <- colnames(metadata)[!colnames(metadata)%in%colnames(sce_final@metadata)]
    row_order <- match(sce_final@metadata$ID, metadata$ID)

    if(length(keep)>0 & length(row_order)>0){
      sce_final@metadata[, keep] <- NA
      sce_final@metadata[, keep] <- metadata[row_order, keep]
    }
  }

  if(!dir.exists(file.path(output, "SEDA"))) dir.create(file.path(output, "SEDA"))
  saveRDS(sce_final, file.path(output, "SEDA", "sce_CytoTree.rds"))
}


