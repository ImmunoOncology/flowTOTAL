#' Run MinMax-Sampling using SEACells
#'
#' @param file_counts File with raw counts
#' @param n_components Number of components
#' @param n_SEACells Number of cells to keep
#' @param ... Other arguments passed to the Python function
#'
#' @return Path to the file with the selected cell
#' @import reticulate
#'
run_min_max_sampling <- function(file_counts, output, n_SEACells, n_waypoint_eigs=5, n_comps_pca=10, build_kernel_on="X_umap") {
  code_args <- paste(
    system.file("python", "your_script.py", package = "flowTOTAL"),
    output, n_SEACells, n_waypoint_eigs, n_comps_pca, build_kernel_on
  )
  reticulate::py_run_file(code_args)
  return(file.path(output, "min-max-sampling-index.txt"))
}
