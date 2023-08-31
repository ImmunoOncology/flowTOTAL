#' Run MinMax-Sampling using SEACells
#'
#' This function performs min-max sampling using SEACells or falls back to a basic sampling method if there's an error with the conda environment.
#'
#' @param file_counts Path to the file containing counts data.
#' @param output Directory where the output files will be saved.
#' @param n_SEACells Number of SEACells.
#' @param n_waypoint_eigs Number of waypoint eigenvectors (default is 5).
#' @param n_comps_pca Number of PCA components (default is 10).
#' @param build_kernel_on Build kernel on which data (default is "X_pca").
#' @param conda_env Conda environment to use for running SEACells.
#'
#' @return Path to the file containing the indices used for min-max sampling.
#'
#' @examples
#' run_min_max_sampling(
#'   file_counts = "counts.txt",
#'   output = "output_folder",
#'   n_SEACells = 100,
#'   n_waypoint_eigs = 5,
#'   n_comps_pca = 10,
#'   build_kernel_on = "X_pca",
#'   conda_env = "my_env"
#' )
#'
#' @import reticulate
#'
run_min_max_sampling <- function(file_counts, output, n_SEACells, n_waypoint_eigs=5, n_comps_pca=10, build_kernel_on="X_pca", conda_env=NULL) {
  tryCatch({
    code_args <- paste(
      system.file("python", "seacells_wrapper.py", package = "flowTOTAL"),
      file_counts, output, n_SEACells, n_waypoint_eigs, n_comps_pca, build_kernel_on
    )
    message("Using ", conda_env)
    reticulate::use_condaenv(condaenv = conda_env, required = TRUE)
    message("Running ", code_args)
    reticulate::py_run_file(code_args)
    message("Finished.")
  },
  error=function(e) {
    message("Error when using conda env. Instead using random.")
    message(e)
    n_row <- nrow(read.delim(file_counts, header = F))
    idt <- sample(1:n_row, n_SEACells)
    write.table(idt, file.path(output, "min-max-sampling-index.txt"), col.names = F, row.names = F, sep = "\t", quote = F)
  })

  return(file.path(output, "min-max-sampling-index.txt"))
}
