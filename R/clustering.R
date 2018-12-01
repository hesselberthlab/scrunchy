#' Perform Kmeans clustering
#'
#' @param fce fce object
#' @param k specify k number of classes
#' @param expt experiment to cluster (either "sce" or "fsce, defaults to "sce")
#' @param dr which dimensionality reduction to use for clustering (defaults to PCA)
#' @param n_dims specify the number of dimensions from "dr" to use for clustering, defaults to all dimensions
#' @param ... additional arguments to pass to [`stats::kmeans()`]
#'
#' @export
run_kmeans <- function(fce,
                       k,
                       expt = "sce",
                       dr = "PCA",
                       n_dims = NULL,
                       ...) {
  ## check inputs
  if(!expt %in% names(assays(fce))) {
    stop("expt not found in fce object")
  }

  if(!dr %in% names(reducedDims(fce[[expt]]))) {
    stop("dr method not found in fce object")
  }

  dr_mat <- reducedDim(fce[[expt]], dr)

  if (!is.null(n_dims)){
    if(n_dims > ncol(dr_mat)){
      stop("n_dims larger than dimensality reduction matrix")
    }
    dr_mat <- dr_mat[, 1:n_dims]
  }

  if (k > nrow(dr_mat)) {
    stop("k is larger than observations in input matrix")
  }

  km_res <- stats::kmeans(dr_mat, centers = k, ...)

  colData(fce[[expt]])$k_cluster <- as.character(km_res$cluster)

  fce
}
