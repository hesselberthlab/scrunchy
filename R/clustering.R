#' K-means clustering
#'
#' @param fsce [`FunctionalSingleCellExperiment`]
#' @param expt Data to use for calculating variable features
#'   (default is `rnaseq`). Must be present in `names(fsce)`.
#' @param k number of classes
#' @param method dimensionality reduction method for clustering (defaults to
#'   PCA)
#' @param n_dims specify the number of dimensions from "dr" to use for
#'   clustering, defaults to all dimensions
#' @param ... additional arguments to pass to [`stats::kmeans()`]
#'
#' @return fsce with `k_cluster` in `expt` colData.
#'
#' @examples
#' # calculate PCA for k-means default method
#' fsce <- calc_pca(fsce_small)
#'
#' fsce <- calc_kmeans(fsce, k = 6)
#' \dontrun{
#' library(SingleCellExperiment)
#' colData(fsce[["rnaseq"]], "k_cluster")
#' }
#'
#' @family clustering
#'
#' @export
calc_kmeans <- function(fsce,
                        expt = "rnaseq",
                        k,
                        method = "PCA",
                        n_dims = NULL,
                        ...) {
  ## check inputs
  if (!expt %in% names(fsce)) {
    stop(glue("expt `{expt}` not found in fsce "), call. = FALSE)
  }

  if (!method %in% names(reducedDims(fsce[[expt]]))) {
    stop(glue("method `{method}` not found in fsce"), call. = FALSE)
  }

  dr_mat <- reducedDim(fsce[[expt]], method)

  if (!is.null(n_dims)) {
    if (n_dims > ncol(dr_mat)) {
      stop("n_dims larger than dimensality reduction matrix", call. = FALSE)
    }
    dr_mat <- dr_mat[, 1:n_dims]
  }

  if (k > nrow(dr_mat)) {
    stop("k is larger than observations in input matrix", call. = FALSE)
  }

  km_res <- stats::kmeans(dr_mat, centers = k, ...)

  colData(fsce[[expt]])$k_cluster <- as.character(km_res$cluster)

  fsce
}
