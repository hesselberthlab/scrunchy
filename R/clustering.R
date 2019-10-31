#' Run k-means clustering algorithm
#'
#' @param fsce An object of class [`FunctionalSingleCellExperiment`]
#' @param expt Data to use for calculating variable features
#'   (default is `rnaseq`). Must be present in `names(fsce)`.
#' @param k number of classes
#' @param method dimensionality reduction method for clustering (defaults to
#'   PCA)
#' @param n_dims specify the number of dimensions from "dr" to use for
#'   clustering, defaults to all dimensions
#' @param seed seed for reproducible result
#' @param ... additional arguments to pass to [`stats::kmeans()`]
#'
#' @return fsce with `k_cluster` in `expt` colData.
#'
#' @examples
#' # calculate PCA for k-means default method
#' fsce <- calc_pca(fsce_small)
#'
#' fsce <- cluster_kmeans(fsce, k = 6)
#'
#' SingleCellExperiment::colData(fsce[["rnaseq"]], "k_cluster")
#'
#' @family clustering functions
#'
#' @export
#'
cluster_kmeans <- function(fsce,
                        expt = "rnaseq",
                        k,
                        method = "PCA",
                        n_dims = NULL,
                        seed = NULL,
                        ...) {
  ## check inputs
  if (!expt %in% names(fsce)) {
    stop(glue("expt `{expt}` not found in fsce "), call. = FALSE)
  }

  if (!method %in% names(reducedDims(fsce[[expt]]))) {
    stop(glue("method `{method}` not found in expt"), call. = FALSE)
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

  set.seed(seed)
  km_res <- stats::kmeans(dr_mat, centers = k, ...)

  colData(fsce[[expt]])$k_cluster <- as.character(km_res$cluster)

  fsce
}

#' Run Leiden clustering algorithm
#'
#' Runs Leiden community detection algorithm to detect clusters.
#'
#' Execute [`install_py_deps()`] to install required python modules `leidenalg`
#' and `igraph`.
#'
#' @param fsce An object of class [`FunctionalSingleCellExperiment`]
#' @param expt Data to use for calculating variable features
#'   (default is `rnaseq`). Must be present in `names(fsce)`.
#' @param method dimensionality reduction method for clustering (defaults to
#'   PCA)
#' @param dims dimensions to use for nearest-neighbor calculation
#' @param prune Pruning parameter for shared nearest-neighbor calculation.
#' @param seed seed for `leidenalg$find_partition()`
#' @param ... Parameters to pass to the Python `leidenalg` function.
#'
#' @source <https://github.com/vtraag/leidenalg>
#' @source <https://github.com/TomKellyGenetics/leiden>
#'
#' @examples
#' \donttest{
#' fsce_small <- cluster_leiden(fsce_small)
#' SingleCellExperiment::colData(fsce_small[["rnaseq"]])
#' }
#'
#' @importFrom RANN nn2
#' @importFrom leiden leiden
#' @importFrom reticulate py_module_available
#'
#' @return fsce with `leiden_cluster` in `expt` colData.
#'
#' @family clustering functions
#'
#' @export
cluster_leiden <- function(fsce,
                           expt = "rnaseq",
                           method = "PCA",
                           dims = 1:5,
                           prune = 1/15,
                           seed = NULL,
                           partition_type = "ModularityVertexPartition",
                           ...){

  if (!reticulate::py_module_available("leidenalg") || !reticulate::py_module_available("igraph")) {
    warning("python modules and `leidenalg` and `igraph` are required", call. = FALSE)
    return(NULL)
  }

  ## check inputs
  if (!expt %in% names(fsce)) {
    stop(glue("expt `{expt}` not found in fsce "), call. = FALSE)
  }

  if (!method %in% names(reducedDims(fsce[[expt]]))) {
    stop(glue("method `{method}` not found in expt"), call. = FALSE)
  }

  mtx <- reducedDim(fsce[[expt]], method)[, dims]

  ## find nearest-neighbors in reduced dim space
  nn <- RANN::nn2(mtx)
  nn_idx <- nn$nn.idx

  adj_mat <- compute_snn_impl(nn_idx, prune)

  ## convert to integers and then matrix
  adj_mat <- as.matrix(ceiling(adj_mat))

  ## Run leidenalg
  parts <- leiden::leiden(
    adj_mat,
    partition_type = partition_type,
    seed = seed,
    ...
  )

  colData(fsce[[expt]])$leiden_cluster <- as.character(parts)

  fsce
}
