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
#' @param ... Parameters to pass to the Python `leidenalg` function.
#'
#' @source <https://github.com/vtraag/leidenalg>
#' @source <https://github.com/TomKellyGenetics/leiden>
#'
#' @examples
#' \donttest{
#' fsce_small <- cluster_leiden(fsce_small)
#' colData(fsce_small[["rnaseq"]])
#' }
#'
#' @importFrom RANN nn2
#' @importFrom igraph graph_from_adjacency_matrix write.graph
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
                           ...){
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

  snn_graph <- igraph::graph_from_adjacency_matrix(
    ceiling(as.matrix(adj_mat))
  )

  ## this pseudo-round-trip via an intermediate file is to convert
  ## the R object to a python object.
  ##
  ## XXX there is likely a way to eliminate the R igraph dependency here.
  ##
  tmp_graph <- fs::file_temp()
  igraph::write.graph(snn_graph, file=tmp_graph, format="graphml")
  snn_graph <- igraph_py$Graph$Read_GraphML(tmp_graph)

  ## Run leidenalg
  parts <- leidenalg_py$find_partition(
    snn_graph,
    leidenalg_py$ModularityVertexPartition,
    ...
  )

  colData(fsce[[expt]])$leiden_cluster <- as.character(parts$membership + 1)

  fsce
}

#' Install python dependencies for `cluster_leiden`.
#'
#' Users must run this prior to using `cluster_leiden()`.
#'
#' @param method method param for [`reticulate::py_install()`]
#' @param conda conda param for [`reticulate::py_install()`]
#'
#' @importFrom reticulate py_install
#' @export
install_py_deps <- function(method = "auto", conda = "auto") {
  reticulate::py_install("igraph", method = method, conda = conda)
  reticulate::py_install("leidenalg", method = method, conda = conda)
}
