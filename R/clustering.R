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
#' @importFrom reticulate r_to_py
#' @importFrom RANN nn2
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

  ## convert to integers and then matrix
  adj_mat <- as.matrix(ceiling(adj_mat))

  ## convert to python numpy.ndarray, then list
  adj_mat_py <- reticulate::r_to_py(adj_mat)
  adj_mat_py <- adj_mat_py$tolist()

  snn_graph <- igraph_py$Graph$Adjacency(adj_mat_py)

  if (!is.null(seed)) {
    seed <- r_to_py(as.integer(seed))
  } else {
    seed <- r_to_py(NULL)
  }

  ## Run leidenalg
  parts <- leidenalg_py$find_partition(
    graph = snn_graph,
    partition_type = leidenalg_py$ModularityVertexPartition,
    seed = seed,
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
  reticulate::py_install("python-igraph", method = method, conda = conda)
  reticulate::py_install("leidenalg", method = method, conda = conda)
}

#' Adds useful labels to colData.
#' Can be used to add cell tyles to cluster numbers
#'
#' @param fsce An object of class [`FunctionalSingleCellExperiment`]
#' @param expt Data to use for calculating variable features
#'   (default is `rnaseq`). Must be present in `names(fsce)`.
#' @param from vector of labels present in colData of `expt`
#' @param to vector of new labels to add to colData of `expt`
#' @param match column name in colData of `expt` where `from` values are found
#' @param new_label new column name of colData of `expt` where `to` values should go
#' @param expt Data to use for match labels
#'   (default is `rnaseq`). Must be present in `names(fsce)`.
#'
#' @return fsce with `new_label` in `expt` colData.
#'
#' @examples
#' # Add cell_type labels to PBMC data
#' labels <- tibble::tribble(
#'   ~ k_cluster, ~ label,
#'   1, "MC",
#'   2, "NK",
#'   3, "NK+T",
#'   4, "MC",
#'   5, "MK",
#'   6, "CD4/8 T")  %>% mutate(k_cluster = as.character(k_cluster))
#' fsce <- add_labels(fsce_small, labels$k_cluster, labels$label)
#'
#' colData(fsce[["rnaseq]])
#'
#' @family clustering functions
#'
#' @export

add_label <- function(fsce, from, to,
                      match = "k_cluster",
                      new_label = "cell_type",
                      expt = 'rnaseq') {

  colData(fsce[[expt]])[new_label] <- plyr::mapvalues(colData(fsce[[expt]])[[match]],
                                                      from,
                                                      to)
  fsce
}
