#' Calculate principal components using irlba
#'
#' @param fsce [`FunctionalSingleCellExperiment`] object
#' @param expt Data to use for calculating variable features
#'   (default is `rnaseq`). Must be present in `names(fsce)`.
#' @param genes vector of genes to include in PCA (default is all genes)
#' @param n_pcs number of principle components to return
#' @param scale If `TRUE`, perform PCA on input on scaled data
#' @param seed seed for reproducible result
#'
#' @return `fsce` with PCA values added to reducedDims
#'
#' @export
calc_pca <- function(fsce,
                     expt = "rnaseq",
                     genes = NULL,
                     n_pcs = 20,
                     scale = TRUE,
                     seed = NULL) {

  ## check inputs
  if (!expt %in% names(fsce)) {
    stop(glue("expt `{expt}` not found in fsce"), call. = FALSE)
  }

  if (!"logcounts" %in% names(assays(fsce[[expt]]))) {
    stop("`logcounts` not found in fsce")
  }

  if (is.null(genes)) {
    genes <- rownames(fsce[[expt]])
  } else {
    genes <- genes[genes %in% rownames(fsce[[expt]])]
  }

  if (length(genes) == 0) {
    stop("input genes not found in fsce", call. = FALSE)
  }

  in_data <- logcounts(fsce[[expt]])[genes, ]

  ## remove rows without counts
  in_data <- in_data[Matrix::rowSums(in_data) > 0, ]

  if (scale) {
    message("scaling data")
    dr_mat <- scale(t(as.matrix(in_data)), center = TRUE, scale = TRUE)
  } else {
    dr_mat <- t(in_data)
  }

  message("calculating pcs")
  n_pcs <- min(c(n_pcs, dim(dr_mat) - 1))

  set.seed(seed)
  pcs <- irlba::prcomp_irlba(dr_mat,
    n = n_pcs,
    center = FALSE,
    scale. = FALSE
  )

  reducedDims(fsce[[expt]]) <- SimpleList(PCA = pcs$x)

  # store gene loadings and variance explained
  gene_loadings <- pcs$rotation
  rownames(gene_loadings) <- colnames(dr_mat)

  pc_out <- list(
    gene_loadings = gene_loadings,
    var_explained = pcs$sdev^2 / pcs$totalvar
  )

  metadata(fsce[[expt]])[["PCA"]] <- pc_out

  fsce
}

#' Generate 2D cell embeddings using UMAP
#'
#' @seealso See <https://umap-learn.readthedocs.io> for a detailed description of
#' parameters.
#'
#' @param fsce [`FunctionalSingleCellExperiment`] object
#' @param expt Data to use for calculating variable features
#'   (default is `rnaseq`). Must be present in `names(fsce)`.
#' @param method dimenality reduction method to use for UMAP (default is "PCA")
#' @param n_dims number of dimensions to pass to UMAP, defaults to all present
#'   in dr matrix
#' @param n_neighbors number of nearest neighbors to use for learning the
#'   manifold. Low values will preserve local structure, at the expense missing
#'   higher order organization. Higher values will capture more global structure
#'   but miss fine grained detail. Defaults to 30.
#' @param min_dist Numeric between 0 and 0.99. min_dist controls how tightly
#'   points can be packed together in 2D space. Lower values will generate more
#'   clumpy projections, but more accurately preserve local structure.
#' @param metric distance metric for UMAP, defaults to pearson.
#' @param seed seed to generate reproducible UMAP projection. Defaults to no
#'   seed.
#' @param ... additional arguments for [`umap::umap()`]
#'
#' @return `fsce` with UMAP values added to reducedDims
#'
#' @export
calc_umap <- function(fsce,
                      expt = "rnaseq",
                      method = "PCA",
                      n_dims = NULL,
                      n_neighbors = 30,
                      min_dist = 0.3,
                      metric = "euclidean",
                      seed = NA, # must be NA for umap defaults
                      ...) {

  ## check inputs
  if (!expt %in% names(fsce)) {
    stop(glue("expt `{expt}` not found in fsce"), call. = FALSE)
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

  umap_res <- umap::umap(
    dr_mat,
    method = "naive",
    metric = metric,
    random_state = seed,
    min_dist = min_dist,
    n_neighbors = n_neighbors,
    n_components = 2,
    ...
  )

  reducedDims(fsce[[expt]])$UMAP <- umap_res$layout

  fsce
}


#' Generate 2D cell embeddings using tSNE
#'
#' @seealso See [`Rtsne::Rtsne()`] for a detailed description of parameters.
#'
#'
#' @param fsce [`FunctionalSingleCellExperiment`] object
#' @param expt Data to use for calculating variable features
#'   (default is `rnaseq`). Must be present in `names(fsce)`.
#' @param method dimensionality reduction method for tSNE
#' @param n_dims number of dimensions for tSNE
#' @param perplexity tSNE perplexity value.
#' @param theta  tSNE theta value.
#' @param seed seed for generate reproducible tSNE projection.
#' @param ...  additional arguments for [`Rtsne::Rtsne()`]
#'
#' @return fce object with tSNE values added to reducedDims
#'
#' @export
calc_tsne <- function(fsce,
                      expt = "rnaseq",
                      method = "PCA",
                      n_dims = NULL,
                      perplexity = 30,
                      theta = 0.5,
                      seed = NULL,
                      ...) {

  ## check inputs
  if (!expt %in% names(fsce)) {
    stop(glue("expt `{expt}` not found in fsce"), call. = FALSE)
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

  set.seed(seed)
  tsne_res <- Rtsne::Rtsne(
    dr_mat,
    perplexity = perplexity,
    theta = theta,
    check_duplicates = FALSE,
    dims = 2,
    pca = FALSE,
    verbose = FALSE,
    pca_center = FALSE,
    ...
  )

  rownames(tsne_res$Y) <- rownames(dr_mat)

  reducedDims(fsce[[expt]])$TSNE <- tsne_res$Y

  fsce
}
