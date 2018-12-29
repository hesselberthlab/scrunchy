#' Retrieve reducedDim data
#'
#' Subset is returned in a tidied, tibble format, suitable for plotting or
#' analysis.
#'
#' Cell IDs are added in column `cell_id` to facilitate joining of subsets.
#'
#' @param fsce [`FunctionalSingleCellExperiment`]
#' @param expt Data to use for calculating variable features
#'   (default is `rnaseq`). Must be present in `names(fsce)`.
#' @param dim_name name of reducedDim to retrieve (default: UMAP)
#' @param dims dimensions to retrieve (default: `c(1, 2)`)
#'
#' @return Data in `tibble` format.
#'
#' @examples
#' subset_dim(fsce_small, expt = "rnaseq", dim_name = "UMAP")
#'
#' subset_dim(fsce_small, expt = "rnaseq", dim_name = "TSNE")
#'
#' subset_dim(fsce_small, expt = "rnaseq", dim_name = "PCA", dims = c(1:3))
#'
#' @export
subset_dim <- function(fsce,
                       expt = "rnaseq",
                       dim_name = "UMAP",
                       dims = c(1, 2)
) {

  if (!expt %in% names(fsce)) {
    stop(glue("expt `{expt}` not found in fsce"), call. = FALSE)
  }

  x <- fsce[[expt]]

  if (!dim_name %in% reducedDimNames(x)) {
    stop(glue("dim names `{dim_name}` not found in expt `{expt}`"), call. = FALSE)
  }

  res <- reducedDim(x, dim_name)[, dims]

  colnames(res) <- paste0(dim_name, dims)

  res <- rownames_to_column(as.data.frame(res), "cell_id")

  as_tibble(res)
}


#' Retrieve counts from an fsce
#'
#' Subset is returned in a tidied, tibble format, suitable for plotting or analysis.
#'
#' Cell IDs are added in column `cell_id` to facilitate joining of subsets.
#'
#' @param fsce [`FunctionalSingleCellExperiment`]
#' @param expt Data to use for calculating variable features
#'   (default is `rnaseq`). Must be present in `names(fsce)`.
#' @param vars variables to retrieve (default: all)
#'
#' @examples
#' subset_counts(fsce_small, expt = "haircut", vars = c("Uracil_45"))
#'
#' @export
subset_counts <- function(fsce, expt, vars = NULL) {

  if (!expt %in% names(fsce)) {
    stop(glue("expt `{expt}` not found in fsce"), call. = FALSE)
  }

  x <- fsce[[expt]]

  if (is.null(vars)) {
    res <- as.matrix(t(counts(x)))
  } else {
    res <- as.matrix(t(counts(x)[vars, , drop = FALSE]))
  }

  res <- rownames_to_column(as.data.frame(res), "cell_id")

  as_tibble(res)
}

#' Retrieve logcounts from an fsce
#'
#' Subset is returned in a tidied, tibble format, suitable for plotting or analysis.
#'
#' Cell IDs are added in column `cell_id` to facilitate joining of subsets.
#'
#' @inheritParams subset_counts
#'
#' @examples
#' subset_logcounts(fsce_small, expt = "haircut", vars = c("Uracil_45"))
#'
#' @export
subset_logcounts <- function(fsce, expt, vars = NULL) {
  if (!expt %in% names(fsce)) {
    stop(glue("expt `{expt}` not found in fsce"), call. = FALSE)
  }

  x <- fsce[[expt]]

  if (is.null(vars)) {
    res <- as.matrix(t(logcounts(x)))
  } else {
    res <- as.matrix(t(logcounts(x)[vars, , drop = FALSE]))
  }

  res <- rownames_to_column(as.data.frame(res), "cell_id")

  as_tibble(res)
}

#' Retrieve colData from an fsce
#'
#' Subset is returned in a tidied, tibble format, suitable for plotting or analysis.
#'
#' Cell IDs are added in column `cell_id` to facilitate joining of subsets.
#'
#' @param fsce [`FunctionalSingleCellExperiment`]
#' @param expt Data to use for calculating variable features
#'   (default is `rnaseq`). Must be present in `names(fsce)`.
#' @param col_names column names to retrieve
#'
#' @return Data in `tibble` format.
#'
#' @examples
#' subset_col(fsce_small, expt = "haircut", col_names = c("Uracil_45"))
#'
#' @export
subset_col <- function(fsce,
                       expt = "rnaseq",
                       col_names = NULL) {

  if (!expt %in% names(fsce)) {
    stop(glue("expt `{expt}` not found in fsce"), call. = FALSE)
  }

  x <- fsce[[expt]]

  if (!col_names %in% colnames(colData(fsce[[expt]]))) {
    stop(glue("colname `{col_name}` not found in expt `{expt}`"), call. = FALSE)
  }

  res <- colData(x, method)

  if (!is.null(cell_ids)) {
    res <- cell_subset(x, cell_ids)
  }

  as_tibble(res)
}

#' Filter a data subset for specific cells
#'
#' Expects a column named `cell_id` in the the data subset.
#'
#' @param x a data subset
#' @param cell_ids vector of cell ids
#'
#' @examples
#' x <- subset_dim(fsce_small, expt = "rnaseq", dim_name = "UMAP")
#'
#' filter_cells(x, c("TGCGGGTTCAACCATG.1", "TGGTTCCTCACTTACT.1"))
#'
#' @export
filter_cells <- function(x, cell_ids) {
  if (!"cell_id" %in% names(x)) {
    stop("expected variable `cell_id` in data subset", call. = FALSE)
  }
  semi_join(x, data_frame(cell_id = cell_ids), by = "cell_id")
}
