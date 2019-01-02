# Tidiers -----------------------------------------------------------

#' Tidy logcounts data
#'
#' @section Tidying:
#' This is a data tidier for a `SingleCellExperient`. Returns
#' data from an `sce` in a tidy format, where variables are columns and
#' observations are rows.
#'
#' @param sce An object of class [`SingleCellExperiment::SingleCellExperiment`].
#'
#' @examples
#' x <- fsce_small[ c("Uracil_45"), , "haircut"]
#' tidy_logcounts(x)
#'
#' @family tidiers
#'
#' @export
tidy_logcounts <- function(sce) {
  es <- as.list(experiments(sce))
  cs <- purrr::map(es, logcounts)

  purrr::map_dfr(cs, counts_tbl, .id = "experiment")
}

#' Tidy counts data
#'
#' @inheritSection tidy_logcounts Tidying
#' @inheritParams tidy_logcounts
#
#' @examples
#' x <- fsce_small[ c("Uracil_45"), , "haircut"]
#' tidy_counts(x)
#'
#' @family tidiers
#'
#' @export
tidy_counts <- function(sce) {
  es <- as.list(experiments(sce))
  cs <- purrr::map(es, counts)

  purrr::map_dfr(cs, counts_tbl, .id = "experiment")
}

#' Tidy reducedDims data
#'
#' @inheritSection tidy_logcounts Tidying
#' @inheritParams tidy_logcounts
#' @param dimnames vector of dimnames to retrieve. If `NULL`, retrieve all
#'   dimnames.
#' @param dims vector of dimensions to retrieve
#'
#' @examples
#' tidy_dims(fsce_small)
#'
#' tidy_dims(fsce_small, dimnames = c("UMAP"))
#'
#' @family tidiers
#'
#' @export
tidy_dims <- function(sce, dimnames = NULL, dims = c(1, 2)) {
  es <- as.list(experiments(sce))

  res <- purrr::map(es, dims_tbl, dimnames, dims)
  unframe(res, name = "experiment")
}

#' Tidy colData
#'
#' @inheritSection tidy_logcounts Tidying
#' @inheritParams tidy_logcounts
#
#' @examples
#' tidy_coldata(fsce_small)
#'
#' @family tidiers
#'
#' @export
tidy_coldata <- function(sce) {
  es <- as.list(experiments(sce))
  cds <- purrr::map(es, colData)
  cds <- purrr::map(cds, as.data.frame)

  res <- purrr::reduce(cds, left_join, by = "cell_id")
  as_tibble(res)
}

# Utilities -----------------------------------------------------------

counts_tbl <- function(mx) {
  mx <- as.matrix(t(mx))
  res <- rownames_to_column(as.data.frame(mx), "cell_id")
  as_tibble(res)
}

dims_tbl <- function(expt, dimnames = NULL, dims = c(1, 2)) {
  rds <- as.list(reducedDims(expt))
  if (length(rds) == 0) {
    return(tibble())
  }

  if (!is.null(dimnames)) {
    rds <- rds[dimnames]
  }

  res <- purrr::imap(rds, dim_tbl, dims)
  purrr::reduce(res, left_join, by = "cell_id")
}

dim_tbl <- function(rd, dimname, dims) {
  if (!is.null(dims)) {
    rd <- rd[, dims]
  } else {
    # default is all dims
    dims <- 1:ncol(rd)
  }

  colnames(rd) <- paste0(dimname, dims)
  rd <- rownames_to_column(as.data.frame(rd), "cell_id")

  as_tibble(rd)
}

#' @noRd
#' @importFrom tidyr unnest
unframe <- function(x, name = "name") {
  unnest(enframe(x, name = name))
}
