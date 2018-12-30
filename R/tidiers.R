# Tidiers -----------------------------------------------------------

#' Tidy logcounts data
#'
#' @section Tidying:
#' This method is a data tidier for a [`FunctionalSingleCellExperient`]. Returns
#' data from an `fsce` in a tidy format, where variables are columns and
#' observations are rows. Output will contain `experiment` and `cell_id` columns
#' that describe the source of the data.
#'
#' @param fsce A [`FunctionalSingleCellExperiment`].
#'
#' @examples
#' x <- fsce_small[ c("Uracil_45"), , "haircut"]
#' tidy_logcounts(x)
#'
#' @export
tidy_logcounts <- function(fsce) {
  es <- as.list(experiments(fsce))
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
#' @export
tidy_counts <- function(fsce) {
  es <- as.list(experiments(fsce))
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
#' @export
tidy_dims <- function(fsce, dimnames = NULL, dims = c(1,2)) {
  es <- as.list(experiments(fsce))

  res <- purrr::map(es, dims_tbl, dimnames, dims)
  purrr::map_dfr(res, unframe, .id = "experiment")
}


# Utilities -----------------------------------------------------------

counts_tbl <- function(mx) {
  mx <- as.matrix(t(mx))
  res <- rownames_to_column(as.data.frame(mx), "cell_id")
  as_tibble(res)
}

dims_tbl <- function(expt, dimnames = NULL, dims = c(1, 2)) {
  rds <- as.list(reducedDims(expt))

  if (!is.null(dimnames)) {
    pred <- purrr::map_lgl(dims, ~ .x %in% reducedDimNames(expt))
    rds <- rds[pred]
  }

  purrr::imap(rds, dim_tbl, dims)
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

unframe <- function(x) {
  unnest(enframe(x))
}