#' Get variable features
#'
#' Identify features with the most variance in the data set.
#'
#' @param fsce fsce object
#' @param expt Data to use for calculating variable features
#'   (default is `rnaseq`). Must be present in `names(fsce)`.
#' @param n number of variable features to return
#'
#' @return Names of variable features
#'
#' @export
get_var_features <- function(fsce, expt = "rnaseq", n = 1000) {

  ## check inputs
  if (!expt %in% names(fsce)) {
    stop(glue("expt `{expt}` not found in fsce"), call. = FALSE)
  }

  if (!"logcounts" %in% names(assays(fsce[[expt]]))) {
    stop(glue("`logcounts` not found for expt `{expt}`"), call. = FALSE)
  }

  logcounts <- logcounts(fsce[[expt]])
  n <- min(n, nrow(logcounts))

  rv <- apply(logcounts, 1, var)

  res <- sort(rv, decreasing = T)[1:n]

  names(res)
}
