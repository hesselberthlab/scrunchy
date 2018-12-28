#' Get variable genes
#'
#' Identify genes the most variance in the data set.
#'
#' @param fce fce object
#' @param expt Data to use for calculating variable genes,
#'  one of either `sce` (the default) or `fsce` (functional data)
#' @param n_genes number of variable genes to return
#'
#' @return Names of variable genes
#'
#' @export
get_var_genes <- function(fce, expt = "sce", n_genes = 1000) {

  ## check inputs
  if (!expt %in% names(fce)) {
    stop("`expt` not found in fsce")
  }

  if (!"logcounts" %in% names(assays(fce[[expt]]))) {
    stop("`logcounts` not found in fsce, run `normalize_counts()`")
  }

  log_counts <- assay(fce[[expt]], "logcounts")
  n_genes <- min(n_genes, nrow(log_counts))

  rv <- apply(log_counts, 1, var)

  top_n_rv <- sort(rv, decreasing = T)[1:n_genes]

  names(top_n_rv)
}
