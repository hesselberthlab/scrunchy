#' Get variable genes
#'
#'
#' @param fce fce object
#' @param expt Data to use for calculating variable genes,
#'  one of either sce (rna data, the default) or fsce (functional data)
#' @param n_genes n variable genes to return
#'
#' @return Character vecctor
#'
#' @export
get_var_genes <- function(fce,
                          expt = "sce",
                          n_genes = 1000) {

  ## check inputs
  if (!expt %in% names(fce)) {
    stop("expt not found in fce object")
  }

  if (!"logcounts" %in% names(assays(fce[[expt]]))) {
    stop("logcounts not found in fce object, run normalize_counts first")
  }
  log_counts <- assay(fce[[expt]], "logcounts")
  n_genes <- min(n_genes, nrow(log_counts))

  rv <- apply(log_counts, 1, var)

  top_n_rv <- sort(rv, decreasing = T)[1:n_genes]

  names(top_n_rv)
}
