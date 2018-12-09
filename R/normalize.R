#' Normalize single cell rna and functional data
#'
#' @param fce fce object generted by [`create_fce()`]
#' @param rna_method normalization method for RNA data
#' @param functional_method normalization method for functional data
#'
#' @return fce object with log normalized counts returned in the logcounts slot
#'
#' @export
normalize_counts <- function(fce,
                             rna_method = "log_normalize",
                             functional_method = "clr") {

  # check inputs
  rna_norm_methods <- c("log_normalize")
  f_norm_methods <- c("clr")
  if (!rna_method %in% rna_norm_methods) {
    stop("rna_method must be one of ", paste0(rna_norm_methods, collapse = ","))
  }

  if (!functional_method %in% f_norm_methods) {
    stop("functional_method must be one of ", paste0(f_norm_methods, collapse = ","))
  }

  if (rna_method == "log_normalize") {
    assays(fce[["sce"]])$logcounts <- log_normalize(counts(fce[["sce"]]))
  }

  if (functional_method == "clr") {
    assays(fce[["fsce"]])$logcounts <- clr_normalize(counts(fce[["fsce"]]))
  }

  fce
}

#' Log normalization
#'
#' @param mat input matrix
#' @param constant multiply normalized counts by this scalar to prevent
#'   small numbers
#'
#' @return matrix of normalized values. Normalization performed by dividing by
#'   column sums (total counts per cell) and scaled by a constant. Log
#'   values are returned with a pseudocount of 1.
#'
log_normalize <- function(mat, constant = 1e4) {
  mat <- constant * (sweep(mat, 2, Matrix::colSums(mat), "/"))
  log1p(mat)
}


#' Centered log-ratio normalization
#'
#' This normalization strategy is used for CITE-seq and other feature data in
#' [Seurat](https://satijalab.org/seurat/).
#'
#' @param mat input matrix
#'
#' @seealso https://stackoverflow.com/questions/2602583/geometric-mean-is-there-a-built-in
#'
#' @return matrix of normalized values. Normalization performed by dividing each
#'   value by the geometric mean of all values for a cell, returning log values.
#'
clr_normalize <- function(mat) {
  apply(mat, 2, function(x) {
    log1p((x) /
      (exp(sum(log1p((x)[x > 0]), na.rm = TRUE) /
        length(x + 1))))
  })
}
