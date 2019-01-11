#' Log-normalize a matrix
#'
#' Normalization performed by dividing by column sums (total counts per cell)
#' and scaled by a constant. Log values are returned with a pseudocount of 1.
#'
#' @param mat input matrix
#' @param constant multiply normalized counts by this scalar to prevent small
#'   numbers
#'
#' @family normalization methods
#'
#' @return matrix of normalized values.
#'
log_normalize <- function(mat, constant = 1e4) {
  mat <- constant * (sweep(mat, 2, Matrix::colSums(mat), "/"))
  log1p(mat)
}


#' Centered log-ratio normalize a matrix
#'
#' Normalization performed by dividing each value by the geometric mean of all
#' values for a cell, returning log values.
#'
#' This normalization strategy is used for CITE-seq and other feature data in
#' [Seurat](https://satijalab.org/seurat/).
#'
#' @param mat input matrix
#'
#' @family normalization methods
#'
#' @seealso <https://stackoverflow.com/questions/2602583/geometric-mean-is-there-a-built-in>
#'
#' @return matrix of normalized values.
#'
clr_normalize <- function(mat) {
  apply(mat, 2, function(x) {
    log1p((x) /
      (exp(sum(log1p((x)[x > 0]), na.rm = TRUE) /
        length(x + 1))))
  })
}
