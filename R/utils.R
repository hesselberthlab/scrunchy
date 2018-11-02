#' Load 10x mtx matrix as sparseMatrix
#'
#' @param path path to 10x directory with mtx matrix
#' @param cell_prefix string to prefix to cell_id  (default = NULL)
#' @param strip_10x_suffix remove numeric suffix added by 10x
#' e.g. change "TTTGTCAAGGTGTGGT-1" to "TTTGTCAAGGTGTGGT" (default = TRUE)
#' @param use_gene_symbols If TRUE use gene symbols as row.names, if false then
#' gene_ids will be used (default = TRUE)
#'
#' @importFrom Matrix readMM
#' @importFrom stringr str_c str_remove
#' @importFrom readr read_tsv
#'
#' @export
read_10x_matrix <- function(path,
                            cell_prefix = NULL,
                            strip_10x_suffix = TRUE,
                            use_gene_symbols = TRUE) {

  fns <- dir(path)
  fns_needed <- c("barcodes.tsv", "genes.tsv", "matrix.mtx")

  if (!all(fns_needed %in% fns)) {
    stop(paste0("missing required 10x file: ", fns_needed[!fns_needed %in% fns], "\n"))
  }

  bcs <- readLines(file.path(path, "barcodes.tsv"))
  genes <- suppressMessages(readr::read_tsv(
    file.path(path, "genes.tsv"),
    col_names = c("gene_id", "gene_symbol")
  ))

  mat <- Matrix::readMM(file.path(path, "matrix.mtx"))

  if (strip_10x_suffix) {
    bcs <- stringr::str_remove(bcs, "-[0-9]+$")
  }

  if (!is.null(cell_prefix)) {
    bcs <- stringr::str_c(cell_prefix, "_", bcs)
  }

  if (use_gene_symbols) {
    genes_out <- genes[["gene_symbol"]]
  } else {
    genes_out <- genes[["gene_id"]]
  }

  genes_out <- make.unique(genes_out)
  rownames(mat) <- genes_out
  colnames(mat) <- bcs

  mat
}

