#' Normalize single cell rna and functional data
#' @param hcut_obj haircut object
#' @param rna_method normalization method for RNA data
#' @param hcut_method normalization method for functional data
#'
#' @return hcut_obj with log normalized counts returned in the logcounts slot
#' @export
normalize <- function(hcut_obj,
                      rna_method = "log_normalize",
                      hcut_method = "clr"){

  # check inputs
  rna_norm_methods <- c("log_normalize")
  hcut_norm_methods <- c("clr")
  if(!rna_method %in% rna_norm_methods){
    stop( "rna_method must be one of ", paste0(rna_norm_methods, collapse = ","))
  }

  if(!hcut_method %in% hcut_norm_methods){
    stop( "hcut_method must be one of ", paste0(hcut_norm_methods, collapse = ","))
  }

  if(rna_method == "log_normalize") {
    assays(hcut_obj[["rna"]])$logcounts <- log_normalize(counts(hcut_obj[["rna"]]))
  }

  if(hcut_method == "clr") {
    assays(hcut_obj[["hcut"]])$logcounts <- clr_normalize(counts(hcut_obj[["hcut"]]))
  }

  hcut_obj
}

#' Simple normalization for scrna-seq
#' @param mat input matrix
#' @param scaling_factor scalar to multiply normalized counts by to prevent small numbers (1e5)
#' @return matrix of normalized values. Normalization performed by dividing by column sums
#' (total counts per cell) and scaled by a scaling factor. Log values are returned with a pseudocount of 1
#' @importFrom Matrix colSums
log_normalize <- function(mat, constant = 1e5){
  mat <- constant * (mat / Matrix::colSums(mat))
  log1p(mat)
}


#' Simple normalization for functional data using centered log ratio
#' @param mat input matrix
#' @return matrix of normalized values. Normalization performed by dividing each function value by
#' the geometric mean of all functional values for a cell, and returning log values
clr_normalize <- function(mat){
  ## norm method from Seurat
  ## geom mean from https://stackoverflow.com/questions/2602583/geometric-mean-is-there-a-built-in
  apply(mat, 2, function(x) {
    log1p((x) /
            (exp(sum(log1p((x)[x > 0]), na.rm = TRUE) /
              length(x + 1))))
  })
}
