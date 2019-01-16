#' Calculate variable features
#'
#' Identify features with the highest variance in the data set.
#'
#' @param fsce [`FunctionalSingleCellExperiment`]
#' @param expt Data to use for calculating variable features
#'   (default is `rnaseq`). Must be present in `names(fsce)`.
#' @param n number of variable features to return
#'
#' @return Names of variable features
#'
#' @examples
#' calc_var_features(fsce_small)[1:5]
#'
#' calc_var_features(fsce_small, expt = "haircut")[1:5]
#' @export
calc_var_features <- function(fsce, expt = "rnaseq", n = 1000) {

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

  res <- sort(rv, decreasing = TRUE)[1:n]

  names(res)
}

#' Calculate cell-cycle assignments
#'
#' Uses [`scran::cyclone()`] to assign cell-cycle phases to each cell.
#' Assignments are available in `colData` in the `cell_cycle` slot.
#'
#' @param fsce [`FunctionalSingleCellExperiment`]
#' @param expt Data to use for calculating variable features
#'   (default is `rnaseq`). Must be present in `names(fsce)`.
#' @param org organism for gene pairs, human or mouse
#' @param ... additional params for [`scran::cyclone()`]
#'
#' @examples
#' fsce_small <- calc_cell_cycle(fsce_small)
#'
#' SingleCellExperiment::colData(fsce_small[["rnaseq"]])
#'
#' @importFrom scran cyclone
#'
#' @export
calc_cell_cycle <- function(fsce, expt = "rnaseq", org = "human", ...) {

  ## check inputs
  if (!expt %in% names(fsce)) {
    stop(glue("expt `{expt}` not found in fsce"), call. = FALSE)
  }

  if (!"counts" %in% names(assays(fsce[[expt]]))) {
    stop(glue("`counts` not found for expt `{expt}`"), call. = FALSE)
  }

  gene_pairs <- readr::read_rds(
    system.file(
      "exdata",
      paste0(org, "_cycle_markers.rds"),
      package="scran",
      mustWork = TRUE
    )
  )

  ## convert to ensembl ids
  mtx <- as.matrix(counts(fsce[[expt]]))
  gene_syms <- tibble(gene_symbol = rownames(mtx))

  ens_ids <- left_join(gene_syms, human_gene_ids, by = "gene_symbol")
  rownames(mtx) <- ens_ids$ensembl_id

  res <- scran::cyclone(mtx, gene_pairs, ...)
  colData(fsce[[expt]])$cell_cycle <- res$phases

  fsce
}
