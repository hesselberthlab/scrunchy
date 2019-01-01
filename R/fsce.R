# fsce creation ------------------------------------------------------------

#' Create a functional single-cell experiment
#'
#' A  [`FunctionalSingleCellExperiment`] is implemented as a
#' `MultiAssayExperiment` with one or more `SingleCellExperiment`s (`sce`). One
#' of these is a single-cell mRNA sequencing experiment that defines mRNA counts
#' for single cells. Other `sce`s include functional data from haircut.
#'
#' @param expt_list a list mapping experiment names to `SingleCellExperiments`
#'   created by an `create_` function.
#'
#' @return `MultiAssayExperiment` containing one or more `SingleCellExperiment`.
#'
#' @examples
#' # this is identical to the `fsce_small` data set:
#' create_fsce(
#'   list(
#'     rnaseq = create_sce_rnaseq(scrunchy_data("mrna/")),
#'     haircut = create_sce_haircut(scrunchy_data("haircut/"))
#'   )
#' )
#' @export
create_fsce <- function(expt_list) {
  expt_list <- ExperimentList(expt_list)

  shared_cells <- find_shared_cells(expt_list)

  if (length(shared_cells) == 0) {
    stop("No common cells between experiments.", call. = FALSE)
  }

  sample_map <- build_sample_map(expt_list)

  FunctionalSingleCellExperiment(
    experiments = expt_list,
    colData = DataFrame(
      row.names = shared_cells,
      cell_id = shared_cells
    ),
    sampleMap = sample_map
  )
}

# Classes -----------------------------------------------------------

#' @export
setClass("FunctionalSingleCellExperiment", contains = "MultiAssayExperiment")

# Constructors ----------------------------------------------

#' Constructor for a `FunctionalSingleCellExperiment` object
#'
#' A thin wrapper around `MultiAssayExperiment::MultiAssayExperiment`.
#'
#' @export
FunctionalSingleCellExperiment <- function(...) {
  sce <- MultiAssayExperiment(...)
  as(sce, "FunctionalSingleCellExperiment")
}

# Utilities ---------------------------------------------------------

build_sample_map <- function(expt_list) {
  cells <- purrr::map(as.list(expt_list), expt_cells)
  expts <- purrr::map(cells, build_cell_map)

  listToMap(expts)
}

extract_cell_ids <- function(x) {
  unlist(purrr::map(strsplit(x, "\\."), 1))
}

build_cell_map <- function(cells) {
  DataFrame(
    primary = cells,
    colname = cells
  )
}

find_shared_cells <- function(expt_list) {
  Reduce(intersect, purrr::map(as.list(expt_list), expt_cells))
}

expt_cells <- function(x) {
  x@int_metadata$cells
}
