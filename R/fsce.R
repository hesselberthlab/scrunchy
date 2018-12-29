# fsce creation ------------------------------------------------------------

#' Create a functional single-cell experiment
#'
#' A functional single-cell experiment (`fsce`) is implemented as a
#' `MultiAssayExperiment` with one or more `SingleCellExperiment`s (`sce`). One
#' of these is a single-cell mRNA sequencing experiment that defines mRNA counts
#' for single cells. Other `sce`s include functional data from haircut.
#'
#' @param expt_list a list mapping experiment names to `SingleCellExperiments`
#'   created by an `sce_` function.
#'
#' @return `MultiAssayExperiment` containing one or more `SingleCellExperiment`.
#'
#' @examples
#'
#' \dontrun{
#' fsce <- create_fsce(
#'   list(
#'     rnaseq  = sce_rnaseq(scrunchy_data("mrna.csv.gz")),
#'     haircut = sce_haircut(scrunchy_data("haircut.csv.gz"))
#'   )
#' )
#'
#' fsce
#' }
#'
#' @export
create_fsce <- function(expt_list) {

  expt_list <- ExperimentList(expt_list)

  shared_cells <- find_shared_cells(expt_list)

  if (length(shared_cells) == 0) {
    stop("No common cells between experiments.", call. = FALSE)
  }

  sample_map <- build_sample_map(expt_list)
  sample_ids <- extract_sample_ids(shared_cells)

  FunctionalSingleCellExperiment(
    experiments = expt_list,
    colData = DataFrame(
        row.names = shared_cells,
        sample_id = sample_ids
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

extract_sample_ids <- function(x) {
  unlist(purrr::map(strsplit(x, "-"), 1))
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
