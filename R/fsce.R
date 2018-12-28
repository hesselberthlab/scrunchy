# Experiment types ------------------------------------------------------------

#' Create single-cell mRNA-seq experiment
#'
#' @param csv path to CSV counts matrix
#'
#' @return `SingleCellExperiment` containing a `sparseMatrix` of counts
#'
#' @examples
#' \dontrun {
#' expt_rnaseq(scrunchy_data("mrna.csv.gz"))
#' }
#'
#' @export
expt_rnaseq <- function(csv) {
  message(glue("Loading sc-rnaseq matrix: {csv}", csv = path_file(csv)))

  x <- load_matrix_csv(csv)
  x <- as(as.matrix(x), "sparseMatrix")

  sce <- SingleCellExperiment(assays = list(counts = x))

  int_metadata(sce) <- list(cells = colnames(x))

  colData(sce) <- DataFrame(
    row.names = colnames(x),
    sample_ids = extract_sample_ids(colnames(x))
  )

  sce
}

#' Create a single-cell Haircut experiment
#'
#' @param csv path to CSV counts matrix
#' @param adducts `data_frame` with positions of hairpin adducts. Expects
#'   two columns named `adduct` and `pos`.
#'
#' @return `SingleCellExperiment` containing a `matrix` of counts
#'
#' @examples
#' \dontrun {
#' expt_haircut(scrunchy_data("haircut.csv.gz"))
#' }
#'
#' @export
expt_haircut <- function(csv, adducts) {
  message(glue("Loading haircut matrix: {csv}", csv = path_file(csv)))

  x <- load_matrix_csv(csv)
  x <- as.matrix(x)

  hairpin_info <- strsplit(rownames(x), "_")

  if (length(hairpin_info[[1]]) != 2) {
    stop("hairpins must be named ADDUCTNAME_POS (e.g., Uracil_1)", call. = FALSE)
  }

  hairpin_id <- purrr::map(hairpin_info, 1)
  hairpin_pos <- purrr::map(hairpin_info, 2)

  sce <- SingleCellExperiment(
    assays = list(counts = x),
    rowData = data_frame(
      row.names = rownames(x),
      hairpin = hairpin_id,
      position = hairpin_pos
    )
  )

  int_metadata(sce) <- list(cells = colnames(x))
  colData(sce) <- DataFrame(
    row.names = colnames(x),
    sample_ids = extract_sample_ids(colnames(x))
  )

  sce
}

# fsce creation ------------------------------------------------------------

#' Create a functional single-cell experiment
#'
#' A functional single-cell experiment (`fsce`) is implemented as a
#' `MultiAssayExperiment` with one or more `SingleCellExperiment`s (`sce`). One
#' of these is a single-cell mRNA sequencing experiment that defines mRNA counts
#' for single cells. Other `sce`s include functional data from haircut.
#'
#' @param expt_list list mapping experiment names to `SingleCellExperiments` created
#'   by an `expt_` function.
#'
#' @return `MultiAssayExperiment` containing one or more `SingleCellExperiment`.
#'
#' @examples
#'
#' \dontrun{
#' fsce <- create_fsce(
#'   list(
#'     rnaseq  = expt_rnaseq(scrunchy_data("mrna.csv.gz")),
#'     haircut = expt_haircut(scrunchy_data("haircut.csv.gz"))
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

  MultiAssayExperiment(
    experiments = expt_list,
    colData = DataFrame(
        row.names = shared_cells,
        sample_id = sample_ids
      ),
    sampleMap = sample_map
  )
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
