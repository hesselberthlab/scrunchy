# Experiment types ------------------------------------------------------------

#' Create a single-cell mRNA-seq experiment
#'
#' @param path path to output matrices or R matrix object.
#' @param norm_method Normalization method for `counts`. Normalized data
#'   is stored in `logcounts`. Set to `NULL` to skip normalization.
#'
#' @return `SingleCellExperiment` containing a `sparseMatrix` of counts
#'
#' @examples
#' create_sce_rnaseq(scrunchy_data("mrna/"))
#'
#' # using pre-loaded matrix
#' mat <- SingleCellExperiment::counts(fsce_small[["rnaseq"]])
#' create_sce_rnaseq(mat)
#'
#' @export
create_sce_rnaseq <- function(path, norm_method = "log_normalize") {

  if(is_any_matrix(path)) {
    x <- path
  } else {
    message(glue("Loading sc-rnaseq matrix files: {path}"))
    x <- read_matrix(path)
  }

  sce <- SingleCellExperiment::SingleCellExperiment(assays = list(counts = x))

  cell_ids <- extract_cell_ids(colnames(x))

  int_metadata(sce)$cells <- cell_ids

  colData(sce) <- DataFrame(
    row.names = cell_ids,
    cell_id = cell_ids
  )

  if (!is.null(norm_method)) {
    logcounts(sce) <- normalize(sce, method = norm_method)
  }

  sce
}

#' Create a single-cell Haircut experiment
#'
#' @param path path to output matrices, or R matrix object.
#' @param norm_method Normalization method for `counts`. Normalized data
#'   is stored in `logcounts`. Set to `NULL` to skip normalization.
#' @param adducts `data_frame` with positions of hairpin adducts. Expects
#'   two columns named `adduct` and `pos`.
#'
#' @return `SingleCellExperiment` containing a `sparseMatrix` of counts

#' @examples
#' create_sce_haircut(scrunchy_data("haircut/"))
#'
#' # using pre-loaded matrix
#' mat <- SingleCellExperiment::counts(fsce_small[["haircut"]])
#' create_sce_haircut(mat)
#'
#' @export
create_sce_haircut <- function(path, norm_method = "clr_normalize", adducts = NULL) {

  if(is_any_matrix(path)) {
    x <- path
  } else {
    message(glue("Loading haircut matrix files: {path}", path = path))
    x <- read_matrix(path)
  }


  hairpin_info <- strsplit(rownames(x), "_")

  if (length(hairpin_info[[1]]) != 2) {
    stop("hairpins must be named ADDUCTNAME_POS (e.g., Uracil_1)", call. = FALSE)
  }

  hairpin_id <- unlist(purrr::map(hairpin_info, 1))
  hairpin_pos <- unlist(purrr::map(hairpin_info, 2))

  sce <- SingleCellExperiment::SingleCellExperiment(
    assays = list(counts = x),
    rowData = DataFrame(
      hairpin = hairpin_id,
      position = hairpin_pos
    )
  )

  cell_ids <- extract_cell_ids(colnames(x))

  int_metadata(sce)$cells <- cell_ids

  colData(sce) <- DataFrame(
    row.names = cell_ids,
    cell_id = cell_ids
  )

  if (!is.null(norm_method)) {
    logcounts(sce) <- normalize(sce, method = norm_method)
  }

  sce
}

#' Create a single-cell feature experiment
#'
#' @param path path to output matrices, or R matrix object.
#' @param norm_method Normalization method for `counts`. Normalized data
#'   is stored in `logcounts`. Set to `NULL` to skip normalization.
#'
#' @return `SingleCellExperiment` containing a `sparseMatrix` of counts
#'
#' @export
create_sce_feature <- function(path, norm_method = "clr_normalize") {

  if(is_any_matrix(path)) {
    x <- path
  } else {
    message(glue("Loading feature matrix files: {path}", path = path))
    x <- read_matrix(path)
  }

  sce <- SingleCellExperiment::SingleCellExperiment(
    assays = list(counts = x),
    rowData = DataFrame(
      id = rownames(x)
    )
  )

  cell_ids <- extract_cell_ids(colnames(x))

  int_metadata(sce)$cells <- cell_ids

  colData(sce) <- DataFrame(
    row.names = cell_ids,
    cell_id = cell_ids
  )

  if (!is.null(norm_method)) {
    logcounts(sce) <- normalize(sce, method = norm_method)
  }

  sce
}

# Generics --------------------------------------------------

#' Normalize data in a SingleCellExperiment
#'
#' @export
setGeneric("normalize", function(x, method, ...) standardGeneric("normalize"))

#' @rdname normalize
#'
#' @param x matrix to normalize
#' @param method normalization method
#' @param ... addition arguments for normalization method
#'
#' @export
setMethod(
  "normalize",
  signature("SingleCellExperiment"),
  function(x, method = "log_normalize", ...) {
    do.call(method, list(counts(x), ...))
  }
)

# Utilities --------------------------------------------------
#' @importFrom methods is
is_sparse <- function(x) {
  is(x, "dgTMatrix")
}

is_any_matrix <- function(x) {
  is_sparse(x) || is.matrix(x)
}
