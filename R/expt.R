# Experiment types ------------------------------------------------------------

#' Create a single-cell mRNA-seq experiment
#'
#' @param csv path to CSV counts matrix
#' @param norm_method Normalization method for `counts`. Normalized data
#'   is stored in `logcounts`. Set to `NULL` to skip normalization.
#'
#' @return `SCERnaSeq` containing a `sparseMatrix` of counts
#'
#' @examples
#' \dontrun {
#' sce_rnaseq(scrunchy_data("mrna.csv.gz"))
#' }
#'
#' @export
sce_rnaseq <- function(csv, norm_method = "log_normalize") {
  message(glue("Loading sc-rnaseq matrix: {csv}", csv = path_file(csv)))

  x <- load_matrix_csv(csv)
  x <- as(as.matrix(x), "sparseMatrix")

  sce <- SCERnaSeq(assays = list(counts = x))

  int_metadata(sce) <- list(cells = colnames(x))

  colData(sce) <- DataFrame(
    row.names = colnames(x),
    sample_ids = extract_sample_ids(colnames(x))
  )

  if (!is.null(norm_method)) {
    logcounts(sce) <- normalize(sce, method = norm_method)
  }

  sce
}

#' Create a single-cell Haircut experiment
#'
#' @param csv path to CSV counts matrix
#' @param norm_method Normalization method for `counts`. Normalized data
#'   is stored in `logcounts`. Set to `NULL` to skip normalization.
#' @param adducts `data_frame` with positions of hairpin adducts. Expects
#'   two columns named `adduct` and `pos`.
#'
#' @return `SCEHaircut` containing a `matrix` of counts
#'
#' @examples
#' \dontrun {
#' sce_haircut(scrunchy_data("haircut.csv.gz"))
#' }
#'
#' @export
sce_haircut <- function(csv, norm_method = "clr_normalize", adducts = NULL) {
  message(glue("Loading haircut matrix: {csv}", csv = path_file(csv)))

  x <- load_matrix_csv(csv)
  x <- as.matrix(x)

  hairpin_info <- strsplit(rownames(x), "_")

  if (length(hairpin_info[[1]]) != 2) {
    stop("hairpins must be named ADDUCTNAME_POS (e.g., Uracil_1)", call. = FALSE)
  }

  hairpin_id <- unlist(purrr::map(hairpin_info, 1))
  hairpin_pos <- unlist(purrr::map(hairpin_info, 2))

  sce <- SCEHaircut(
    assays = list(counts = x),
    rowData = DataFrame(
      hairpin = hairpin_id,
      position = hairpin_pos
    )
  )

  int_metadata(sce) <- list(cells = colnames(x))

  colData(sce) <- DataFrame(
    row.names = colnames(x),
    sample_ids = extract_sample_ids(colnames(x))
  )

  if (!is.null(normalize)) {
    logcounts(sce) <- normalize(sce, method = norm_method)
  }

  sce
}

# Classes ---------------------------------------------------

#' @export
setClass("SCERnaSeq", contains = "SingleCellExperiment")

#' @export
setClass("SCEHaircut", contains = "SingleCellExperiment")

# Constructors ----------------------------------------------

#' Constructor for a `SCERnaSeq` object
#'
#' @export
SCERnaSeq <- function(...) {
  sce <- SingleCellExperiment(...)
  as(sce, "SCERnaSeq")
}

#' Constructor for a `SCEHaircut` object
#'
#' @export
SCEHaircut <- function(...) {
  sce <- SingleCellExperiment(...)
  as(sce, "SCEHaircut")
}

# Generics --------------------------------------------------

#' @export
setGeneric("normalize", function(x, method, ...) standardGeneric("normalize"))

setMethod("normalize",
  signature("SCERnaSeq"),
  function(x, method = "log_normalize", ...) {
    logcounts(x) <- do.call(method, list(counts(x), ...))
  }
)

setMethod("normalize",
  signature("SCEHaircut"),
  function(x, method = "clr_normalize", ...) {
    logcounts(x) <- do.call(method, list(counts(x), ...))
  }
)
