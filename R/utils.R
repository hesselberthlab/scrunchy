
# Utilities ---------------------------------------------------------

#' Provide path to scrunchy internal data
#'
#' @param x file name
#'
#' @export
scrunchy_data <- function(x) {
  system.file("extdata", x, package = "scrunchy", mustWork = TRUE)
}

#' Load mtx formatted 10x or haircut matrices
#'
#' @param path path to directory with mtx matrix
#' @param cell_prefix string to prefix to cell_id  (default = NULL)
#' @param strip_10x_suffix remove numeric suffix added by 10x e.g. change
#'   "TTTGTCAAGGTGTGGT-1" to "TTTGTCAAGGTGTGGT" (default = TRUE)
#' @param use_gene_symbols If TRUE use gene symbols as row.names, if false then
#'   gene_ids will be used (default = TRUE)
#' @param matrix_fn filename for matrix
#' @param features_fn filename for genes file
#' @param barcodes_fn filename for barcodes file
#'
#' @export
read_matrix <- function(path,
                        cell_prefix = NULL,
                        strip_10x_suffix = TRUE,
                        use_gene_symbols = TRUE,
                        matrix_fn = "matrix.mtx.gz",
                        features_fn = "features.tsv.gz",
                        barcodes_fn = "barcodes.tsv.gz") {
  filenames <- list(
    matrix = matrix_fn,
    features = features_fn,
    barcodes = barcodes_fn
  )

  filenames <- map(filenames, ~ path_join(c(path, .x)))

  files_exist <- file_exists(unlist(filenames))

  if (!all(files_exist)) {
    # check for gzipped equivalents
    gzipped_fns <- paste0(filenames, ".gz")
    is_gzipped <- file_exists(gzipped_fns)

    # rename if found
    if (any(is_gzipped)) {
      filenames[is_gzipped] <- gzipped_fns[is_gzipped]
    } else {
      stop(paste(
        "missing required files:",
        filenames[!files_exist], "\n"
      ),
      call. = FALSE
      )
    }
  }

  # assign column names based on feature file type
  n_fcols <- count_cols(filenames$features)
  if (n_fcols == 3) {
    col_args <- fcols_10x_v3
  } else if (n_fcols == 2) {
    col_args <- fcols_10x_v2
  } else if (n_fcols == 1) {
    col_args <- fcols_haircut
  } else {
    stop("unknown feature file format", call. = FALSE)
  }

  features <- suppressMessages(do.call(
    readr::read_tsv,
    c(
      file = filenames$features,
      col_args
    )
  ))

  bcs <- readLines(filenames$barcodes)

  mat <- Matrix::readMM(filenames$matrix)

  if (strip_10x_suffix) {
    bcs <- gsub("-[0-9]+$", "", bcs)
  }

  if (!is.null(cell_prefix)) {
    bcs <- paste(cell_prefix, bcs, sep = "_")
  }

  if (use_gene_symbols & ("gene_symbol" %in% colnames(features))) {
    features_ids <- features[["gene_symbol"]]
  } else {
    # default to first column
    features_ids <- features[[1]]
  }

  features_ids <- make.unique(features_ids)
  rownames(mat) <- features_ids
  colnames(mat) <- bcs

  mat
}

#' Return number of cols from first line of a file
#' @noRd
count_cols <- function(file,
                       tokenizer_fun = tokenizer_tsv()) {
  readr::count_fields(file,
    tokenizer = tokenizer_fun, n_max = 1
  )
}


# Column definitions for features.tsv files -----------------------------

#' @noRd
fcols_10x_v3 <- list(
  col_types = "ccc",
  col_names = c("gene_id", "gene_symbol", "type")
)
#' @noRd
fcols_10x_v2 <- list(
  col_types = "cc",
  col_names = c("gene_id", "gene_symbol")
)
#' @noRd
fcols_haircut <- list(
  col_types = "c",
  col_names = "feature"
)

#' Convert umitools flat format tsv to sparseMatrix .mtx format
#'
#' @param count_file path to umitools output file
#' @param output_path path for output files. matrix.mtx.gz, barcodes.tsv.gz and
#'   features.tsv.gz will be generated at the supplied path, or by default
#'   created in same directory as the `count_file`.
#' @param ... additional arguments to pass to [`readr::read_tsv()`]
#'
#' @importFrom R.utils gzip
#' @importFrom readr read_tsv
#' @importFrom Matrix readMM
#' @importFrom fs path_dir
#'
#' @export
umitools_to_mtx <- function(count_file,
                            output_path = NULL,
                            ...) {
  dat <- readr::read_tsv(count_file, ...)

  barcodes <- unique(dat$cell)
  genes <- unique(dat$gene)

  dat$gene_idx <- match(dat$gene, genes)
  dat$cell_idx <- match(dat$cell, barcodes)

  mat <- Matrix::sparseMatrix(
    i = dat$gene_idx,
    j = dat$cell_idx,
    x = dat$count
  )

  if (is.null(output_path)) {
    output_path <- fs::path_dir(count_file)
  }

  if (!dir.exists(output_path)) {
    dir.create(output_path, recursive = TRUE)
  }

  Matrix::writeMM(mat, file.path(output_path, "matrix.mtx"))
  R.utils::gzip(file.path(output_path, "matrix.mtx"),
    overwrite = TRUE, remove = TRUE
  )

  readr::write_lines(genes, file.path(output_path, "features.tsv.gz"))
  readr::write_lines(barcodes, file.path(output_path, "barcodes.tsv.gz"))
}


#' Write a sparseMatrix to disk
#'
#' @param mat sparseMatrix
#' @param output_path path for output files. matrix.mtx.gz, barcodes.tsv.gz and
#'   features.tsv.gz will be generated at the supplied path
#'
#' @export
write_matrix <- function(mat, output_path) {

  output_files <- c(
    matrix = "matrix.mtx.gz",
    barcodes = "barcodes.tsv.gz",
    features = "features.tsv.gz"
  )
  output_files <- map(output_files,
                      ~path_join(c(output_path, .x)))

  if (any(file_exists(unlist(output_files)))) {
    warning(paste(
      unlist(output_files),
      "already exist(s), matrices will be overwritten"
    ),
    call. = FALSE
    )
  } else {
    dir.create(output_path, showWarnings = FALSE)
  }


  readr::write_lines(
    colnames(mat),
    output_files$barcodes
  )

  readr::write_lines(
    rownames(mat),
    output_files$features
  )

  uncomp_matrix_fn <- gsub(".gz", "", output_files$matrix)

  Matrix::writeMM(mat, uncomp_matrix_fn)

  R.utils::gzip(uncomp_matrix_fn,
                overwrite = TRUE, remove = TRUE)
}

#' Filter and write a sparseMatrix keeping only specified barcodes
#'
#' sparseMatrix will be read from disk, filtered, and then written to disk.
#'
#' @param matrix_path  path to directory with mtx matrix
#' @param barcodes_path path to barcodes.tsv.gz file generated by 10x count
#' @param output_path path for output files. matrix.mtx.gz, barcodes.tsv.gz and
#'   features.tsv.gz will be generated at the supplied path
#' @param strip_10x_suffix if TRUE, remove 10x suffix from barcodes read from
#'   `barcodes_path` (default = TRUE)
#'
#' @export
filter_matrix <- function(matrix_path,
                          barcodes_path,
                          output_path,
                          strip_10x_suffix = TRUE) {

  mat <- read_matrix(matrix_path,
                     strip_10x_suffix = strip_10x_suffix,
                     use_gene_symbols = TRUE)
  bcs <- readr::read_lines(barcodes_path)

  if(strip_10x_suffix){
    bcs <- gsub("-[0-9]+$", "", bcs)
  }

  shared_bcs <- intersect(colnames(mat), bcs)
  message(glue("there are {n_cells} barcodes remaining in the filtered data",
               n_cells = length(shared_bcs)))

  mat <- mat[, shared_bcs]

  write_matrix(mat, output_path)
}


#' Adds useful labels to colData.
#' Can be used to add cell tyles to cluster numbers
#'
#' @param fsce An object of class [`FunctionalSingleCellExperiment`]
#' @param labels dataframe of new labels. Must contain at least one column of matching variables (e.g. cell_id or k_cluster)
#' @param by column name labels to match in colData of `expt`. If NULL, will match by all matching column names
#' @param expt Data to use for match labels
#'   (default is `rnaseq`). Must be present in `names(fsce)`.
#'
#' @return fsce with all `labels` in `expt` colData.
#'
#' @examples
#' # Add cell_type labels to PBMC data
#'
#' labels <- data.frame(k_cluster = as.factor(c("1", "2", "3", "4", "5", "6")),
#' label = c("MC", "NK", "NK+T", "MC", "MK", "CD4/8 T"))
#'
#' fsce <- add_label(fsce_small, labels)
#'
#' SingleCellExperiment::colData(fsce[["rnaseq"]])
#'
#' @export

add_label <- function(fsce,
                      labels,
                      by = NULL,
                      expt = 'rnaseq') {

  df <- as.data.frame(colData(fsce[[expt]])) %>%
    left_join(labels, by = by)

  colData(fsce[[expt]]) <- DataFrame(df, row.names = df$cell_id)

  fsce
}

