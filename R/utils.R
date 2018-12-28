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

  filenames <- list(matrix = matrix_fn,
                    features = features_fn,
                    barcodes = barcodes_fn)

  fns <- dir(path, full.names = TRUE)

  filenames <- map(filenames, ~fs::path_join(c(path, .x)))

  if (!all(filenames %in% fns)) {
    # check for gzipped equivalents
    gzipped_fns <- paste0(filenames, ".gz")
    is_gzipped <- gzipped_fns %in% fns

    # rename if found
    if(any(is_gzipped)) {
      filenames[is_gzipped] <- gzipped_fns[is_gzipped]
    } else {
      stop(paste0("missing required files: ",
                  unlist(filenames[!filenames %in% fns]),
                  "\n"))
    }
  }

  # assign column names based on feature file type
  n_fcols <- count_cols(filenames$features)
  if (n_fcols == 3) {
    col_args = fcols_10x_v3
  } else if (n_fcols == 2) {
    col_args = fcols_10x_v2
  } else if (n_fcols == 1) {
    col_args = fcols_haircut
  } else {
    stop("unknown feature file format", call. = FALSE)
  }

  features <- suppressMessages(do.call(readr::read_tsv,
                                       c(file = filenames$features,
                                         col_args)))

  bcs <- readLines(filenames$barcodes)

  mat <- Matrix::readMM(filenames$matrix)

  if (strip_10x_suffix) {
    bcs <- stringr::str_remove(bcs, "-[0-9]+$")
  }

  if (!is.null(cell_prefix)) {
    bcs <- stringr::str_c(cell_prefix, "_", bcs)
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
                       tokenizer_fun = tokenizer_tsv()){
  readr::count_fields(file,
                      tokenizer = tokenizer_fun, n_max = 1)
}


# Column definitions for features.tsv files -----------------------------

#' @noRd
fcols_10x_v3 <- list(col_types = 'ccc',
                     col_names = c("gene_id", "gene_symbol", "type"))
#' @noRd
fcols_10x_v2 <- list(col_types = "cc",
                     col_names = c("gene_id", "gene_symbol"))
#' @noRd
fcols_haircut <- list(col_types = "c",
                      col_names = "feature")

#' Convert umitools flat format tsv to sparseMatrix .mtx format
#'
#' @param count_file path to umitools output file
#' @param output_path path for output files. matrix.mtx.gz, barcodes.tsv.gz and features.tsv.gz
#' will be generated at the supplied path, or by default created in same directory as
#' the `count_file`.
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

  mat <- Matrix::sparseMatrix(i = dat$gene_idx,
                              j = dat$cell_idx,
                              x = dat$count)

  if(is.null(output_path)) {
    output_path <- fs::path_dir(count_file)
  }

  if(!dir.exists(output_path)) {
    dir.create(output_path, recursive = TRUE)
  }

  Matrix::writeMM(mat, file.path(output_path, "matrix.mtx"))
  R.utils::gzip(file.path(output_path, "matrix.mtx"),
                overwrite = TRUE, remove = TRUE)

  readr::write_lines(genes, file.path(output_path, "features.tsv.gz"))
  readr::write_lines(barcodes, file.path(output_path, "barcodes.tsv.gz"))
}

#' Create a functional cell experiment (fce) object as a MultiAssayExperiment
#'
#' @param rna_data UMI count matrix
#' @param functional_data Functional count matrix
#' @param id_from_name Extract sample id from cell name (default = TRUE) If
#'   false then the sample_id field will be populated with NA.
#' @param id_delim delimiter to split cell name (default = ".")
#' @param id_fields index(es) of fields to extract from name for determining
#'   sample id from name (defaults to c(2))
#' @param adduct_positions optional data.frame with positions of first and
#'   second adducts in each hairpin, must be a three column data.frame with the
#'   first column named hairpin containing with entries that match the adducts
#'   from haircut_data (i.e. Uracil from entry Uracil_1), a second column named
#'   adduct_position1 with the first adduct position, and a third column named
#'   adduct_position2 with the second adduct position
#'
#' @return fce object of class MultiAssayExperiment containing
#'   SingleCellExperiments. mRNA data is stored in slot "sce", and functional
#'   data is stored in slot "fsce"
#'
#' @export
create_fce <- function(rna_data,
                       functional_data,
                       id_from_name = TRUE,
                       id_delim = ".",
                       id_fields = 2,
                       adduct_positions = NULL) {

  # keep umis as sparseMatrix
  rna_mat <- as(as.matrix(rna_data), "sparseMatrix")

  # haircut data is not sparse
  f_mat <- as.matrix(functional_data)

  sce <- SingleCellExperiment(assays = list(
    counts = rna_mat
  ))

  rna_cells <- colnames(rna_mat)
  f_cells <- colnames(f_mat)

  # get info for adduct id and position
  # for now assume rownames of input f_matrix have this structure
  # ADDUCTID_POS
  hairpin_fields <- stringr::str_split(rownames(f_mat),
    "_",
    simplify = T
  )

  if (ncol(hairpin_fields) != 2) {
    stop("expecting hairpins to be named as ADDUCTNAME_POS, e.g. Uracil_1")
  }

  hairpin_id <- hairpin_fields[, 1]
  hairpin_pos <- hairpin_fields[, 2]

  adduct_data <- data.frame(
    row.names = rownames(f_mat),
    hairpin = hairpin_id,
    position = hairpin_pos,
    stringsAsFactors = FALSE
  )

  if (!is.null(adduct_positions)) {
    adduct_data <- tibble::rownames_to_column(adduct_data, "hairpin_pos")
    adduct_data <- dplyr::left_join(adduct_data, adduct_positions, by = c("hairpin"))
    adduct_data <- tibble::column_to_rownames(adduct_data, "hairpin_pos")

    if (nrow(adduct_data) != nrow(f_mat)) {
      stop("unable to add adduct_data to rowData for hairpin object")
    }
  }

  fsce <- SingleCellExperiment(
    assays = list(
      counts = f_mat
    ),
    rowData = adduct_data
  )

  expr_list <- ExperimentList(list(
    sce = sce,
    fsce = fsce
  ))

  ## make map between cell ids and experimental assay
  nshared_cells <- intersect(rna_cells, f_cells)

  if (length(nshared_cells) == 0) {
    warning("No cell ids were shared between the haircut data and rna data")
  }

  rna_cell_map <- data.frame(
    primary = rna_cells,
    colname = rna_cells,
    stringsAsFactors = FALSE
  )

  f_cell_map <- data.frame(
    primary = f_cells,
    colname = f_cells,
    stringsAsFactors = FALSE
  )

  maplist <- list(
    sce = rna_cell_map,
    fsce = f_cell_map
  )

  sampMap <- listToMap(maplist)

  ## build colData
  cell_ids <- unique(
    rna_cells,
    f_cells
  )

  if (id_from_name) {
    sample_ids <- stringr::str_split(cell_ids,
      stringr::fixed(id_delim),
      simplify = TRUE
    )

    if (!all(id_fields %in% 1:ncol(sample_ids))) {
      warning(paste0(
        "supplied id_fields not found in cell_ids\n",
        "replacing all sample_ids with NA"
      ))
      sample_ids <- rep(NA, length(cell_ids))
    } else {
      sample_ids <- sample_ids[, id_fields]
    }
  } else {
    sample_ids <- rep(NA, length(cell_ids))
  }

  colDat <- data.frame(
    row.names = cell_ids,
    sample_id = sample_ids
  )

  fce <- MultiAssayExperiment::MultiAssayExperiment(
    experiments = expr_list,
    colData = colDat,
    sampleMap = sampMap
  )

  ## set colData for each assay
  colData(fce[["sce"]]) <- colData(fce)[rna_cells, , drop = F]
  colData(fce[["fsce"]]) <- colData(fce)[f_cells, , drop = F]

  fce
}

#' Provide path to scrunchy internal data
#'
#' @param x file name
#'
#' @export
scrunchy_data <- function(x) {
  system.file("extdata", x, package = "scrunchy", mustWork = TRUE)
}

#' Load CSV format files
#'
#' @param x file name
#'
#' @export
load_csv <- function(x) {
  read.csv(x, sep = ",", header = TRUE, row.names = 1)
}
