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
#' @importFrom dplyr left_join
#' @importFrom tibble rownames_to_column column_to_rownames
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

#' Create a haircut object as a MultiAssayExperiment
#'
#' @param rna_data UMI count matrix
#' @param haircut_data Functional count matrix
#' @param id_from_name Extract sample id from cell name (default = TRUE)
#' If false then the sample_id field will be populated with NA.
#' @param id_delim delimiter to split cell name (default = ".")
#' @param id_fields index(es) of fields to extract from name for determining sample id from name
#' (defaults to c(2))
#' @param adduct_positions optional data.frame with positions of first and second adducts in each hairpin,
#'  must be a three column data.frame with the first column named hairpin containing with entries that match the adducts
#'  from haircut_data (i.e. Uracil from entry Uracil_1), a second column named adduct_position1 with the first adduct position,
#'  and a third column named adduct_position2 with the second adduct position
#'
#' @importFrom stringr fixed str_split
#' @import SingleCellExperiment
#' @import MultiAssayExperiment
#' @export
create_haircut <- function(rna_data,
                           haircut_data,
                           id_from_name = TRUE,
                           id_delim = ".",
                           id_fields = 2,
                           adduct_positions = NULL){

  # keep umis as sparseMatrix
  rna_mat <- as(as.matrix(rna_data), "sparseMatrix")

  # haircut data is not sparse
  hcut_mat <- as.matrix(haircut_data)

  sce_mrna <- SingleCellExperiment(assays = list(
    counts = rna_mat))

  rna_cells <- colnames(rna_mat)
  hcut_cells <- colnames(hcut_mat)

  # get info for adduct id and position
  # for now assume rownames of input hcut_matrix have this structure
  # ADDUCTID_POS
  hairpin_fields <- stringr::str_split(rownames(hcut_mat),
                                       "_",
                                       simplify = T)

  if(ncol(hairpin_fields) != 2){
    stop("expecting hairpins to be named as ADDUCTNAME_POS, e.g. Uracil_1")
  }

  hairpin_id <- hairpin_fields[, 1]
  hairpin_pos <- hairpin_fields[, 2]

  adduct_data <- data.frame(
    row.names = rownames(hcut_mat),
    hairpin = hairpin_id,
    position = hairpin_pos,
    stringsAsFactors = FALSE
  )

  if (!is.null(adduct_positions)){
    adduct_data <- tibble::rownames_to_column(adduct_data, "hairpin_pos")
    adduct_data <- dplyr::left_join(adduct_data, adduct_positions, by = c("hairpin"))
    adduct_data <- tibble::column_to_rownames(adduct_data, "hairpin_pos")

    if(nrow(adduct_data) != nrow(hcut_mat)){
      stop("unable to add adduct_data to rowData for hairpin object")
    }

  }

  sce_hcut <- SingleCellExperiment(assays = list(
    counts = hcut_mat),
    rowData = adduct_data)

  expr_list <- ExperimentList(list(rna = sce_mrna,
                                   hcut = sce_hcut))

  ## make map between cell ids and experimental assay
  nshared_cells <- intersect(rna_cells, hcut_cells)

  if(length(nshared_cells) == 0){
    warning("No cell ids were shared between the haircut data and rna data")
  }

  rna_cell_map <- data.frame(
    primary = rna_cells,
    colname = rna_cells,
    stringsAsFactors = FALSE
  )

  hcut_cell_map <- data.frame(
    primary = hcut_cells,
    colname = hcut_cells,
    stringsAsFactors = FALSE
  )

  maplist <- list(rna = rna_cell_map,
                  hcut = hcut_cell_map)

  sampMap <- listToMap(maplist)

  ## build colData
  cell_ids <- unique(rna_cells,
                     hcut_cells)

  if (id_from_name){
    sample_ids <- stringr::str_split(cell_ids,
                                     stringr::fixed(id_delim),
                                     simplify = TRUE)

    if(!all(id_fields %in% 1:ncol(sample_ids))){
      warning(paste0("supplied id_fields not found in cell_ids\n",
                     "replacing all sample_ids with NA"))
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

  hce <- MultiAssayExperiment::MultiAssayExperiment(experiments = expr_list,
                                                    colData = colDat,
                                                    sampleMap = sampMap)

  ## set colData for each assay
  colData(hce[["rna"]]) <- colData(hce)[rna_cells, , drop = F]
  colData(hce[["hcut"]]) <- colData(hce)[hcut_cells, , drop = F]

  hce
}
