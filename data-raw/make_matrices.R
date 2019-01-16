library(scrunchy)
library(here)
library(tidyverse)

# code for producing the small mrna and haircut matrices found in inst/extdata/

# download and read in matrices ------------------------
files <- c(
  "matrix.mtx.gz",
  "barcodes.tsv.gz",
  "features.tsv.gz")

load_matrix_files <- function(base_url, output_dir, files) {

  if(dir.exists(output_dir)) {
    stop(paste(output_dir, "already exists, please specify a temporary directory name"))
  }

  urls <- paste0(base_url, files)

  output_files <- file.path(output_dir, files)
  dir.create(output_dir)

  purrr::walk2(urls, output_files, download.file)

  mat <- read_matrix(output_dir)

  # delete files
  #unlink(output_dir, recursive = T)

  mat
}


write_matrix <- function(mat, output_dir, feature_file = NULL) {
  if(dir.exists(output_dir)) {
    stop(paste(output_dir, "already exists, please specify a new directory name"))
  }

  dir.create(output_dir)

  readr::write_lines(colnames(mat),
                     file.path(output_dir, "barcodes.tsv.gz"))

  # write out feature file with all attributes
  if(!is.null(feature_file)){
    features <- readr::read_tsv(feature_file,
                        col_names = c("id", "symbol", "type"))
    features$symbol <- make.unique(features$symbol)

    features <- left_join(tibble(symbol = rownames(mat)),
                  features, by = "symbol")

    features <- features[, c("id", "symbol", "type")]

    readr::write_tsv(features,
                     file.path(output_dir, "features.tsv.gz"),
                     col_names = FALSE)
  } else {
    readr::write_lines(rownames(mat),
                       file.path(output_dir, "features.tsv.gz"))
  }

  Matrix::writeMM(mat, file.path(output_dir, "matrix.mtx"))

  R.utils::gzip(file.path(output_dir, "matrix.mtx"))
}

mrna_mat <- load_matrix_files("http://amc-sandbox.ucdenver.edu/User33/hcut/mrna/",
                              "tmp_mrna_full",
                              files)

# subsample to 250 cells, remove 0 count genes, and write to disk ------------------------
set.seed(42)

selected_cells <- sample(colnames(mrna_mat), 250)

small_mat <- mrna_mat[, selected_cells]

small_mat <- small_mat[Matrix::rowSums(small_mat) > 0, ]

#add back suffix to make match input
colnames(small_mat) <- paste0(colnames(small_mat), "-1")

write_matrix(small_mat, "tmp_mrna_small", "tmp_mrna_full/features.tsv.gz")

# download, read, and subsample haircut umi matrix ------------------

hc_mat <- load_matrix_files("http://amc-sandbox.ucdenver.edu/User33/hcut/haircut/unfiltered/",
                            "tmp_hc_full",
                             files)

write_matrix(hc_mat[, selected_cells], "tmp_hc_small")

# move files to inst/extdata -----------------------------------------

hc_dir <- file.path(here(), "inst", "extdata", "haircut")
mrna_dir <- file.path(here(), "inst", "extdata", "mrna")

dir.create(hc_dir)
dir.create(mrna_dir)

file.copy(file.path("tmp_mrna_small", files), mrna_dir, overwrite = TRUE)
file.copy(file.path("tmp_hc_small", files), hc_dir, overwrite = TRUE)

unlink("tmp_mrna_small", recursive = TRUE)
unlink("tmp_hc_small", recursive = TRUE)
unlink("tmp_mrna_full", recursive = TRUE)
unlink("tmp_hc_full", recursive = TRUE)
