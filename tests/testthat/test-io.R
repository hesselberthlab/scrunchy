context("test-io")

# Reading -----------------------------------------------------------

mrna_matrix_path <- scrunchy_data("mrna")
hc_matrix_path <- scrunchy_data("haircut")

test_that("sparse matrix import works", {

  mrna_mat <- read_matrix(mrna_matrix_path)
  expect_true(is(mrna_mat, "sparseMatrix"))
  expect_equal(dim(mrna_mat), c(9462, 250))

  hc_mat <- read_matrix(hc_matrix_path)
  expect_true(is(hc_mat, "sparseMatrix"))
  expect_equal(dim(hc_mat), c(426, 250))

})

test_that("10x_suffixes are stripped correctly",{
  mrna_mat <- read_matrix(mrna_matrix_path)
  expect_false(any(grepl("-[0-9]+$", colnames(mrna_mat))))

  mrna_mat <- read_matrix(mrna_matrix_path,
                          strip_10x_suffix = FALSE)
  expect_true(all(grepl("-[0-9]+$", colnames(mrna_mat))))
})

test_that("string prefixing colnames works",{

  mrna_mat <- read_matrix(mrna_matrix_path,
                          cell_prefix = "foo")
  expect_true(all(grepl("^foo_", colnames(mrna_mat))))
})


test_that("gene symbols or ids can be imported",{

  mrna_mat <- read_matrix(mrna_matrix_path,
                          use_gene_symbols = FALSE)
  expect_true(all(grepl("^ENSG", rownames(mrna_mat))))
})

test_that("split_matrix = TRUE splits v3 matrices by feature type",{

  mrna_mat <- read_matrix(mrna_matrix_path,
                          split_matrix = TRUE)
  expect_equal(names(mrna_mat), "Gene Expression")
  expect_equal(typeof(mrna_mat), "list")
  expect_equal(length(mrna_mat), 1L)
})

# Writing -----------------------------------------------------------

mrna_mat <- read_matrix(mrna_matrix_path)
outfiles <- c(
  "matrix.mtx.gz",
  "features.tsv.gz",
  "barcodes.tsv.gz"
)

test_that("writing sparseMatrix works",{
  tmp <- tempdir()
  write_matrix(mrna_mat, tmp)
  outfile_paths <- purrr::map_chr(outfiles,
                                  ~fs::path_join(c(tmp, .x)))
  on.exit(unlink(outfile_paths))

  expect_true(all(fs::file_exists(outfile_paths)))

  mat <- read_matrix(tmp)
  expect_true(is(mat, "sparseMatrix"))

})


test_that("filtering sparseMatrix works",{

  # select 10 barcodes to keep
  ten_bcs <- tibble::tibble(bcs = colnames(mrna_mat)[1:10])
  tmp_bcs <- tempfile()
  on.exit(unlink(tmp_bcs))

  readr::write_tsv(ten_bcs, tmp_bcs, col_names = FALSE)

  tmp_path <- tempdir()
  outfile_paths <- purrr::map_chr(outfiles,
                                  ~fs::path_join(c(tmp_path, .x)))
  on.exit(unlink(outfile_paths))

  expect_message(filter_matrix(mrna_matrix_path,
                 tmp_bcs,
                 tmp_path),
                "there are 10 barcodes remaining in the filtered data")

  expect_true(all(fs::file_exists(outfile_paths)))
  mat <- read_matrix(tmp_path)

  expect_true(is(mat, "sparseMatrix"))
  expect_equal(ncol(mat), 10L)
})

