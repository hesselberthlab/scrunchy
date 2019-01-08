context("test-dims")

seed <- 42
fsce <- create_fsce(
  list(
    haircut = create_sce_haircut(scrunchy_data("haircut")),
    rnaseq = create_sce_rnaseq(scrunchy_data("mrna"))
  )
)

test_that("dims can be calculated", {
  res <- calc_pca(fsce, n_pcs = 20, seed = seed)
  expect_equal(
    colnames(reducedDim(res[["rnaseq"]])),
    paste0("PC", 1:20)
  )
  res <- calc_umap(res, seed = seed)
  expect_equal(
    dim(reducedDims(res[["rnaseq"]])$UMAP),
    c(250, 2)
  )
  res <- calc_tsne(res, seed = seed)
  expect_equal(
    dim(reducedDims(res[["rnaseq"]])$TSNE),
    c(250, 2)
  )
})

test_that("input files are checked", {
  expect_error(
    calc_pca(fsce_small, "missing"),
    "expt `missing` not found in fsce"
  )
  expect_error(
    calc_umap(fsce_small, "missing"),
    "expt `missing` not found in fsce"
  )
  expect_error(
    calc_tsne(fsce_small, "missing"),
    "expt `missing` not found in fsce"
  )
})
