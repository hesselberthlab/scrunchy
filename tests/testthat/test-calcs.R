context("test-calcs")

test_that("inputs are checked", {

  fsce_empty <- FunctionalSingleCellExperiment(
    experiments = list(sce_empty = SingleCellExperiment())
  )

  ## calc_var_features
  expect_error(
    calc_var_features(fsce_empty),
    "not found in fsce"
  )
  expect_error(
    calc_var_features(fsce_empty, expt = "sce_empty"),
    "not found for expt"
  )

  ## calc_pca
  expect_error(
    calc_pca(fsce_empty),
    "not found in fsce"
  )
  expect_error(
    calc_pca(fsce_empty, expt = "sce_empty"),
    "`logcounts` not found in fsce"
  )

  ## calc_umap
  expect_error(
    calc_umap(fsce_empty),
    "not found in fsce"
  )
  expect_error(
    calc_umap(fsce_empty, expt = "sce_empty"),
    "method `PCA` not found in fsce"
  )

  ## calc_tsne
  expect_error(
    calc_tsne(fsce_empty),
    "not found in fsce"
  )
  expect_error(
    calc_tsne(fsce_empty, expt = "sce_empty"),
    "method `PCA` not found in fsce"
  )

})

test_that("functions are reproducible with a seed", {
  var_genes <- calc_var_features(fsce_small, "rnaseq", n = 50)

  seed <- 47681

  expect_equal(
    reducedDim(calc_pca(fsce_small, genes = var_genes, seed = seed)[["rnaseq"]], "PCA"),
    reducedDim(calc_pca(fsce_small, genes = var_genes, seed = seed)[["rnaseq"]], "PCA")
  )
  expect_equal(
    reducedDim(calc_umap(fsce_small, n_dims = 2, seed = seed)[["rnaseq"]], "UMAP"),
    reducedDim(calc_umap(fsce_small, n_dims = 2, seed = seed)[["rnaseq"]], "UMAP")
  )
  expect_equal(
    reducedDim(calc_tsne(fsce_small, n_dims = 2, seed = seed)[["rnaseq"]], "TSNE"),
    reducedDim(calc_tsne(fsce_small, n_dims = 2, seed = seed)[["rnaseq"]], "TSNE")
  )
})
