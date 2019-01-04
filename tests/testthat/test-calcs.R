context("test-calcs")

test_that("calc functions are reproducible with a seed", {
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
  expect_equal(
    colData(calc_kmeans(fsce_small, k = 6, seed = seed)[["rnaseq"]]),
    colData(calc_kmeans(fsce_small, k = 6, seed = seed)[["rnaseq"]])
  )
})
