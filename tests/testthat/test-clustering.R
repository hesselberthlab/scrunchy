context("test-clustering")

test_that("clusters are calculated", {

  expt <- "rnaseq"
  method <- "PCA"

  ## k-means
  k <- 6
  colname <- "k_cluster"
  res <- cluster_kmeans(
    fsce_small,
    expt = expt, method = method, k = k
  )
  expect_length(colData(res[[expt]])[[colname]], 250)

  ## leiden
  colname <- "leiden_cluster"
  res <- cluster_leiden(
    fsce_small,
    expt = expt, method = method
  )
  expect_length(colData(res[[expt]])[[colname]], 250)

})

test_that("functions are reproducible with a seed", {
  seed <- 47681
  expect_equal(
    colData(cluster_kmeans(fsce_small, k = 6, seed = seed)[["rnaseq"]]),
    colData(cluster_kmeans(fsce_small, k = 6, seed = seed)[["rnaseq"]])
  )
})

test_that("inputs are checked", {

  fsce_empty <- FunctionalSingleCellExperiment(
    experiments = list(sce_empty = SingleCellExperiment())
  )

  ## cluster_kmeans
  expect_error(
    cluster_kmeans(fsce_empty),
    "not found in fsce"
  )
  expect_error(
    cluster_kmeans(fsce_empty, expt = "sce_empty"),
    "method `PCA` not found in expt"
  )

  ## cluster_leiden
  expect_error(
    cluster_leiden(fsce_empty),
    "not found in fsce"
  )
  expect_error(
    cluster_leiden(fsce_empty, expt = "sce_empty"),
    "method `PCA` not found in expt"
  )
})
