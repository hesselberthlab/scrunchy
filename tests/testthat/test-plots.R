context("test-plots")

test_that("plot_heatmap can accept sparseMatrix input", {
  i <- c(1,3:8); j <- c(2,9,6:10); x <- 7 * (1:7)
  m <- Matrix::sparseMatrix(i, j, x = x)

  p <- plot_heatmap(m)
  expect_is(p, "Heatmap")
})

test_that("incorrect labels throw an error", {
  expect_error(
    plot_dims(fsce_tidy, UMAP1, UMAP2, k_cluster, labels = LETTERS[1:10]),
    "must match factors"
  )
})
