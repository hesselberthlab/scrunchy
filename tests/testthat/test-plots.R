context("test-plots")

test_that("incorrect labels throw an error", {
  expect_error(
    plot_dims(fsce_tidy, UMAP1, UMAP2, k_cluster, labels = LETTERS[1:10]),
    "must match factors"
  )
})
