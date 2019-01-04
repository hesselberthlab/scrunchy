context("test-stats")

test_that("cross_groups returns unique combinations", {
  x <- select(fsce_tidy, k_cluster, Uracil_45, riboG_44)
  res <- calc_group_stats(x, group = k_cluster)
  expect_equal(dim(res), c(30, 4))
})
