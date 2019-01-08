context("test-stats")

test_that("aov results can be tidied", {
  x <- fsce_tidy[c("k_cluster", "Uracil_45", "riboG_44")]
  x$k_cluster <- as.factor(x$k_cluster)

  res <- stat_anova_grouped(x, k_cluster, tidy = TRUE)
  expect_equal(dim(res), c(4, 7))
})

test_that("aov post-hoc results can be tidied", {
  x <- fsce_tidy[c("k_cluster", "Uracil_45", "riboG_44")]
  x$k_cluster <- as.factor(x$k_cluster)

  ares <- stat_anova_grouped(x, k_cluster)
  res <- stat_anova_tukey(res, tidy = TRUE)

  expect_equal(dim(res), c(30, 4))
})

test_that("factors can be completed", {
  x <- fsce_tidy[c("k_cluster", "Uracil_45", "riboG_44")]
  res <- stat_activity_grouped(x, group = k_cluster)
  expect_equal(nrow(res), 30)

  res <- stat_activity_grouped(x, group = k_cluster, complete = TRUE)
  expect_equal(nrow(res), 72)
})

test_that("cross_groups returns unique combinations", {
  x <- select(fsce_tidy, k_cluster, Uracil_45, riboG_44)
  res <- stat_activity_grouped(x, group = k_cluster)
  expect_equal(dim(res), c(30, 4))
})
