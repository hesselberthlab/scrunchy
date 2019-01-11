context("test-stats")

# ANOVA -------------------------------------------------------------

test_that("anova results are ok", {
  x <- fsce_tidy[c("k_cluster", "Uracil_45")]
  x$k_cluster <- as.factor(x$k_cluster)

  ## test aov results
  res <- stat_anova_grouped(x, k_cluster)
  expect_length(res, 1)
  expect_true(all(class(res[[1]]) == c("aov", "lm")))

  ## test tidied results
  res_tidy <- tidy_stats_grouped(res)
  expect_equal(dim(res_tidy), c(2,7))
  expect_is(res_tidy, "tbl_df")
})

test_that("anova post-hoc results are ok", {
  x <- fsce_tidy[c("k_cluster", "Uracil_45")]
  x$k_cluster <- as.factor(x$k_cluster)

  res_anova <- stat_anova_grouped(x, group = k_cluster)

  ## test tukey result
  res <- stat_anova_tukey(res_anova, group = k_cluster)
  expect_length(res, 1)
  expect_is(res[[1]], "glht")

  ## test tidied results
  res_tidy <- tidy_stats_grouped(res)
  expect_equal(dim(res_tidy), c(15, 4))
  expect_is(res_tidy, "tbl_df")
})

# Activities -------------------------------------------------------------

test_that("fold-changes are calculated", {
  x <- fsce_tidy[c("k_cluster", "Uracil_45")]
  res <- stat_activity_grouped(x, group = k_cluster)

  expect_true("ratio" %in% names(res))
})

test_that("qvalues are calculated", {
  x <- fsce_tidy[c("k_cluster", "Uracil_45", "riboG_44")]
  res <- stat_activity_grouped(x, group = k_cluster)

  expect_true("q.value" %in% names(res))
  expect_true(all(res$p.value <= res$q.value))
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
  expect_equal(dim(res), c(30, 6))
})
