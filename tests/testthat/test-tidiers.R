context("test-tidiers")

# haircut subset
fsce_hc <- fsce_small[c("Uracil_45"), , ]
# rnaseq subset
fsce_rs <- fsce_small[c("TP53"), , ]

# tidy_logcounts ----------------------------------------------------

test_that("tidy_logcounts results have expected shapes", {
  res_hc <- tidy_logcounts(fsce_hc)
  res_rs <- tidy_logcounts(fsce_rs)

  expect_equal(dim(res_hc), c(250, 3))
  expect_equal(dim(res_rs), c(250, 3))
})

# tidy_counts ----------------------------------------------------

test_that("tidy_counts results have expected shapes", {
  res_hc <- tidy_counts(fsce_hc)
  res_rs <- tidy_counts(fsce_rs)

  expect_equal(dim(res_hc), c(250, 3))
  expect_equal(dim(res_rs), c(250, 3))
})

# tidy_dims ----------------------------------------------------

test_that("tidy_dims results have expected shapes", {
  res_hc <- tidy_dims(fsce_hc)
  res_rs <- tidy_dims(fsce_rs)

  expect_equal(dim(res_hc), c(0, 2))
  expect_equal(dim(res_rs), c(750, 9))
})

# test_that("tidy_dims can filter", {
#   expect_equal(
#     dim(tidy_dims(fsce, dims = c(1))),
#     c(750, 9)
#   )
# })
