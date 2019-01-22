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

# tidy_coldata ----------------------------------------------------

test_that("tidy_coldata results have expected shapes", {
  res_hc <- tidy_coldata(fsce_hc)
  res_rs <- tidy_coldata(fsce_rs)

  expect_equal(dim(res_hc), c(250, 1))
  expect_equal(dim(res_rs), c(250, 4))
})

# tidy_dims ----------------------------------------------------

test_that("tidy_dims results have expected shapes", {
  res_hc <- tidy_dims(fsce_hc)
  res_rs <- tidy_dims(fsce_rs)

  expect_equal(dim(res_hc), c(0, 1))
  expect_equal(dim(res_rs), c(250, 8))
})

test_that("tidy_dims can filter", {
  expect_equal(
    dim(tidy_dims(fsce_small, dimnames = c("UMAP"))),
    c(250, 4)
  )

  res <- tidy_dims(fsce_small, dimnames = c("PCA"), dims = c(1, 3))
  expect_equal(dim(res), c(250, 4))
  expect_true(all(c("PCA1", "PCA3") %in% names(res)))
})
