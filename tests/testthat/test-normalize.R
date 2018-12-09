context("test-normalize")

set.seed(1234)

nr <- 100
nc <- 100

x <- matrix(
  sample(1:1e5, nc * nr, replace = TRUE),
  ncol = nc
)

test_that("mats can be log normlized", {
  expect_is(log_normalize(x), "matrix")
})

test_that("mats can be clr normlized", {
  expect_is(clr_normalize(x), "matrix")
})
