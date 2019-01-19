context("test-classes")

test_that("fsce class can be created", {
  fsce <- create_fsce(
    list(
      rnaseq = create_sce_rnaseq(scrunchy_data("mrna")),
      haircut = create_sce_haircut(scrunchy_data("haircut"))
    )
  )

  ## assays
  expect_equal(names(assays(fsce)), c("rnaseq", "haircut"))

  ## experiments
  expect_length(experiments(fsce), 2)

  e1 <- experiments(fsce_small)[[1]]
  e2 <- experiments(fsce_small)[[2]]

  expect_equal(dim(e1), c(9462, 250))
  expect_equal(dim(e2), c(426, 250))
})
