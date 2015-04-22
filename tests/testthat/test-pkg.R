context("metab model init")

test_that("can create basic metab object", {
  expect_is(metab_model(),'metab_model')

})