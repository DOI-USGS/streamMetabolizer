context("mm_determine_cores")

test_that("falls back to 1 core when detectCores() is non-finite", {
  local_mocked_bindings(detectCores = function() NA_real_, .package = "parallel")
  expect_equal(mm_determine_cores(n_cores=5), 1)

  local_mocked_bindings(detectCores = function() Inf, .package = "parallel")
  expect_equal(mm_determine_cores(n_cores=5), 1)
})

test_that("caps the requested core count at the number detected", {
  local_mocked_bindings(detectCores = function() 4, .package = "parallel")
  expect_equal(mm_determine_cores(n_cores=10), 4)
  expect_equal(mm_determine_cores(n_cores=2), 2)
  expect_equal(mm_determine_cores(n_cores=4), 4)
})

test_that("emits a status message only when verbose=TRUE", {
  local_mocked_bindings(detectCores = function() 4, .package = "parallel")

  expect_message(
    result <- mm_determine_cores(n_cores=10, n_chains=3, verbose=TRUE),
    "requesting 3 chains on 4 of 4 available cores",
    fixed=TRUE)
  expect_equal(result, 4)

  expect_no_message(mm_determine_cores(n_cores=10, n_chains=3, verbose=FALSE))
})

test_that("verbose message matches runstan_bayes()'s original wording exactly", {
  local_mocked_bindings(detectCores = function() 8, .package = "parallel")

  expect_message(
    mm_determine_cores(n_cores=4, n_chains=4, verbose=TRUE),
    "MCMC (Stan): requesting 4 chains on 4 of 8 available cores",
    fixed=TRUE)
})

test_that("omits the chains clause when n_chains is not supplied", {
  local_mocked_bindings(detectCores = function() 4, .package = "parallel")

  expect_message(
    mm_determine_cores(n_cores=10, verbose=TRUE),
    "MCMC (Stan): requesting 4 of 4 available cores",
    fixed=TRUE)
})
