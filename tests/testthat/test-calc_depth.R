context("calc_depth")

test_that("calc_depth checks & returns values & units as expected", {
  # basic numbers, with & without defaults
  expect_is(calc_depth(Q=7.3), "numeric")
  expect_equal(calc_depth(Q=7.3, c=1, f=0.5), 7.3^0.5)
  expect_equal(calc_depth(Q=rep(100, 3), c=seq(0.4, 0.6, 0.1), f=c(0.25, 0.29, 0.33)), c(1.26, 1.9, 2.74), tol=0.1)
  
  # with units - throws warnings rather than errors in several cases
  library(unitted)
  expect_error(expect_warning(calc_depth(Q=1:10, f=u(0.36)), "not unitted"), "Unexpected units")
  expect_warning(calc_depth(Q=u(1:10, "m^3 s^-1"), c=0.36), "unknown depth units")
  expect_warning(calc_depth(Q=u(1:10, "m^3 s^-1"), f=0.36), "not unitted")
  
  # use units of c, whatever those are, but throw a warning if they're odd
  expect_equal(get_units(calc_depth(Q=u(1:10, "m^3 s^-1"), c=u(40,"cm"))), "cm", "retains units of c")
  expect_warning(calc_depth(Q=u(1:10, "m^3 s^-1"), c=u(40,"kebabs")), "unknown depth units")
  
})
