context("calc_DO_deficit works as expected w/ unitted")
library(unitted)

test_that("proper units checked and returned for DO.deficit", {
  
  DO.obs <<- u(c(7,7.5,7),'mg L^-1')
  temp.water <<- u(c(25,24.5,18.9), 'degC')
  pressure.air <<- u(c(900,903,910), 'mb')
  salinity.water <<- u(2.43, 'PSU')
  DO.deficit <- calc_DO_deficit(DO.obs, temp.water, pressure.air, salinity.water)
  
  expect_is(DO.deficit,'unitted')
  expect_equal(get_units(DO.deficit),'mg L^-1')
})
test_that("proper unit failures for DO.deficit", {
  expect_error(calc_DO_deficit(DO.obs, temp.water, pressure.air, v(salinity.water)))
  expect_error(calc_DO_deficit(u(DO.obs,'dogs'), temp.water, pressure.air, salinity.water))
  expect_error(calc_DO_deficit(DO.obs, u(temp.water,'dogs'), pressure.air, salinity.water))
})

test_that("unit fallback with partial specs", {
  DO.deficit <- calc_DO_deficit(v(DO.obs), u(temp.water,'dogs'), pressure.air, salinity.water)
  expect_is(DO.deficit,'numeric')
})