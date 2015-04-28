library(unitted)

context("calc_DO_at_sat works as expected w/ unitted")
test_that("proper units checked and returned for calc_DO_at_sat", {
  expect_equal(calc_DO_at_sat(temp=u(21,"degC"), press=u(1.1,"mb"), sal=u(0,"PSU")), u(0.0087838, "L^-1 mg"), tol=0.0001, info="with units")
  expect_equal(calc_DO_at_sat(temp=21, press=1.1, sal=0), 0.0087838, tol=0.0001, info="with no units")
  expect_error(calc_DO_at_sat(temp=u(77,"degF"), press=u(1.1,"atm"), sal=0), "Unexpected units", info="with wrong units")  
})


context("calc_DO_deficit works as expected w/ unitted")
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
  expect_warning(expect_error(calc_DO_deficit(DO.obs, temp.water, pressure.air, v(salinity.water))))
  expect_error(calc_DO_deficit(u(DO.obs,'dogs'), temp.water, pressure.air, salinity.water))
  expect_error(calc_DO_deficit(DO.obs, u(temp.water,'dogs'), pressure.air, salinity.water))
  expect_equal(is.numeric(calc_DO_deficit(DO.obs=c(7,7.5,7), temp=21, press=1.1, sal=0)), TRUE, info="with no units")
  expect_error(calc_DO_deficit(DO.obs=u(7,"mg L^-1"), temp=4, press=1.1, sal=0), info="notice mismatched units in subtraction")  
  expect_error(calc_DO_deficit(DO.obs=u(7,"g m^-3"), temp=u(4,"degC"), press=u(1100,"mb"), sal=u(0,"PSU")), info="notice mismatched units in subtraction")  
  expect_equal(get_units(calc_DO_deficit(DO.obs=u(7,"mg L^-1"), temp=u(4,"degC"), press=u(1100,"mb"), sal=u(0,"PSU"))), "mg L^-1", info="notice mismatched units in subtraction")  
})

test_that("unit fail with partial specs", {
  expect_error(calc_DO_deficit(v(DO.obs), u(temp.water,'dogs'), pressure.air, salinity.water))
})