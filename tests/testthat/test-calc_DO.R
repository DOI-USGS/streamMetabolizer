context("calc_DO")

library(unitted)

DO.obs <<- u(c(7,7.5,7),'mgO2 L^-1')
temp.water <<- u(c(25,24.5,18.9), 'degC')
pressure.air <<- u(c(900,903,910), 'mb')
salinity.water <<- u(2.43, 'PSU')

test_that("proper units checked and returned for calc_DO_sat", {
  expect_equal(calc_DO_sat(temp=u(21,"degC"), press=u(1013.25,"mb"), sal=u(0,"PSU")), u(8.914559, "L^-1 mgO2"), tol=0.0001, info="with units")
  expect_equal(calc_DO_sat(temp=21, press=1013.25, sal=0), 8.914559, tol=0.0001, info="with no units")
  expect_error(calc_DO_sat(temp=u(77,"degF"), press=u(1.1,"atm"), sal=0), "Unexpected units", info="with wrong units")  
})

test_that("proper unit failures for DO.deficit", {
  expect_error(calc_DO_sat(u(temp.water,'dogs'), pressure.air, salinity.water), 'Unexpected units')
  expect_error(expect_warning(calc_DO_sat(v(temp.water), pressure.air, salinity.water), 'not unitted'), 'Unexpected units')
  expect_equal(is.numeric(calc_DO_sat(temp=21, press=1.1, sal=0)), TRUE)
  expect_equal(get_units(calc_DO_sat(temp=u(4,"degC"), press=u(1100,"mb"), sal=u(0,"PSU"))), "mgO2 L^-1")  
})
