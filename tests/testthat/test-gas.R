context("calc_schmidt works as expected w/ unitted")

test_that("proper units checked and returned for schmidt", {
  temp.water <- 25
  gas <- "O2"
  expect_is(calc_schmidt(temp.water, gas), 'numeric')
  temp.water <- unitted::u(25, 'degC')
  expect_is(calc_schmidt(temp.water, gas), 'unitted')
  # is unitted, but wrong units
  expect_error(calc_schmidt(u(25,'degF'), gas))
  
})

test_that("calc_schmidt fails w/ bad args",{
  expect_error(calc_schmidt(25, 'notagas'))
  expect_warning(calc_schmidt(54, 'CO2'))
  expect_warning(calc_schmidt(-12, 'CO2'))
})