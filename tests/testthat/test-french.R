context("french creek data")

test_that("French Creek data are similar for streamMetabolizer & Bob Hall's code", {
  
  # load both datasets
  fx <- streamMetabolizer:::load_french_creek(attach.units=FALSE)
  fy <- streamMetabolizer:::load_french_creek_std()
  expect_equal(dim(fx), dim(fy))
  expect_equal(names(fx), names(fy))
  expect_equal(sapply(unitted::v(fx), class), sapply(fy, class))
  expect_equal(unitted::v(fx$DO.obs), fy$DO.obs)
  expect_equal(unitted::v(fx$local.time), fy$local.time)
  
  # combine for further comparison
  fxy <- dplyr::full_join(unitted::v(fx), fy, by="local.time")
  expect_equal(nrow(fx), nrow(fxy))
  expect_equal(fxy$local.time.x, fxy$local.time.y)
  
  # check values that should be completely equal
  expect_equal(fxy$DO.obs.x, fxy$DO.obs.y)
  expect_equal(fxy$depth.x, fxy$depth.y)
  expect_equal(fxy$temp.water.x, fxy$temp.water.y)
  
  # check values that should be pretty much equal - DO.sat
  expect_more_than(cor(fxy$DO.sat.x, fxy$DO.sat.y), 0.99)
  expect_less_than(abs(coef(lm(fxy$DO.sat.y ~ fxy$DO.sat.x))['(Intercept)']), 0.1)
  expect_less_than(abs(coef(lm(fxy$DO.sat.y ~ fxy$DO.sat.x))['fxy$DO.sat.x'] - 1), 0.1)
  # library(ggplot2); ggplot(fxy, aes(x=DO.sat.x, y=DO.sat.y)) + geom_abline() + geom_point(alpha=0.5) + theme_bw()
  
  # check values that should be pretty much equal - light
  expect_more_than(cor(fxy$light.x, fxy$light.y), 0.9999)
  expect_less_than(abs(coef(lm(fxy$light.y ~ fxy$light.x))['(Intercept)']), 0.4)
  expect_less_than(abs(coef(lm(fxy$light.y ~ fxy$light.x))['fxy$light.x'] - 1), 0.01)
  # library(ggplot2); ggplot(fxy, aes(x=light.x, y=light.y)) + geom_abline() + geom_point(alpha=0.5) + theme_bw()
  # library(ggplot2); ggplot(fxy, aes(x=local.time)) + geom_line(aes(y=light.x), color='red') + geom_line(aes(y=light.y), color='blue') + theme_bw()
  # library(ggplot2); ggplot(fxy[7100:7350,], aes(x=local.time)) + geom_line(aes(y=light.x), color='red') + geom_line(aes(y=light.y), color='blue') + theme_bw()
  
})