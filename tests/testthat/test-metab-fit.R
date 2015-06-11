context("output of metabolism models")

# use a subset of data from Bob
library(plyr); library(dplyr)
library(unitted)
french <- streamMetabolizer:::load_french_creek()
#french$date.time <- convert_GMT_to_solartime(french$date.time, longitude=-106.48059) # seems to already be local time, close enough to solar
french$DO.sat <- calc_DO_at_sat(temp.water=french$temp.water, pressure.air=u(1000, "mb"))
french$light <- convert_SW_to_PAR(calc_solar_insolation(date.time=v(french$date.time), latitude=41.22668, attach.units=TRUE))
french <- french[c("date.time","DO.obs","DO.sat","depth","temp.water","light")]


test_that("metabolism models run & produce reasonable output", {
  
  # metab_mle
  mm <- metab_mle(data=v(french))
  expect_is(mm, "metab_mle")
  expect_is(slot(mm, "fit"), "data.frame") # specific to this model
  expect_is(slot(mm, "args"), "list")
  expect_is(slot(mm, "data"), "data.frame")
  expect_is(slot(mm, "pkg_version"), "character")
})

test_that("metabolism predictions (predict_metab, predict_DO) make sense", {
  
  # metab_mle
  mm <- metab_mle(data=v(french))
  metab <- predict_metab(mm)
  expect_equal(metab$GPP + metab$ER, metab$NEP)
  DO_preds <- predict_DO(mm)
  DO_preds_Aug24<- filter(DO_preds, format(date.time, "%Y-%m-%d") == "2012-08-24")
  expect_true(all(abs(DO_preds_Aug24$DO.obs - DO_preds_Aug24$DO.mod) < 0.15), "DO.mod tracks DO.obs with not too much error")
  # library(ggplot2); ggplot(DO_preds, aes(x=date.time)) + geom_line(aes(y=DO.obs), color="blue") + geom_line(aes(y=DO.mod), color="red")
  
})

test_that("metab_models can be saved", {
  
  # save and reload metab_mle
  mm <- metab_mle(data=v(french))
  save(mm, file="test_temp.RData")
  mm_orig <- mm
  load("test_temp.RData")
  file.remove("test_temp.RData")
  expect_equal(mm_orig, mm)
  expect_equal(get_fit(mm_orig), get_fit(mm))
  
})