context("output of metabolism models")

test_that("metabolism models run & produce reasonable output", {
  
  # use a subset of data from Bob
  library(plyr); library(dplyr)
  library(unitted)
  french <- streamMetabolizer:::load_french_broad()
  french$DO.sat <- calc_DO_at_sat(temp.water=french$temp.water, pressure.air=u(1000, "mb"))
  french$light <- convert_SW_to_PAR(calc_solar_insolation(date.time=v(french$date.time), latitude=35.95, attach.units=TRUE))
  french <- french[c("date.time","DO.obs","DO.sat","depth","temp.water","light")]
  
  # metab_simple
  mm <- metab_simple(data=v(french))
  expect_is(mm, "metab_simple")
  expect_is(slot(mm, "fit"), "list") # specific to this model
  expect_equal(names(slot(mm, "fit")), c("minimum","estimate","gradient","code","iterations")) # specific to this model
  expect_is(slot(mm, "args"), "list")
  expect_is(slot(mm, "data"), "data.frame")
  expect_is(slot(mm, "pkg_version"), "character")
})

test_that("metabolism predictions (predict_metab, predict_DO) make sense", {
  
  # use a subset of data from Bob
  library(plyr); library(dplyr)
  library(unitted)
  french <- streamMetabolizer:::load_french_broad()
  french$DO.sat <- calc_DO_at_sat(temp.water=french$temp.water, pressure.air=u(1000, "mb"))
  french$light <- convert_SW_to_PAR(calc_solar_insolation(date.time=v(french$date.time), latitude=35.95, attach.units=TRUE))
  french <- french[c("date.time","DO.obs","DO.sat","depth","temp.water","light")]
  
  # metab_simple
  mm <- metab_simple(data=v(french))
  metab <- predict_metab(mm)
  expect_equal(metab$GPP + metab$ER, metab$NEP)
  DO_preds <- predict_DO(mm)
  expect_true(all(abs(DO_preds$DO.obs - DO_preds$DO.mod) < 0.15), "DO.mod tracks DO.obs with not too much error")
  #library(ggplot2)
  #ggplot(DO_preds, aes(x=date.time)) + geom_line(aes(y=DO.obs), color="blue") + geom_line(aes(y=DO.mod), color="red")
  
})