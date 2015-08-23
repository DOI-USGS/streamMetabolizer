context("metab_model_helpers")

test_that("mm_data()", {
  expect_equal(mm_data(), get_data(metab_model()), info="mm_data() and metab_model()@data are the same")
  expect_equal(mm_data(local.time), mm_data()["local.time"], info="mm_data can be used for to select columns")
  expect_equal(mm_data(depth, temp.water, local.time), mm_data()[c("depth","temp.water","local.time")], info="mm_data can be used for to select columns")
})

test_that("mm_validate_data works", {
  expect_error(mm_validate_data(dplyr::select(mm_data(),-temp.water), "metab_mle"), "missing these columns: temp.water")
  expect_error(mm_validate_data(dplyr::mutate(mm_data(),temp.air=9), "metab_mle"), "should omit these extra columns: solar.time, temp.air")
  expect_error(mm_validate_data(dplyr::mutate(mm_data(),solar.time=NULL, temp.water=u(temp.water, "degF"), DO.obs=u(DO.obs, "mgO L^-1")), "metab_mle"), "unexpected units")
  expect_equal(mm_validate_data(dplyr::mutate(mm_data(),solar.time=NULL), "metab_mle"), dplyr::mutate(mm_data(),solar.time=NULL))
})

test_that("mm_is_valid_day works", {
  
  # use a subset of data from Bob
  library(plyr); library(dplyr)
  library(unitted)
  french <- streamMetabolizer:::load_french_creek()
  french$DO.sat <- calc_DO_at_sat(temp.water=french$temp.water, pressure.air=u(1000, "mb"))
  # assuming this is the french creek near buffalo
  french$utc.time <- convert_localtime_to_GMT(force_tz(french$local.time,"MST"))
  french$local.time <- force_tz(convert_GMT_to_localtime(french$utc.time, longitude=-106.753099, latitude=44.362594), "UTC")
  french$solar.time <- convert_GMT_to_solartime(french$utc.time, longitude=-106.753099, time.type='apparent solar')
  french$light <- convert_SW_to_PAR(calc_solar_insolation(solar.time=v(french$solar.time), latitude=44.362594, attach.units=TRUE))
  french <- french[c("local.time","DO.obs","DO.sat","depth","temp.water","light")]
  
  good_day <- u(filter(v(french), v(local.time) >= as.POSIXct("2012-08-24 22:30:00", tz="UTC"),
                       v(local.time) <= as.POSIXct("2012-08-26 06:00:00", tz="UTC")), get_units(french))
  bad_day <- u(filter(v(french), v(local.time) >= as.POSIXct("2012-08-28 22:30:00", tz="UTC"),
                       v(local.time) <= as.POSIXct("2012-08-30 06:00:00", tz="UTC")), get_units(french))
  
  # test and pass
  expect_equal(mm_is_valid_day(good_day), character(0))
  
  # test faulty timestep
  dateless_day <- good_day; dateless_day$local.time <- replace(dateless_day$local.time, 2:(nrow(dateless_day)-1), NA)
  expect_equal(mm_is_valid_day(dateless_day), c("can't measure timesteps", "NAs in local.time"))
  
  # test full_day
  expect_equal(mm_is_valid_day(good_day, day_start=0), "data don't start when expected")
  expect_equal(mm_is_valid_day(good_day, day_end=25), "data don't end when expected")
  
  # test timestep lengths
  irregular_day <- good_day[-c(3,20,99),]
  expect_equal(mm_is_valid_day(irregular_day), "uneven timesteps")
  
  # test column completeness
  expect_equal(mm_is_valid_day(bad_day), c("NAs in DO.obs","NAs in DO.sat","NAs in temp.water"))
})

test_that("mm_filter_valid_days works", {
  
  library(plyr); library(dplyr)
  library(unitted); library(lubridate)
  french <- streamMetabolizer:::load_french_creek()
  french$DO.sat <- calc_DO_at_sat(temp.water=french$temp.water, pressure.air=u(1000, "mb"))
  french$utc.time <- convert_localtime_to_GMT(force_tz(french$local.time,"MST"))
  french$local.time <- force_tz(convert_GMT_to_localtime(french$utc.time, longitude=-106.753099, latitude=44.362594), "UTC")
  french$solar.time <- convert_GMT_to_solartime(french$utc.time, longitude=-106.753099, time.type='apparent solar')
  french$light <- convert_SW_to_PAR(calc_solar_insolation(solar.time=v(french$solar.time), latitude=44.362594, attach.units=TRUE))
  french <- french[c("local.time","DO.obs","DO.sat","depth","temp.water","light")]
  
  mm_filter_valid_days(french)
})
