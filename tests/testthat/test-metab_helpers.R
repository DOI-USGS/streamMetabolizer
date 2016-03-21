context("metab_model_helpers")

test_that("mm_data works", {
  
  # runs and can be used to select columns
  expect_is(mm_data(), "unitted_data.frame")
  expect_equal(nrow(mm_data()), 1)
  expect_equivalent(mm_data(solar.time), mm_data()["solar.time"])
  expect_equivalent(mm_data(depth, temp.water, solar.time), mm_data()[c("depth","temp.water","solar.time")])
  
  # 'optional' attribute is set sensibly
  expect_equal(attr(mm_data(NULL), 'optional'), 'all')
  expect_equal(attr(mm_data(depth, temp.water, solar.time), 'optional'), 'none')
  expect_equal(attr(mm_data(solar.time, DO.obs, optional='DO.obs'), 'optional'), 'DO.obs')
  expect_equal(attr(mm_data(solar.time, DO.obs, optional=c('DO.obs','solar.time')), 'optional'), 'all')

  # 'optional' attribute is checked
  expect_error(mm_data(solar.time, DO.obs, optional='DO.sat'), "'arg' should be one of")
  expect_error(mm_data(solar.time, DO.obs, optional=c('DO.obs','all')), "should be length 1")
})

test_that("mm_validate_data works", {
  
  # runs and accepts the defaults without errors
  ignore <- mm_validate_data(eval(formals(metab_mle)$data), eval(formals(metab_mle)$data_daily), 'metab_mle')
  ignore <- mm_validate_data(eval(formals(metab_night)$data), eval(formals(metab_night)$data_daily), 'metab_night')
  ignore <- mm_validate_data(eval(formals(metab_bayes)$data), eval(formals(metab_bayes)$data_daily), 'metab_bayes')
  ignore <- mm_validate_data(eval(formals(metab_Kmodel)$data), eval(formals(metab_Kmodel)$data_daily), 'metab_Kmodel')
  
  # accepts NULL for the fully optional data.frames
  ignore <- mm_validate_data(eval(formals(metab_mle)$data), NULL, 'metab_mle')
  ignore <- mm_validate_data(eval(formals(metab_night)$data), NULL, 'metab_night')
  ignore <- mm_validate_data(eval(formals(metab_bayes)$data), NULL, 'metab_bayes')
  ignore <- mm_validate_data(NULL, eval(formals(metab_Kmodel)$data_daily), 'metab_Kmodel')
  
  # returns a list
  val_out <- mm_validate_data(eval(formals(metab_mle)$data), eval(formals(metab_mle)$data_daily), 'metab_mle')
  expect_is(val_out, 'list')
  expect_equal(names(val_out), c('data','data_daily'))
  expect_is(val_out[[1]], 'data.frame')
  
  # notices missing, extra, badly unitted columns in data; accepts non-unitted data
  library(unitted)
  ok_data <- eval(formals(metab_mle)$data)
  expect_error(mm_validate_data(NULL, NULL, "metab_mle"), "data is NULL but required")
  expect_error(mm_validate_data(data.frame(), NULL, "metab_mle"), "data is missing these columns: solar.time, DO.obs, DO.sat, depth, temp.water, light")
  expect_error(mm_validate_data(dplyr::select(ok_data,-temp.water), NULL, "metab_mle"), "data is missing these columns: temp.water")
  expect_error(mm_validate_data(dplyr::mutate(ok_data,temp.air=9), mm_data('temp.air'), "metab_mle"), "data should omit these extra columns: temp.air")
  expect_error(mm_validate_data(dplyr::mutate(ok_data,temp.water=u(temp.water, "degF"), DO.obs=u(DO.obs, "mgO L^-1")), NULL, "metab_mle"), "unexpected units in data:")
  expect_is(mm_validate_data(unitted::v(ok_data), NULL, "metab_mle"), 'list')

  # notices missing, extra, badly unitted columns in data_daily
  ok_data_daily <- eval(formals(metab_mle)$data_daily)
  expect_is(mm_validate_data(ok_data, NULL, "metab_mle"), 'list')
  expect_error(mm_validate_data(ok_data, data.frame(), "metab_mle"), "found 0 possible timestamp columns")
  expect_error(mm_validate_data(ok_data, dplyr::mutate(ok_data_daily, temp.air=9), "metab_mle"), "data_daily should omit these extra columns: temp.air")
  expect_error(mm_validate_data(ok_data, dplyr::mutate(ok_data_daily, K600=u(K600, "mgO L^-1")), "metab_mle"), "unexpected units in data_daily:")
  expect_is(mm_validate_data(ok_data, unitted::v(ok_data_daily), "metab_mle"), 'list') # should we accept units for one but not other? for now, we do.
  
})

test_that("mm_is_valid_day works", {
  
  # use a subset of data from Bob
  library(unitted)
  french <- streamMetabolizer:::load_french_creek()
 
  good_day <- u(dplyr::filter(v(french), solar.time >= as.POSIXct("2012-08-24 22:30:00", tz="UTC"),
                              solar.time <= as.POSIXct("2012-08-26 06:00:00", tz="UTC")), get_units(french))
  bad_day <- dplyr::mutate(
    good_day,
    DO.obs=replace(DO.obs, 40, u(NA, get_units(DO.obs))),
    DO.sat=replace(DO.sat, 51, u(NA, get_units(DO.sat))),
    temp.water=replace(temp.water, 20:42, u(NA, get_units(temp.water))))
  
  # test and pass
  expect_true(mm_is_valid_day(good_day, day_start=-1.5, day_end=30))
  
  # test faulty timestep
  dateless_day <- good_day; dateless_day$solar.time <- replace(dateless_day$solar.time, 2:(nrow(dateless_day)-1), NA)
  expect_equal(mm_is_valid_day(dateless_day), c("NAs in solar.time"))
  
  # test full_day
  expect_equal(mm_is_valid_day(good_day, day_start=-10, day_end=30), "data don't start when expected")
  expect_equal(mm_is_valid_day(good_day, day_start=0, day_end=30), "data don't start when expected")
  expect_equal(mm_is_valid_day(good_day, day_start=-1.5, day_end=25), "data don't end when expected")
  expect_equal(mm_is_valid_day(good_day, day_start=-1.5, day_end=35), "data don't end when expected")
  
  # test timestep lengths
  irregular_day <- good_day[-c(3,20,99),]
  expect_equal(mm_is_valid_day(irregular_day, day_start=-1.5, day_end=30), "uneven timesteps")
  
  # test column completeness
  expect_equal(mm_is_valid_day(bad_day, day_start=-1.5, day_end=30), c("NAs in DO.obs","NAs in DO.sat","NAs in temp.water"))

})

test_that("mm_filter_valid_days works", {
  
  library(dplyr); library(unitted)
  
  # catch missorted data
  french <- data_metab('10', res='30', flaws='missorted', day_start=6, day_end=14)
  expect_error(mm_filter_valid_days(french, day_start=10, day_end=12), "min timestep is <= 0")
  
  # filter to specified hours
  french <- data_metab('10', res='30', day_start=6, day_end=14)
  french_filt1 <- mm_filter_valid_days(french, day_start=10, day_end=12)
  expect_equal(names(french_filt1), c('data','data_daily','removed'))
  expect_equal(nrow(french_filt1$data), 4*10)
  
  # filter to specified days
  french_daily <- data.frame(date=as.Date(sprintf("2012-09-%2d",15:30)), K600=7)
  french_filt2 <- mm_filter_valid_days(french, data_daily=french_daily, day_start=10, day_end=12)
  expect_equal(nrow(french_filt2$data), 10*4)
  expect_equal(nrow(french_filt2$removed), 6)
  expect_equal(sort(as.Date(unique(c(french_filt2$removed$date, french_filt2$data_daily$date)))), french_daily$date)
})

test_that("mm_filter_dates works", {
  
  start_time <- Sys.time()
  start_date <- as.Date(start_time)
  udat <- data.frame(solar.time=start_time + as.difftime(1:100, units='hours'), value=1:100)
  ddat <- data.frame(date=start_date + as.difftime(1:100, units='days'), value=1:100)
  # no filter with defaults
  expect_equal(streamMetabolizer:::mm_filter_dates(udat), udat)
  expect_equal(streamMetabolizer:::mm_filter_dates(ddat), ddat)
  # dates are inclusive
  expect_equal(nrow(streamMetabolizer:::mm_filter_dates(udat, date_start=start_date+as.difftime(1, units='days'), date_end=start_date+as.difftime(1, units='days'))), 24)
  expect_equal(nrow(streamMetabolizer:::mm_filter_dates(udat, date_start=start_date+as.difftime(1, units='days'), date_end=start_date+as.difftime(2, units='days'))), 48)
  expect_equal(nrow(streamMetabolizer:::mm_filter_dates(udat, date_start=start_date+as.difftime(2, units='days'), date_end=start_date+as.difftime(1, units='days'))), 0)
  expect_equal(nrow(streamMetabolizer:::mm_filter_dates(ddat, date_start=start_date+as.difftime(20, units='days'), date_end=start_date+as.difftime(20, units='days'))), 1)
  expect_equal(nrow(streamMetabolizer:::mm_filter_dates(ddat, date_start=start_date+as.difftime(22, units='days'), date_end=start_date+as.difftime(28, units='days'))), 7)
  expect_equal(nrow(streamMetabolizer:::mm_filter_dates(ddat, date_start=start_date+as.difftime(22, units='days'), date_end=start_date+as.difftime(8, units='days'))), 0)
})
