
# Build a minimal, valid two-station data.frame. Defaults give a 5-minute
# timestep (0.0034722 days) and a 0.01-day travel time, so
# max_lag = round(0.01 / 0.0034722) = 3 timesteps of required upstream lead-in.
make_2station_data <- function(n=10, timestep_min=5, travel_time=0.01) {
  data.frame(
    solar.time = as.POSIXct("2050-06-01 00:00:00", tz="UTC") +
      as.difftime((seq_len(n) - 1) * timestep_min, units="mins"),
    DO.obs.up = rep(9, n),
    DO.sat.up = rep(10, n),
    DO.obs.down = rep(8.8, n),
    DO.sat.down = rep(9.9, n),
    light = rep(300, n),
    depth = rep(0.5, n),
    temp.water = rep(20, n),
    travel.time = rep(travel_time, n)
  )
}

test_that("mm_validate_data catches missing required columns", {
  dat <- dplyr::select(make_2station_data(), -DO.obs.up)
  expect_error(metab_2station(data=dat), "missing these columns")
})

test_that("travel.time <= 0 triggers an error", {
  dat <- make_2station_data(travel_time=0)
  expect_error(metab_2station(data=dat), "travel.time must be > 0")

  dat <- make_2station_data(travel_time=-0.01)
  expect_error(metab_2station(data=dat), "travel.time must be > 0")
})

test_that("travel.time >= 1 triggers an error with a units hint", {
  dat <- make_2station_data(travel_time=1)
  expect_error(metab_2station(data=dat), "travel.time must be < 1 day.*incorrect units")

  dat <- make_2station_data(travel_time=1.5)
  expect_error(metab_2station(data=dat), "travel.time must be < 1 day.*incorrect units")
})

test_that("insufficient lead-in data triggers an error", {
  # 2 rows but max_lag=3 timesteps of upstream lead-in are needed
  dat <- make_2station_data(n=2)
  expect_error(metab_2station(data=dat), "insufficient lead-in data")
})

test_that("a minimal valid data.frame passes all data-format checks", {
  dat <- make_2station_data(n=10)
  # metab_2station is a stub, so a fully valid data.frame should make it all
  # the way through the column/units/travel.time/lead-in checks and fail only
  # on the final "not yet implemented" stop -- not on any validation check.
  expect_error(metab_2station(data=dat), "metab_2station not yet implemented")
})


# mm_ts_prep_data() -----------------------------------------------------

# Build a two-day, unit-labeled data.frame with a known, traceable
# DO.obs.up/DO.sat.up series (sequential integers) so the shift can be
# checked by exact value, plus a leading lead-in block. 5-minute timestep
# (0.0034722 days) and 0.01-day travel.time give
# max_lag = round(0.01 / 0.0034722) = 3 lead-in timesteps.
# Day 1 = 10 rows (3 lead-in + 7 modeled), Day 2 = 7 rows (all modeled), so
# both modeled days end up with n_obs = 7 rows.
make_ts_data <- function(n_leadin=3, n_day1=10, n_day2=7, travel_time=0.01, unitted=FALSE) {
  n_total <- n_day1 + n_day2
  solar.time <- c(
    as.POSIXct("2050-06-01 00:00:00", tz="UTC") + as.difftime((seq_len(n_day1) - 1) * 5, units="mins"),
    as.POSIXct("2050-06-02 00:00:00", tz="UTC") + as.difftime((seq_len(n_day2) - 1) * 5, units="mins"))
  dat <- data.frame(
    solar.time = solar.time,
    DO.obs.up = seq_len(n_total),        # traceable: value == original row index
    DO.sat.up = seq_len(n_total) + 100,  # traceable, offset so it's distinguishable from DO.obs.up
    DO.obs.down = seq_len(n_total) + 1000, # traceable, offset so it's distinguishable from up/sat values
    DO.sat.down = rep(9.9, n_total),
    light = rep(300, n_total),
    depth = rep(0.5, n_total),
    temp.water = rep(20, n_total),
    travel.time = rep(travel_time, n_total)
  )
  if(unitted) {
    units_template <- get_units(mm_data(
      solar.time, DO.obs.up, DO.sat.up, DO.obs.down, DO.sat.down, light, depth, temp.water, travel.time))
    dat <- u(dat, unname(units_template[names(dat)]))
  }
  dat
}

test_that("upstream DO is shifted by the correct lag", {
  dat <- make_ts_data()
  out <- mm_ts_prep_data(dat)

  # max_lag=3, so modeled row i (original index i) uses upstream data from
  # original row (i - 3). day 1's 7 modeled rows are original rows 4:10, so
  # they pick up DO.obs.up from original rows 1:7; day 2's 7 modeled rows are
  # original rows 11:17, picking up DO.obs.up from original rows 8:14.
  expect_equal(out$DO_obs_up[,1], as.numeric(1:7))
  expect_equal(out$DO_obs_up[,2], as.numeric(8:14))
  expect_equal(out$DO_sat_up[,1], as.numeric(1:7) + 100)
  expect_equal(out$DO_sat_up[,2], as.numeric(8:14) + 100)
})

test_that("lead-in rows are excluded from the output matrices", {
  dat <- make_ts_data(n_leadin=3, n_day1=10, n_day2=7)
  out <- mm_ts_prep_data(dat)

  # 17 total rows in, 3 are lead-in-only, so 14 modeled rows should remain
  expect_equal(out$n_obs * out$n_days, nrow(dat) - 3)
  # none of the lead-in DO.obs.up values (1, 2, 3) should appear as a
  # DOWNSTREAM-paired value, i.e., the first modeled column should start at
  # the shifted value 1, not 1:3 appearing as downstream/lead-in rows
  expect_false(any(dat$DO.obs.down[1:3] %in% unlist(out$DO_obs_down)))
})

test_that("output matrices have n_obs x n_days dimensions", {
  dat <- make_ts_data()
  out <- mm_ts_prep_data(dat)

  expect_equal(out$n_obs, 7)
  expect_equal(out$n_days, 2)
  for(varname in c('DO_obs_up','DO_sat_up','DO_obs_down','DO_sat_down','light','depth','temp_water','travel_time')) {
    expect_equal(dim(out[[varname]]), c(7, 2), info=varname)
  }
})

test_that("all required Stan data block variables are present", {
  dat <- make_ts_data()
  out <- mm_ts_prep_data(dat)

  expected_names <- c(
    'n_obs','n_days','DO_obs_up','DO_sat_up','DO_obs_down','DO_sat_down',
    'light','depth','temp_water','travel_time','K600_lnorm_meanlog','K600_lnorm_sdlog')
  expect_true(all(expected_names %in% names(out)))

  # placeholder K600 lognormal priors, per PR D-3
  expect_equal(out$K600_lnorm_meanlog, log(3.48))
  expect_equal(out$K600_lnorm_sdlog, 0.5)
})

test_that("units are stripped from all numeric outputs", {
  dat <- make_ts_data(unitted=TRUE)
  expect_true(is.unitted(dat))

  out <- mm_ts_prep_data(dat)
  for(varname in c('DO_obs_up','DO_sat_up','DO_obs_down','DO_sat_down','light','depth','temp_water','travel_time')) {
    expect_false(is.unitted(out[[varname]]), info=varname)
  }
  expect_false(is.unitted(out$K600_lnorm_meanlog))
  expect_false(is.unitted(out$K600_lnorm_sdlog))
  expect_false(is.unitted(out$n_obs))
  expect_false(is.unitted(out$n_days))
})


# mm_parse_name() for two-station models ---------------------------------

test_that("mm_parse_name recognizes the b2_ prefix for two-station models", {
  parsed <- mm_parse_name('b2_np_oi_tr_plrckm.stan')

  expect_equal(parsed$type, 'bayes_2station')
  # the rest of the name is shared syntax with one-station bayes models and
  # should parse the same way regardless of the b vs. b2 prefix
  expect_equal(parsed$pool_K600, 'none')
  expect_true(parsed$err_obs_iid)
  expect_false(parsed$err_proc_acor)
  expect_false(parsed$err_proc_iid)
  expect_false(parsed$err_proc_GPP)
  expect_equal(parsed$ode_method, 'trapezoid')
  expect_equal(parsed$GPP_fun, 'linlight')
  expect_equal(parsed$ER_fun, 'constant')
  expect_equal(parsed$deficit_src, 'DO_mod')
  expect_equal(parsed$engine, 'stan')

  # a one-station name with the same suffix should still parse as plain 'bayes'
  expect_equal(mm_parse_name('b_np_oi_tr_plrckm.stan')$type, 'bayes')
})


# predict_DO.metab_2station stub -----------------------------------------

test_that("predict_DO dispatches to the metab_2station stub and errors as expected", {
  # metab_2station() is itself a stub (see tests above) and never returns a
  # fitted model object, so there's no way yet to construct a real
  # metab_2station instance. Mock the class attribute alone to exercise S3
  # dispatch of predict_DO() to predict_DO.metab_2station().
  mm <- structure(list(), class = c('metab_2station', 'metab_model'))

  expect_error(
    predict_DO(mm),
    "predict_DO for two-station models not yet implemented .* requires D-4 completion"
  )
})
