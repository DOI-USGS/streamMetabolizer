
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
  expect_error(metab_bayes_2s(data=dat), "missing these columns")
})

test_that("travel.time <= 0 triggers an error", {
  dat <- make_2station_data(travel_time=0)
  expect_error(metab_bayes_2s(data=dat), "travel.time must be > 0")

  dat <- make_2station_data(travel_time=-0.01)
  expect_error(metab_bayes_2s(data=dat), "travel.time must be > 0")
})

test_that("travel.time > 8/24 days (8 hours) triggers an error with a units/limit hint", {
  dat <- make_2station_data(travel_time=1)
  expect_error(metab_bayes_2s(data=dat), "travel.time must be <= 8/24 days.*incorrect units")

  dat <- make_2station_data(travel_time=1.5)
  expect_error(metab_bayes_2s(data=dat), "travel.time must be <= 8/24 days.*incorrect units")
})

test_that("insufficient lead-in data triggers an error", {
  # 2 rows but max_lag=3 timesteps of upstream lead-in are needed
  dat <- make_2station_data(n=2)
  expect_error(metab_bayes_2s(data=dat), "insufficient lead-in data")
})


# prepdata_bayes_2s() -----------------------------------------------------

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
  out <- prepdata_bayes_2s(dat)

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
  out <- prepdata_bayes_2s(dat)

  # 17 total rows in, 3 are lead-in-only, so 14 modeled rows should remain
  expect_equal(out$n_obs * out$n_days, nrow(dat) - 3)
  # none of the lead-in DO.obs.up values (1, 2, 3) should appear as a
  # DOWNSTREAM-paired value, i.e., the first modeled column should start at
  # the shifted value 1, not 1:3 appearing as downstream/lead-in rows
  expect_false(any(dat$DO.obs.down[1:3] %in% unlist(out$DO_obs_down)))
})

test_that("output matrices have n_obs x n_days dimensions", {
  dat <- make_ts_data()
  out <- prepdata_bayes_2s(dat)

  expect_equal(out$n_obs, 7)
  expect_equal(out$n_days, 2)
  for(varname in c('DO_obs_up','DO_sat_up','DO_obs_down','DO_sat_down','light','depth','temp_water','travel_time')) {
    expect_equal(dim(out[[varname]]), c(7, 2), info=varname)
  }
})

test_that("all required Stan data block variables are present", {
  dat <- make_ts_data()
  # K600_lnorm_meanlog/sdlog are owned by specs() (see PR D-6/I1); pass
  # distinctive marker values here to confirm prepdata_bayes_2s() just reads
  # them through from specs rather than computing its own defaults
  out <- prepdata_bayes_2s(dat, specs=list(K600_lnorm_meanlog=1.23, K600_lnorm_sdlog=4.56))

  expected_names <- c(
    'n_obs','n_days','DO_obs_up','DO_sat_up','DO_obs_down','DO_sat_down',
    'light','depth','temp_water','travel_time','K600_lnorm_meanlog','K600_lnorm_sdlog')
  expect_true(all(expected_names %in% names(out)))

  expect_equal(out$K600_lnorm_meanlog, 1.23)
  expect_equal(out$K600_lnorm_sdlog, 4.56)
})

test_that("units are stripped from all numeric outputs", {
  dat <- make_ts_data(unitted=TRUE)
  expect_true(is.unitted(dat))

  out <- prepdata_bayes_2s(dat, specs=list(K600_lnorm_meanlog=2.484907, K600_lnorm_sdlog=1.0))
  for(varname in c('DO_obs_up','DO_sat_up','DO_obs_down','DO_sat_down','light','depth','temp_water','travel_time')) {
    expect_false(is.unitted(out[[varname]]), info=varname)
  }
  expect_false(is.unitted(out$K600_lnorm_meanlog))
  expect_false(is.unitted(out$K600_lnorm_sdlog))
  expect_false(is.unitted(out$n_obs))
  expect_false(is.unitted(out$n_days))
})


# mm_lag_light_2s() ---------------------------------------------------------

test_that("mm_lag_light_2s computes the traceable within-day light proportion", {
  # 1 lead-in row + one 06:00-06:00 day (4 rows at a 6-hourly timestep, the
  # same day window mm_align_2s() uses); travel.time=6hr -> lag=1, so each
  # modeled row's window is itself plus the immediately preceding row
  solar.time <- as.POSIXct("2050-06-01 00:00:00", tz="UTC") +
    as.difftime((0:4) * 6, units="hours")
  light <- c(10, 1, 2, 3, 4)
  travel.time <- rep(0.25, 5)

  out <- mm_lag_light_2s(solar.time, light, travel.time)

  # row 1 (2050-06-01 00:00) lacks lead-in (shift_idx = 0) and is also the
  # lone row of the (irrelevant, partial) 06:00-day 2050-05-31; rows 2:5
  # (06:00, 12:00, 18:00, and the next day's 00:00) fill the 06:00-day
  # labeled 2050-06-01 and are hand-computed:
  # day total (rows 2:5) = 1+2+3+4 = 10
  # row2: sum(light[1:2])/10 = 11/10; row3: sum(light[2:3])/10 = 3/10
  # row4: sum(light[3:4])/10 = 5/10;  row5: sum(light[4:5])/10 = 7/10
  expect_true(is.na(out[1]))
  expect_equal(out[2:5], c(1.1, 0.3, 0.5, 0.7))
})

test_that("insufficient lead-in yields NA, not a value from a truncated window", {
  # 2 lead-in hours + one full 06:00-06:00 day (24 hourly rows, 06:00
  # 2050-06-01 through 06:00 2050-06-02)
  n <- 26
  solar.time <- as.POSIXct("2050-06-01 04:00:00", tz="UTC") +
    as.difftime((0:(n - 1)) * 1, units="hours")
  light <- rep(1, n)
  travel.time <- rep(2/24, n) # 2 hours -> lag=2

  out <- mm_lag_light_2s(solar.time, light, travel.time)

  # rows 1:2 (04:00, 05:00) have shift_idx < 1 (no lead-in); confirmed via
  # mm_lag_2s directly rather than reimplementing the lead-in test here
  has_leadin <- mm_lag_2s(solar.time, travel.time)$has_leadin
  expect_false(any(has_leadin[1:2]))
  expect_true(all(is.na(out[!has_leadin])))
  expect_false(any(is.na(out[has_leadin])))
})

test_that("an incomplete 06:00-06:00 day (including a trailing partial day) gets NA, not a distorted fraction", {
  # 2 lead-in hours + two full 06:00-06:00 days (24 hourly rows each) + a
  # 5-row trailing partial day; constant light makes every complete-day,
  # lead-in-eligible row's proportion the same constant, so any deviation
  # would reveal a partial day leaking into the sum
  n <- 55
  solar.time <- as.POSIXct("2050-06-01 04:00:00", tz="UTC") +
    as.difftime((0:(n - 1)) * 1, units="hours")
  light <- rep(1, n)
  travel.time <- rep(2/24, n) # 2 hours -> lag=2, window width 3

  out <- mm_lag_light_2s(solar.time, light, travel.time)

  # rows 1:2 lack lead-in; rows 3:50 fall in the two complete 06:00-days
  # (proportion = window width / day length = 3/24); rows 51:55 are the
  # partial trailing day
  expect_true(all(is.na(out[1:2])))
  expect_equal(out[3:50], rep(3/24, 48))
  expect_true(all(is.na(out[51:55])))
})

test_that("mm_align_2s() and mm_lag_light_2s() agree on day boundaries (regression test for #475 day-window mismatch)", {
  # A day mm_align_2s() marks complete must never have NA light (see the
  # day-sum-denominator roxygen section of mm_lag_light_2s()).
  #
  # Fixture: hourly timestep, 2050-06-01 04:00 through 2050-06-03 08:00.
  # mm_align_2s() calls 06:00-days 2050-06-01 and 2050-06-02 complete (24/24
  # rows each); 2050-06-01's day spans the early-morning hours of the
  # midnight-incomplete calendar day 2050-06-02 (only 00:00-08:00 present,
  # since the dataset ends there) -- exactly the boundary this test guards.
  solar.time <- seq(
    as.POSIXct("2050-06-01 04:00:00", tz="UTC"),
    as.POSIXct("2050-06-03 08:00:00", tz="UTC"),
    by="1 hour")
  n <- length(solar.time)
  light <- rep(100, n)
  travel.time <- rep(2/24, n)

  aln <- mm_align_2s(data.frame(solar.time=solar.time, travel.time=travel.time), max_travel_time_hours=10)
  expect_equal(sort(as.character(unique(aln$date))), c("2050-06-01", "2050-06-02"))

  light_lag <- mm_lag_light_2s(solar.time, light, travel.time)
  # every row mm_align_2s() considers part of a complete day must have a
  # real (non-NA) light proportion -- the boundaries now agree by default
  expect_false(any(is.na(light_lag[aln$keep])))

  # and prepdata_bayes_2s()'s defensive check (independent of how light was
  # computed) should stay silent on this now-consistent data
  dat <- data.frame(
    solar.time = solar.time,
    DO.obs.up = rep(9, n), DO.sat.up = rep(10, n),
    DO.obs.down = rep(8.8, n), DO.sat.down = rep(9.9, n),
    light = light_lag, depth = rep(0.5, n), temp.water = rep(20, n),
    travel.time = travel.time)
  expect_silent(prepdata_bayes_2s(dat, specs=list(K600_lnorm_meanlog=2.484907, K600_lnorm_sdlog=1.0), aln=aln))
})

test_that("prepdata_bayes_2s() errors clearly if a complete day's light is NA (defensive check)", {
  # Directly exercises the guard in prepdata_bayes_2s() -- independent of
  # mm_lag_light_2s() -- so it still catches a reintroduced day-window
  # mismatch (or any other source of NA light) even if the two functions
  # drift apart again in the future
  n <- 26
  solar.time <- as.POSIXct("2050-06-01 04:00:00", tz="UTC") +
    as.difftime((0:(n - 1)) * 1, units="hours")
  dat <- data.frame(
    solar.time = solar.time,
    DO.obs.up = seq_len(n), DO.sat.up = seq_len(n) + 100,
    DO.obs.down = seq_len(n) + 1000, DO.sat.down = rep(9.9, n),
    light = rep(300, n), depth = rep(0.5, n), temp.water = rep(20, n),
    travel.time = rep(2/24, n))

  # sanity: this fixture fits cleanly with valid light
  expect_silent(prepdata_bayes_2s(dat))

  # inject NA into a modeled row of the (otherwise complete) 06:00-day
  dat$light[10] <- NA
  expect_error(prepdata_bayes_2s(dat), "light is NA for.*day.*mm_align_2s.*complete")
})

test_that("mm_lag_light_2s errors clearly on an irregular timestep grid", {
  n <- 24
  solar.time <- as.POSIXct("2050-06-01 00:00:00", tz="UTC") +
    as.difftime((0:(n - 1)) * 1, units="hours")
  light <- rep(1, n)
  travel.time <- rep(2/24, n)

  # regular grid succeeds
  expect_silent(mm_lag_light_2s(solar.time, light, travel.time))

  # a single mid-series gap breaks regularity
  solar.time.gap <- solar.time
  solar.time.gap[10] <- solar.time.gap[10] + as.difftime(37, units="mins")
  expect_error(
    mm_lag_light_2s(solar.time.gap, light, travel.time),
    "regular timestep grid.*#475")
})


# mm_parse_name() for two-station models ---------------------------------

test_that("mm_parse_name recognizes the b2_ prefix for two-station models", {
  parsed <- mm_parse_name('b2_np_oi_tr_plrckm.stan')

  expect_equal(parsed$type, 'bayes_2s')
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


# mm_name() / mm_valid_names() / specs() for bayes_2s ---------------

test_that("mm_name(type='bayes_2s') returns the single two-station model name", {
  expect_equal(mm_name(type='bayes_2s'), 'b2_np_oi_tr_plrckm.stan')
})

test_that("mm_valid_names('bayes_2s') returns the single two-station model name", {
  expect_equal(mm_valid_names('bayes_2s'), 'b2_np_oi_tr_plrckm.stan')
})

test_that("specs(mm_name('bayes_2s')) has the expected params_in/params_out/split_dates", {
  sp <- specs(mm_name('bayes_2s'))

  expect_equal(
    sp$params_in,
    c('GPP_daily_mu', 'GPP_daily_sigma', 'ER_daily_mu', 'ER_daily_sigma',
      'K600_lnorm_meanlog', 'K600_lnorm_sdlog'))
  expect_equal(sp$params_out, c('GPP_daily', 'ER_daily', 'K600_daily', 'sigma', 'metab'))
  expect_false(sp$split_dates)
  expect_equal(sp$engine, 'stan')
})

test_that("specs()'s max_travel_time_hours default stays in sync with mm_max_travel_time_default", {
  # specs()'s default is a hand-copied literal, not a reference to
  # mm_max_travel_time_default (see the roxygen note in mm_lag_2s.R for why);
  # this is a tripwire in case the two are ever edited independently
  sp <- specs(mm_name('bayes_2s'))
  expect_equal(sp$max_travel_time_hours, mm_max_travel_time_default)
})


# metab_bayes_2s() fitting, predict_metab(), predict_DO() ------------------

# Subset two_station_example to just a few modeled days for a faster test
# fit. Naively slicing rows doesn't work: max_lag (the number of upstream
# lead-in rows required) is recomputed from whatever travel.time values are
# present in the slice, so an arbitrary row range can leave a partial first
# date once prepdata_bayes_2s() trims max_lag rows off the front -- the same
# lead-in-sizing logic used in data-raw/two_station_example.R is needed here
# too.
subset_2station_data <- function(full_data, n_modeled_days) {
  solar_time <- v(full_data$solar.time)
  timestep_days <- stats::median(as.numeric(diff(solar_time), units='days'))

  all_dates <- unique(as.Date(solar_time))
  modeled_dates <- all_dates[2:(1 + n_modeled_days)]
  modeled_start <- as.POSIXct(paste0(modeled_dates[1], ' 00:00:00'), tz='UTC')
  modeled_end <- as.POSIXct(paste0(modeled_dates[length(modeled_dates)], ' 23:45:00'), tz='UTC')

  candidate_start <- modeled_start - as.difftime(1, units='days')
  candidate <- full_data[solar_time >= candidate_start & solar_time <= modeled_end, ]
  max_lag <- max(round(v(candidate$travel.time) / timestep_days))
  lead_in_start <- modeled_start - as.difftime(max_lag * timestep_days, units='days')

  full_data[solar_time >= lead_in_start & solar_time <= modeled_end, ]
}

test_that("metab() fits a two-station model and predict_metab()/predict_DO() work", {
  skip_on_cran()
  skip_if_not_installed('rstan')

  small_dat <- subset_2station_data(two_station_example, n_modeled_days=3)

  sp <- specs(
    mm_name('bayes_2s'),
    n_chains=1, n_cores=1, burnin_steps=100, saved_steps=100, verbose=FALSE)

  mm <- metab(specs=sp, data=small_dat)
  expect_s4_class(mm, 'metab_bayes_2s')

  pm <- predict_metab(mm)
  expect_s3_class(pm, 'data.frame')
  expect_true(all(c('GPP','ER','K600') %in% names(pm)))
  expect_equal(nrow(pm), 3)

  pdo <- predict_DO(mm)
  expect_s3_class(pdo, 'data.frame')
  expect_true(all(c('DO.obs.down','DO.mod.down') %in% names(pdo)))
})

test_that("a failed Stan run (mode==2L) warns and continues, matching runstan_bayes()'s pattern, rather than erroring out", {
  skip_if_not_installed('rstan')

  # stand in for rstan::stan()'s return value on a failed run: only the
  # 'mode' slot is inspected by metab_bayes_2s() before deciding to skip
  # post-processing, so a minimal S4 object with that slot is sufficient
  setClass('fake_failed_stanfit', representation(mode='integer'))
  fake_stanfit <- methods::new('fake_failed_stanfit', mode=2L)
  testthat::local_mocked_bindings(stan=function(...) fake_stanfit, .package='rstan')

  dat <- make_ts_data()
  sp <- specs(
    mm_name('bayes_2s'),
    n_chains=1, n_cores=1, burnin_steps=10, saved_steps=10, verbose=FALSE)

  expect_warning(
    mm <- metab_bayes_2s(specs=sp, data=dat),
    'Modeling failed')

  expect_s4_class(mm, 'metab_bayes_2s')
  fit <- mm@fit
  expect_true(nrow(fit$daily) > 0)
  expect_true(all(is.na(fit$daily$GPP_daily_50pct)))
  expect_true(all(is.na(fit$daily$ER_daily_50pct)))
  expect_true(all(is.na(fit$daily$K600_daily_50pct)))
  expect_true(all(fit$daily$valid_day))
  expect_null(fit$inst)
  expect_equal(length(fit$errors), 0)
  expect_true(length(fit$warnings) > 0)
  expect_true(any(grepl('fake_failed_stanfit', fit$warnings)))
})
