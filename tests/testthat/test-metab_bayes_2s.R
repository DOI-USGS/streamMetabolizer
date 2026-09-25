test_that("mm_validate_data catches missing required columns", {
  dat <- dplyr::select(make_2station_data(), -DO.obs.up)
  expect_error(mm_align_data_2s(dat), "missing these columns")
})

test_that("travel.time <= 0 triggers an error", {
  dat <- make_2station_data(travel_time=0)
  expect_error(mm_align_data_2s(dat), "travel.time must be > 0")

  dat <- make_2station_data(travel_time=-0.01)
  expect_error(mm_align_data_2s(dat), "travel.time must be > 0")
})

test_that("metab_bayes_2s() rejects data that has not been aligned, naming both ways to align it", {
  dat <- make_2station_data()
  expect_error(metab_bayes_2s(data=dat), "mm_align_data_2s\\(\\).*mm_as_aligned_2s\\(\\)")
  # a class-stripping conversion loses the aligned status too
  expect_error(
    metab_bayes_2s(data=as.data.frame(suppressMessages(mm_align_data_2s(dat)))),
    "must be aligned first")
})

test_that("metab_bayes_2s() rejects NA travel.time in aligned data with the affected dates", {
  # an aligned frame edited after aligning
  aligned <- suppressMessages(mm_align_data_2s(make_2station_data()))
  bad_date <- aligned$date[1]
  aligned$travel.time[3] <- NA
  err <- expect_error(
    metab_bayes_2s(specs=revise(fast_2station_specs(), day_tests=c()), data=aligned),
    paste0("travel.time is NA on ", bad_date, "; fill it or remove those days before fitting"))
  expect_false(grepl("mm_[a-z_0-9]+\\(|missing value", conditionMessage(err)))
})

test_that("mm_format_days_hours() shows days with hours in parentheses, singular at exactly one", {
  fmt <- streamMetabolizer:::mm_format_days_hours
  expect_equal(fmt(0.5), "0.5 days (12 hours)")
  expect_equal(fmt(10/24), "0.42 days (10 hours)")
  expect_equal(fmt(15/24), "0.63 days (15 hours)")
  expect_equal(fmt(1/24), "0.04 days (1 hour)")
  expect_equal(fmt(1), "1 day (24 hours)")
  expect_equal(fmt(0.125), "0.13 days (3 hours)")
})

test_that("mm_align_2s drops a day whose travel.time exceeds the ceiling, leaving other days intact", {
  dat <- make_2day_ceiling_data()

  alignment <- expect_message(
    mm_align_2s(dat),
    "dropping 1 day\\(s\\) whose travel.time exceeds the 0\\.42 days \\(10 hours\\) ceiling: 2050-06-02 \\(0\\.63 days \\(15 hours\\)\\)")

  # the offending day is gone entirely, not merely marked invalid; the good
  # day (2050-06-01) is unaffected
  expect_equal(as.character(unique(alignment$date)), "2050-06-01")
  expect_equal(alignment$n_days, 1)
  expect_equal(alignment$n_obs, 288)
})

test_that("metab_bayes_2s() drops a day exceeding the travel-time ceiling and fits the remaining day normally", {
  skip_on_cran()
  skip_if_not_installed('rstan')

  dat <- make_2day_ceiling_data()
  sp <- specs(
    mm_name('bayes_2s'),
    n_chains=1, n_cores=1, burnin_steps=100, saved_steps=100, verbose=FALSE)

  # the drop is announced when aligning, and reported again by the fit
  aligned <- expect_message(
    mm_align_data_2s(dat),
    "dropping 1 day\\(s\\) whose travel.time exceeds the 0\\.42 days \\(10 hours\\) ceiling: 2050-06-02")
  expect_equal(as.character(attr(aligned, 'removed')$date), "2050-06-02")
  mm <- suppressMessages(metab_bayes_2s(specs=sp, data=aligned))

  # day 2 was excluded before fitting, so it is reported as an invalid day
  # naming the ceiling rather than as a failed fit; day 1 is unaffected
  expect_length(mm@fit$errors, 0)
  pm <- predict_metab(mm)
  expect_equal(nrow(pm), 2)
  expect_equal(as.character(pm$date), c("2050-06-01", "2050-06-02"))
  expect_true(all(is.finite(unlist(pm[1, c('GPP','ER','K600')]))))
  expect_true(all(is.na(unlist(pm[2, c('GPP','ER','K600')]))))
  expect_equal(mm@fit$daily$valid_day, c(TRUE, FALSE))
  expect_match(mm@fit$daily$errors[2], "travel.time exceeds the 0\\.42 days \\(10 hours\\) ceiling")
})


# prepdata_bayes_2s() -----------------------------------------------------

test_that("upstream DO is shifted by the correct lag", {
  dat <- make_ts_data()
  out <- prepdata_bayes_2s(suppressMessages(mm_align_data_2s(dat)))

  # max_lag=3, so modeled row i (original index i) uses upstream data from
  # original row (i - 3). day 1's 24 modeled rows are original rows 4:27, so
  # they pick up DO.obs.up from original rows 1:24; day 2's 24 modeled rows
  # are original rows 28:51, picking up DO.obs.up from original rows 25:48.
  expect_equal(out$DO_obs_up[,1], as.numeric(1:24))
  expect_equal(out$DO_obs_up[,2], as.numeric(25:48))
  expect_equal(out$DO_sat_up[,1], as.numeric(1:24) + 100)
  expect_equal(out$DO_sat_up[,2], as.numeric(25:48) + 100)
})

test_that("lead-in rows are excluded from the output matrices", {
  dat <- make_ts_data(n_leadin=3, n_day1=24, n_day2=24)
  out <- prepdata_bayes_2s(suppressMessages(mm_align_data_2s(dat)))

  # 51 total rows in, 3 are lead-in-only, so 48 modeled rows should remain
  expect_equal(out$n_obs * out$n_days, nrow(dat) - 3)
  # none of the lead-in DO.obs.up values (1, 2, 3) should appear as a
  # DOWNSTREAM-paired value, i.e., the first modeled column should start at
  # the shifted value 1, not 1:3 appearing as downstream/lead-in rows
  expect_false(any(dat$DO.obs.down[1:3] %in% unlist(out$DO_obs_down)))
})

test_that("output matrices have n_obs x n_days dimensions", {
  dat <- make_ts_data()
  out <- prepdata_bayes_2s(suppressMessages(mm_align_data_2s(dat)))

  expect_equal(out$n_obs, 24)
  expect_equal(out$n_days, 2)
  for(varname in c('DO_obs_up','DO_sat_up','DO_obs_down','DO_sat_down','light','depth','temp_water','travel_time')) {
    expect_equal(dim(out[[varname]]), c(24, 2), info=varname)
  }
})

test_that("all required Stan data block variables are present", {
  dat <- make_ts_data()
  # K600_lnorm_meanlog/sdlog are owned by specs() (see PR D-6/I1); pass
  # distinctive marker values here to confirm prepdata_bayes_2s() just reads
  # them through from specs rather than computing its own defaults
  out <- prepdata_bayes_2s(
    suppressMessages(mm_align_data_2s(dat)), specs=list(K600_lnorm_meanlog=1.23, K600_lnorm_sdlog=4.56))

  expected_names <- c(
    'n_obs','n_days','DO_obs_up','DO_sat_up','DO_obs_down','DO_sat_down',
    'light','depth','temp_water','travel_time','K600_lnorm_meanlog','K600_lnorm_sdlog')
  expect_true(all(expected_names %in% names(out)))

  expect_equal(out$K600_lnorm_meanlog, 1.23)
  expect_equal(out$K600_lnorm_sdlog, 4.56)
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

  alignment <- mm_align_2s(data.frame(solar.time=solar.time, travel.time=travel.time), max_travel_time_days=10/24)
  expect_equal(sort(as.character(unique(alignment$date))), c("2050-06-01", "2050-06-02"))

  light_lag <- mm_lag_light_2s(solar.time, light, travel.time)
  # every row mm_align_2s() considers part of a complete day must have a
  # real (non-NA) light proportion -- the boundaries now agree by default
  expect_false(any(is.na(light_lag[alignment$keep])))

  # and so the aligned frame built from that light has none either
  dat <- data.frame(
    solar.time = solar.time,
    DO.obs.up = rep(9, n), DO.sat.up = rep(10, n),
    DO.obs.down = rep(8.8, n), DO.sat.down = rep(9.9, n),
    light = light_lag, depth = rep(0.5, n), temp.water = rep(20, n),
    travel.time = travel.time)
  aligned <- suppressMessages(mm_align_data_2s(dat, max_travel_time_days=10/24))
  expect_false(anyNA(aligned$light))
  expect_silent(prepdata_bayes_2s(aligned, specs=list(K600_lnorm_meanlog=2.484907, K600_lnorm_sdlog=1.0)))
})

test_that("a joint fit stops on NA data that day_tests let through, naming the column and the day", {
  # one complete hourly 06:00-day, 2 hours of travel time
  n <- 26
  solar.time <- as.POSIXct("2050-06-01 04:00:00", tz="UTC") +
    as.difftime((0:(n - 1)) * 1, units="hours")
  dat <- data.frame(
    solar.time = solar.time,
    DO.obs.up = seq_len(n), DO.sat.up = seq_len(n) + 100,
    DO.obs.down = seq_len(n) + 1000, DO.sat.down = rep(9.9, n),
    light = rep(300, n), depth = rep(0.5, n), temp.water = rep(20, n),
    travel.time = rep(2/24, n))
  sp <- revise(fast_2station_specs(), day_tests=c())

  # a modeled row's own column
  dat_light <- dat
  dat_light$light[10] <- NA
  err <- expect_error(
    suppressMessages(metab_bayes_2s(specs=sp, data=suppressMessages(mm_align_data_2s(dat_light)))),
    "NAs in light on 1 day\\(s\\): 2050-06-01")
  expect_match(conditionMessage(err), "drop or fill those days before fitting")
  expect_false(grepl("mm_[a-z_0-9]+\\(|prepdata", conditionMessage(err)))

  # an upstream column, drawn from a row that sits before the day starts
  dat_up <- dat
  dat_up$DO.obs.up[2] <- NA
  expect_error(
    suppressMessages(metab_bayes_2s(specs=sp, data=suppressMessages(mm_align_data_2s(dat_up)))),
    "NAs in DO.obs.up on 1 day\\(s\\): 2050-06-01")
})

test_that("mm_lag_light_2s errors clearly when solar.time isn't on a snap-to-bin grid", {
  n <- 24
  solar.time <- as.POSIXct("2050-06-01 00:00:00", tz="UTC") +
    as.difftime((0:(n - 1)) * 1, units="hours")
  light <- rep(1, n)
  travel.time <- rep(2/24, n)

  # already on the grid succeeds
  expect_silent(mm_lag_light_2s(solar.time, light, travel.time))

  # a single off-grid timestamp (not a real gap -- an unsnapped offset)
  # fails mm_lag_2s()'s snap-to-bin precondition check
  solar.time.offgrid <- solar.time
  solar.time.offgrid[10] <- solar.time.offgrid[10] + as.difftime(37, units="mins")
  expect_error(
    mm_lag_light_2s(solar.time.offgrid, light, travel.time),
    "not on a regular timestep grid.*mm_format_data_2s")
})

test_that("a real mid-series gap (on-grid, but a bin missing) is handled, not rejected", {
  # 24 hourly rows with row 10 (bin 09:00) dropped entirely -- a true gap,
  # still exactly on the hourly grid, unlike the off-grid case above.
  # travel.time=4hr -> lag=4, so the gap is far enough back to fall both at
  # a window boundary (for one row) and strictly inside a window (for
  # another), exercising both of mm_lag_2s()'s consumers' failure modes
  n <- 24
  solar.time <- as.POSIXct("2050-06-01 00:00:00", tz="UTC") +
    as.difftime((0:(n - 1)) * 1, units="hours")
  light <- rep(1, n)
  travel.time <- rep(4/24, n)

  keep <- seq_len(n) != 10 # drop the 09:00 bin
  solar.time.gap <- solar.time[keep]
  light.gap <- light[keep]
  travel.time.gap <- travel.time[keep]

  # no error: a missing bin is not a grid violation
  expect_silent(lagged <- mm_lag_2s(solar.time.gap, travel.time.gap))

  # point-lookup case: 13:00's target (13:00 - 4h = 09:00) is exactly the
  # missing bin -- no real upstream match, so no silent (wrong) pairing
  row_13 <- which(solar.time.gap == as.POSIXct("2050-06-01 13:00:00", tz="UTC"))
  expect_false(lagged$has_leadin[row_13])

  # window-sum case: 12:00's target (08:00) exists, so this row still has
  # lead-in, but its window (08:00 through 12:00) spans the missing 09:00
  # bin. The array-position range shift_idx:i still correctly resolves to
  # only the 4 rows that actually exist (08:00, 10:00, 11:00, 12:00), not
  # 5 -- a gap inside the window contributes fewer terms, not NA and not a
  # wrong pairing
  row_12 <- which(solar.time.gap == as.POSIXct("2050-06-01 12:00:00", tz="UTC"))
  expect_true(lagged$has_leadin[row_12])
  window <- light.gap[lagged$shift_idx[row_12]:row_12]
  expect_equal(length(window), 4)
  expect_equal(
    solar.time.gap[lagged$shift_idx[row_12]:row_12],
    as.POSIXct(c("2050-06-01 08:00:00", "2050-06-01 10:00:00", "2050-06-01 11:00:00", "2050-06-01 12:00:00"), tz="UTC"))
})

test_that("mm_lag_light_2s runs without erroring on the real, gappy two_station_raw_example (regression test for #475)", {
  # previously blocked outright by the old regularity guard (any gap ->
  # hard error); now gap-safe via timestep-bin matching in mm_lag_2s()
  data("two_station_raw_example", envir=environment())
  up <- two_station_raw_example$upstream
  up_ts <- as.POSIXct(v(up$timestamp), tz="UTC")
  travel.time <- rep(0.1, length(up_ts)) # real travel.time lives on downstream; a stand-in is fine here

  expect_no_error(mm_lag_light_2s(up_ts, rep(1, length(up_ts)), travel.time))
})

test_that("mm_align_2s drops a day mixing gap-affected and clean rows wholesale, leaving a fully clean day intact", {
  # lead-in (4 rows before day 1's 06:00) + day 1 (96 rows @ 15 min, one
  # dropped from the middle -- a real gap) + day 2 (96 rows, fully clean).
  # travel.time=2 timesteps (30 min) -> lag=2
  timestep_min <- 15
  leadin <- as.POSIXct("2050-06-01 05:00:00", tz="UTC") + as.difftime((0:3) * timestep_min, units="mins")
  day1 <- as.POSIXct("2050-06-01 06:00:00", tz="UTC") + as.difftime((0:95) * timestep_min, units="mins")
  day2 <- as.POSIXct("2050-06-02 06:00:00", tz="UTC") + as.difftime((0:95) * timestep_min, units="mins")
  solar.time <- c(leadin, day1, day2)
  solar.time <- solar.time[solar.time != day1[48]] # drop one row from day 1's middle

  travel.time <- rep(2 * timestep_min / 1440, length(solar.time))
  data <- data.frame(solar.time=solar.time, travel.time=travel.time)

  alignment <- mm_align_2s(data, max_travel_time_days=10/24)

  # day 1 (gap-affected: the dropped row plus the one row whose target bin
  # was the dropped row) is 94/96 complete and gets dropped wholesale, even
  # though most of its rows individually had real upstream matches; day 2
  # (fully clean) is unaffected and modeled in full
  expect_equal(as.character(unique(alignment$date)), "2050-06-02")
  expect_equal(alignment$n_days, 1)
  expect_equal(alignment$n_obs, 96)
})


# mm_snap_to_bin_2s() --------------------------------------------------------

test_that("mm_snap_to_bin_2s snaps two phase-shifted patterns onto one common grid", {
  # a 15-min nominal grid, plus two synthetic deployments each internally
  # regular but offset from that grid -- one at :01/:16/:31/:46 (matching
  # the design doc's example pattern), the other at :04/:19/:34/:49
  marks <- as.POSIXct("2050-06-01 00:00:00", tz="UTC") + as.difftime((0:19) * 15, units="mins")
  series_a <- marks + as.difftime(1, units="mins")
  series_b <- marks + as.difftime(4, units="mins")

  snapped_a <- mm_snap_to_bin_2s(series_a)
  snapped_b <- mm_snap_to_bin_2s(series_b)

  # both phase patterns round (well within the +/-7.5 min half-timestep
  # tolerance) onto the exact same grid, and thus onto each other
  expect_equal(snapped_a, marks)
  expect_equal(snapped_b, marks)
  expect_equal(snapped_a, snapped_b)
})

test_that("mm_snap_to_bin_2s turns a true gap between deployments into empty bins, not a merge", {
  # deployment 1: 5 rows, 15-min cadence, phase +1 min
  seg1 <- as.POSIXct("2050-06-01 00:01:00", tz="UTC") + as.difftime((0:4) * 15, units="mins")
  # a true gap of 3 hours (12 timesteps), then deployment 2: 5 rows, 15-min
  # cadence, phase +4 min -- concatenated deployments, per the design doc
  seg2 <- seg1[5] + as.difftime(3, units="hours") + as.difftime(4, units="mins") +
    as.difftime((0:4) * 15, units="mins")
  solar.time <- c(seg1, seg2)

  snapped <- mm_snap_to_bin_2s(solar.time)

  # snapping only rounds existing rows -- it neither drops nor fabricates any
  expect_equal(length(snapped), length(solar.time))
  expect_false(any(duplicated(snapped)))

  # the true gap survives as a 3-hour (12-timestep) jump between the last row
  # of deployment 1 and the first row of deployment 2, not a 1-timestep step
  expect_equal(as.numeric(difftime(snapped[6], snapped[5], units="hours")), 3)

  # the bins inside that gap are simply absent -- not silently interpolated
  gap_bins <- as.POSIXct("2050-06-01 00:00:00", tz="UTC") +
    as.difftime(c(75, 90, 195), units="mins") # 01:15, 01:30, ..., 03:45
  expect_false(any(gap_bins %in% snapped))
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


# metab_bayes_2s() fitting, predict_metab(), predict_DO() ------------------

test_that("metab() fits a two-station model and predict_metab()/predict_DO() work", {
  skip_on_cran()
  skip_if_not_installed('rstan')

  small_dat <- subset_2station_data(two_station_example, n_modeled_days=3)

  sp <- specs(
    mm_name('bayes_2s'),
    n_chains=1, n_cores=1, burnin_steps=100, saved_steps=100, verbose=FALSE)

  mm <- suppressMessages(metab(specs=sp, data=suppressMessages(mm_align_data_2s(small_dat))))
  expect_s4_class(mm, 'metab_bayes_2s')

  pm <- predict_metab(mm)
  expect_s3_class(pm, 'data.frame')
  expect_true(all(c('GPP','ER','K600') %in% names(pm)))
  # 3 modeled days, plus the partial day formed by the lead-in block, which
  # is reported as an invalid day rather than omitted
  expect_equal(nrow(pm), 4)
  expect_equal(mm@fit$daily$valid_day, c(FALSE, TRUE, TRUE, TRUE))
  expect_equal(sum(is.finite(pm$GPP)), 3)

  pdo <- predict_DO(mm)
  expect_s3_class(pdo, 'data.frame')
  expect_true(all(c('DO.obs.down','DO.mod.down') %in% names(pdo)))
})

test_that("a failed Stan run (mode==2L) warns and continues, matching runstan_bayes()'s pattern, rather than erroring out", {
  skip_if_not_installed('rstan')

  # stand in for rstan::sampling()'s return value on a failed run: only the
  # 'mode' slot is inspected by metab_bayes_2s() before deciding to skip
  # post-processing, so a minimal S4 object with that slot is sufficient.
  # mm_compile_stan_model() is also mocked so this test doesn't depend on a
  # real compile or on the .stanrds cache file's state on disk.
  setClass('fake_failed_stanfit', representation(mode='integer'))
  fake_stanfit <- methods::new('fake_failed_stanfit', mode=2L)
  testthat::local_mocked_bindings(
    mm_compile_stan_model=function(...) list(stan_mobj=NULL, compile_time=system.time({}), compile_log=NULL))
  testthat::local_mocked_bindings(sampling=function(...) fake_stanfit, .package='rstan')

  dat <- make_ts_data()
  sp <- specs(
    mm_name('bayes_2s'),
    n_chains=1, n_cores=1, burnin_steps=10, saved_steps=10, verbose=FALSE)

  expect_warning(
    mm <- metab_bayes_2s(specs=sp, data=suppressMessages(mm_align_data_2s(dat))),
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


# bayes_perday_2s() per-day fitting ---------------------------------------

test_that("a per-day slice of aligned data preps the same Stan matrices as the joint fit's column", {
  # The correctness test for the slicing itself, independent of Stan: each
  # day's one-column data list must match that day's column of the all-days
  # data list element for element.
  sp <- fast_2station_specs()
  dat <- subset_2station_days(two_station_example, 4)
  aligned <- suppressMessages(mm_align_data_2s(dat))
  joint <- prepdata_bayes_2s(aligned, specs=sp)

  # every day's upstream values reach back before its own first row, so each
  # slice must carry values drawn from the previous day's rows
  alignment <- suppressMessages(mm_align_2s(v(dat)))
  for(dt in as.list(unique(alignment$date))) {
    rows <- which(alignment$date == dt)
    expect_lt(min(alignment$shift_idx[rows]), min(alignment$keep[rows]))
  }

  dates <- unique(aligned$date)
  mats <- c('DO_obs_up','DO_sat_up','DO_obs_down','DO_sat_down',
            'light','depth','temp_water','travel_time')
  for(i in seq_along(dates)) {
    day <- prepdata_bayes_2s(aligned[aligned$date == dates[i], ], specs=sp)
    expect_equal(day$n_days, 1L)
    expect_equal(day$n_obs, joint$n_obs)
    for(m in mats) {
      expect_identical(as.vector(day[[m]]), as.vector(joint[[m]][, i]),
                       info=paste(m, 'on', dates[i]))
    }
  }
})

test_that("bayes_perday_2s() fits one day at a time and returns the joint fit's output shape", {
  skip_on_cran()
  skip_if_not_installed('rstan')

  sp <- fast_2station_specs()
  dat <- suppressMessages(mm_align_data_2s(subset_2station_days(two_station_example, 3)))

  res <- suppressMessages(bayes_perday_2s(dat, specs=sp))

  # (a) one result per day
  expect_equal(nrow(res$daily), 3)
  expect_length(res$dates_fit, 3)
  expect_length(res$dates_failed, 0)
  expect_equal(res$daily$date, unique(dat$date))
  expect_length(res$mcmcs, 3)

  # (b) plausible per-day estimates: GPP positive, ER negative, K600 positive
  expect_true(all(res$daily$GPP_daily_50pct > 0))
  expect_true(all(res$daily$ER_daily_50pct < 0))
  expect_true(all(res$daily$K600_daily_50pct > 0))
  expect_true(all(is.finite(res$daily$GPP_daily_50pct)))

  # inst covers every modeled timestep of every day, in solar.time order
  expect_equal(nrow(res$inst), 3 * 96)
  expect_named(res$inst, c('solar.time','DO.obs.down','DO.mod.down'))
  expect_false(is.unsorted(res$inst$solar.time))

  # the shape downstream consumers see must not depend on per-day vs. joint
  joint <- suppressMessages(metab_bayes_2s(specs=sp, data=dat))
  expect_true(all(names(joint@fit$daily) %in% names(res$daily)))
  expect_named(res$inst, names(joint@fit$inst))
  expect_equal(nrow(res$daily), nrow(joint@fit$daily))

  # run-level warnings/errors stay empty: predict_metab() blanks out every
  # date's estimates when they aren't
  expect_length(res$errors, 0)
  expect_length(res$warnings, 0)
})

# Stand-in for runstan_bayes() that records each Stan data list it is handed
# and returns a failed run, so the per-day path can be followed without Stan
mock_runstan_recorder <- function(env) {
  env$calls <- list()
  function(data_list, ...) {
    env$calls[[length(env$calls) + 1]] <- data_list
    list(compile_time=system.time({}), mcmcfit=NULL, daily=NULL, log=NULL, compile_log=NULL)
  }
}

test_that("metab_bayes_2s(split_dates=TRUE) finds and names each date's NA columns, and only those dates skip Stan", {
  skip_if_not_installed('rstan')
  rec <- new.env()
  testthat::local_mocked_bindings(runstan_bayes=mock_runstan_recorder(rec))

  fx <- corrupt_middle_day(col='DO.obs.up', value=NA_real_)
  sp <- revise(fast_2station_specs(), split_dates=TRUE, day_tests=c())
  mm <- suppressWarnings(suppressMessages(metab_bayes_2s(specs=sp, data=fx$data)))
  daily <- get_fit(mm)$daily

  # metab_bayes_2s() itself found the NA date and named its column, with no
  # hint supplied by the caller
  expect_equal(daily$date, fx$dates)
  bad <- daily[daily$date == fx$bad_date, ]
  expect_match(bad$errors, '^NAs in DO.obs.up; ')
  expect_match(bad$errors, "include 'complete_data' in specs\\$day_tests")
  expect_false(grepl('light|depth', bad$errors))
  expect_true(all(daily$errors[daily$date != fx$bad_date] == ''))

  # and only the other dates reached Stan, each with its own rows
  expect_length(rec$calls, 2)
  good_dates <- setdiff(as.character(fx$dates), as.character(fx$bad_date))
  for(i in seq_along(good_dates)) {
    expected <- prepdata_bayes_2s(fx$data[fx$data$date == as.Date(good_dates[i]), ], specs=sp)
    expect_identical(rec$calls[[i]]$DO_obs_up, expected$DO_obs_up)
  }
})

test_that("metab_bayes_2s(split_dates=TRUE) names each NA date's own columns and still fits the good dates", {
  skip_on_cran()
  skip_if_not_installed('rstan')

  # two NA dates with different NA columns, between and after good dates
  dat <- suppressMessages(mm_align_data_2s(subset_2station_days(two_station_example, 4)))
  dates <- unique(dat$date)
  dat$light[which(dat$date == dates[2])[10]] <- NA
  dat$depth[which(dat$date == dates[3])[20]] <- NA
  sp <- revise(fast_2station_specs(), split_dates=TRUE, day_tests=c())

  mm <- suppressMessages(metab_bayes_2s(specs=sp, data=dat))
  daily <- get_fit(mm)$daily
  expect_equal(daily$date, dates)

  # each bad date names only its own column, with the complete_data hint
  expect_match(daily$errors[2], "^NAs in light; ")
  expect_false(grepl("depth", daily$errors[2]))
  expect_match(daily$errors[3], "^NAs in depth; ")
  expect_false(grepl("light", daily$errors[3]))
  for(i in 2:3) expect_match(daily$errors[i], "include 'complete_data' in specs\\$day_tests")
  expect_true(all(is.na(daily$GPP_daily_50pct[2:3])))

  # the good dates fit normally
  expect_true(all(daily$errors[c(1, 4)] == ''))
  expect_true(all(is.finite(daily$GPP_daily_50pct[c(1, 4)])))
  # the whole run completed, and only the good dates contribute
  # instantaneous predictions
  inst <- get_fit(mm)$inst
  expect_equal(nrow(inst), 2 * 96)
  expect_setequal(unique(mm_date_2s(inst$solar.time)), dates[c(1, 4)])
})

test_that("bayes_perday_2s() reports a date named in na_errors as failed without fitting it", {
  skip_if_not_installed('rstan')
  rec <- new.env()
  testthat::local_mocked_bindings(runstan_bayes=mock_runstan_recorder(rec))

  dat <- suppressMessages(mm_align_data_2s(subset_2station_days(two_station_example, 3)))
  dates <- unique(dat$date)
  na_errors <- setNames('NAs in light; test message', as.character(dates[2]))

  res <- suppressMessages(bayes_perday_2s(dat, specs=fast_2station_specs(), na_errors=na_errors))

  expect_equal(res$daily$date, dates)
  expect_equal(res$daily$errors[2], 'NAs in light; test message')
  expect_true(res$daily$valid_day[2])
  expect_true(is.na(res$daily$GPP_daily_50pct[2]))
  expect_true(dates[2] %in% res$dates_failed)

  # the named date never reached Stan; the other two did
  expect_length(rec$calls, 2)
  expect_false(any(vapply(rec$calls, function(dl) identical(
    dl$DO_obs_down, prepdata_bayes_2s(dat[dat$date == dates[2], ])$DO_obs_down), logical(1))))
})

test_that("a failed Stan run on one day is recorded as that day's warning, leaving other days fit", {
  skip_on_cran()
  skip_if_not_installed('rstan')

  # mirrors the joint fit's treatment of a mode==2L stanfit (a warning, not
  # an error), but scoped to the one day that failed
  sp <- fast_2station_specs()
  dat <- suppressMessages(mm_align_data_2s(subset_2station_days(two_station_example, 3)))
  bad_date <- unique(dat$date)[2]

  setClass('fake_failed_stanfit', representation(mode='integer'))
  fake_stanfit <- methods::new('fake_failed_stanfit', mode=2L)
  real_sampling <- rstan::sampling
  call_n <- 0
  testthat::local_mocked_bindings(
    sampling=function(...) {
      call_n <<- call_n + 1
      if(call_n == 2) fake_stanfit else real_sampling(...)
    }, .package='rstan')

  res <- suppressMessages(bayes_perday_2s(dat, specs=sp))

  expect_equal(res$dates_failed, bad_date)
  expect_length(res$dates_fit, 2)
  bad <- res$daily[res$daily$date == bad_date, ]
  expect_true(is.na(bad$GPP_daily_50pct))
  expect_match(bad$warnings, 'fake_failed_stanfit')
  expect_true(all(!is.na(res$daily$GPP_daily_50pct[res$daily$date != bad_date])))
  expect_equal(nrow(res$inst), 2 * 96)
})

test_that("bayes_perday_2s() compiles the Stan model once and hits the cache for every later day", {
  skip_on_cran()
  skip_if_not_installed('rstan')

  sp <- fast_2station_specs()
  sp$model_path <- mm_locate_filename(sp$model_name)
  stanrds_path <- gsub('\\.stan$', '.stanrds', sp$model_path)

  # warm the cache first, so this test measures the loop's behavior rather
  # than whichever test happened to run first
  invisible(mm_compile_stan_model(sp$model_path))
  expect_true(file.exists(stanrds_path))
  mtime_before <- file.info(stanrds_path)$mtime

  dat <- suppressMessages(mm_align_data_2s(subset_2station_days(two_station_example, 3)))
  res <- suppressMessages(bayes_perday_2s(dat, specs=sp))

  # no day recompiled: the cache file is untouched and no day was charged
  # any compile time
  expect_identical(file.info(stanrds_path)$mtime, mtime_before)
  expect_equal(unname(res$compile_time[['elapsed']]), 0)
})


# split_dates=TRUE wiring ------------------------------------------------

test_that("specs() still defaults split_dates to FALSE but accepts TRUE for two-station", {
  expect_false(specs(mm_name('bayes_2s'))$split_dates)
  # setting it in specs() warns (revise() is the preferred route) exactly as
  # it does for one-station, but the value is honored
  expect_true(suppressWarnings(specs(mm_name('bayes_2s'), split_dates=TRUE))$split_dates)
  expect_true(revise(specs(mm_name('bayes_2s')), split_dates=TRUE)$split_dates)
})

test_that("bayes_1fit_2s() formats with nosplit even when specs$split_dates is TRUE", {
  # split_dates selects runstan_bayes()'s output formatter, not its looping.
  # format_mcmc_mat_split() would collapse the summary to a single row and
  # strip a '[1]' suffix, mangling metab[i,t] and leaving no 'daily' element
  # -- which would silently turn every per-day fit into an apparent failure.
  skip_on_cran()
  skip_if_not_installed('rstan')

  dat <- suppressMessages(mm_align_data_2s(subset_2station_days(two_station_example, 2)))
  sp <- fast_2station_specs()
  sp$split_dates <- TRUE
  sp$model_path <- mm_locate_filename(sp$model_name)

  fit1 <- suppressMessages(bayes_1fit_2s(dat, specs=sp))

  expect_false(is.null(fit1$daily))
  expect_equal(nrow(fit1$daily), 2)
  expect_true(all(c('date','GPP_daily_50pct','ER_daily_50pct','K600_daily_50pct') %in% names(fit1$daily)))
  expect_equal(nrow(fit1$inst), 2 * 96)
})

test_that("metab_bayes_2s(split_dates=TRUE) fits per date and returns the joint mode's output shape", {
  skip_on_cran()
  skip_if_not_installed('rstan')

  dat <- suppressMessages(mm_align_data_2s(subset_2station_days(two_station_example, 3)))
  dates <- unique(dat$date)
  sp <- fast_2station_specs()

  split <- suppressMessages(metab_bayes_2s(specs=revise(sp, split_dates=TRUE), data=dat))
  joint <- suppressMessages(metab_bayes_2s(specs=revise(sp, split_dates=FALSE), data=dat))

  expect_s4_class(split, 'metab_bayes_2s')
  expect_equal(get_fit(split)$daily$date, dates)
  expect_equal(nrow(predict_metab(split)), 3)
  expect_equal(nrow(predict_DO(split)), 3 * 96)
  expect_named(predict_DO(split), c('solar.time','DO.obs.down','DO.mod.down'))
  expect_equal(nrow(get_params(split)), 3)

  # downstream consumers must see the same shape either way
  expect_named(predict_metab(split), names(predict_metab(joint)))
  expect_named(get_params(split), names(get_params(joint)))
  expect_true(all(names(get_fit(joint)$daily) %in% names(get_fit(split)$daily)))

  # per-day estimates are plausible and independent, not copies of each other
  pm <- predict_metab(split)
  expect_true(all(pm$GPP > 0)); expect_true(all(pm$ER < 0)); expect_true(all(pm$K600 > 0))
  expect_false(anyDuplicated(pm$GPP) > 0)
})

test_that("metab_bayes_2s(split_dates=TRUE) reports every date as failed when none reaches Stan", {
  skip_on_cran()
  skip_if_not_installed('rstan')

  dat <- subset_2station_days(two_station_example, 3)
  alignment <- suppressMessages(mm_align_2s(v(dat)))
  dates <- unique(alignment$date)
  # a contiguous 4-hour NA run in each day's upstream DO, left unfilled, and
  # day_tests=c() lets it through to the fit
  for(dt in as.list(dates)) {
    up <- alignment$shift_idx[alignment$date == dt][40] + 0:15
    dat$DO.obs.up[up] <- u(NA_real_, get_units(dat$DO.obs.up))
  }
  sp <- revise(fast_2station_specs(), split_dates=TRUE, day_tests=c())

  # no date produced predictions, so the run-level "Modeling failed" warning
  # fires, as it does for a failed joint fit
  expect_warning(
    mm <- suppressMessages(metab_bayes_2s(specs=sp, data=suppressMessages(mm_align_data_2s(dat)))),
    'Modeling failed')

  expect_s4_class(mm, 'metab_bayes_2s')
  daily <- get_fit(mm)$daily
  expect_equal(daily$date, dates)
  expect_true(all(is.na(daily$GPP_daily_50pct)))
  expect_true(all(grepl('NAs in DO.obs.up', daily$errors)))
  expect_null(get_fit(mm)$inst)
})

test_that("the mcmc slot holds one stanfit jointly and a date-named list per day, as for one-station", {
  skip_on_cran()
  skip_if_not_installed('rstan')

  dat <- suppressMessages(mm_align_data_2s(subset_2station_days(two_station_example, 3)))
  dates <- unique(dat$date)
  sp <- fast_2station_specs()

  joint <- suppressMessages(metab_bayes_2s(specs=revise(sp, split_dates=FALSE, keep_mcmcs=TRUE), data=dat))
  expect_s4_class(get_mcmc(joint), 'stanfit')

  split <- suppressMessages(metab_bayes_2s(specs=revise(sp, split_dates=TRUE, keep_mcmcs=TRUE), data=dat))
  mc <- get_mcmc(split)
  expect_type(mc, 'list')
  expect_length(mc, 3)
  expect_equal(names(mc), as.character(dates))
  expect_true(all(sapply(mc, function(x) inherits(x, 'stanfit'))))

  # keep_mcmcs=FALSE keeps nothing, and says so the same way either mode does
  none <- suppressMessages(metab_bayes_2s(specs=revise(sp, split_dates=TRUE, keep_mcmcs=FALSE), data=dat))
  expect_null(get_mcmc(none))
  expect_equal(nrow(predict_metab(none)), 3)
})

test_that("keep_mcmcs/keep_mcmc_data accept a vector of dates in per-day mode only", {
  skip_on_cran()
  skip_if_not_installed('rstan')

  dat <- suppressMessages(mm_align_data_2s(subset_2station_days(two_station_example, 3)))
  dates <- unique(dat$date)
  sp <- fast_2station_specs()
  pick <- dates[c(1, 3)]

  mm <- suppressMessages(metab_bayes_2s(
    specs=revise(sp, split_dates=TRUE, keep_mcmcs=pick, keep_mcmc_data=dates[2]), data=dat))

  mc <- get_mcmc(mm)
  expect_equal(names(mc), as.character(dates))
  expect_equal(names(mc)[!sapply(mc, is.null)], as.character(pick))

  md <- get_mcmc_data(mm)
  expect_equal(names(md)[!sapply(md, is.null)], as.character(dates[2]))
  kept <- md[[as.character(dates[2])]]
  expect_true(all(c('n_obs','n_days','DO_obs_up') %in% names(kept)))
  expect_equal(kept$n_days, 1L)

  # a date vector is meaningless when there is only one fit
  expect_error(
    metab_bayes_2s(specs=revise(sp, split_dates=FALSE, keep_mcmcs=pick), data=dat),
    'keep_mcmcs must be a single logical value')
  expect_error(
    metab_bayes_2s(specs=revise(sp, split_dates=FALSE, keep_mcmc_data=pick), data=dat),
    'keep_mcmc_data must be a single logical value')
})

#### dropped days are reported, not omitted ####

test_that("mm_align_2s() reports the days it drops, and why", {
  # the travel-time ceiling
  alignment <- suppressMessages(mm_align_2s(make_2day_ceiling_data()))
  expect_equal(nrow(alignment$removed), 1)
  expect_named(alignment$removed, c('date','errors'))
  expect_equal(as.character(alignment$removed$date), "2050-06-02")
  expect_match(alignment$removed$errors, "travel.time exceeds the 0\\.42 days \\(10 hours\\) ceiling \\(0\\.63 days \\(15 hours\\)\\)")

  # a day that doesn't fill its 06:00-06:00 window. Shortening the second day
  # of a two-day fixture, since dropping the only day is an error, not a result
  dat <- make_2day_ceiling_data()
  dat$travel.time <- 0.01 # put day 2 back under the ceiling
  dat <- dat[1:(nrow(dat) - 20), ]
  alignment3 <- suppressMessages(mm_align_2s(dat))
  expect_equal(nrow(alignment3$removed), 1)
  expect_equal(as.character(alignment3$removed$date), "2050-06-02")
  expect_match(alignment3$removed$errors, "does not fill the 06:00-06:00 window \\(268 of 288")

  # nothing dropped, nothing reported -- in particular the fixture's
  # deliberate 3-row lead-in block (2050-05-31) is not reported as a lost day
  alignment4 <- suppressMessages(mm_align_2s(make_2station_data()))
  expect_equal(nrow(alignment4$removed), 0)
  expect_false("2050-05-31" %in% as.character(alignment4$removed$date))
})

test_that("mm_align_2s() reports a mid-record day with no upstream lead-in at all", {
  # A handful of rows stranded inside an upstream outage, as happens mid-record
  # in real data: the lag is 3 timesteps and the bins they reach back to all
  # fall in the gap, so no row has an upstream match and the whole day would
  # otherwise disappear -- no message, no row, nothing
  day1 <- make_2station_data() # 3 lead-in rows + all of 2050-06-01
  stranded <- day1[1:3, ]
  stranded$solar.time <- as.POSIXct("2050-06-02 12:00:00", tz="UTC") +
    as.difftime((0:2) * 5, units="mins")
  dat <- rbind(day1, stranded)

  alignment <- expect_message(
    mm_align_2s(dat),
    "dropping 1 day\\(s\\) with no upstream lead-in at all: 2050-06-02 \\(3 row\\(s\\) supplied\\)")

  # the good day is unaffected
  expect_equal(as.character(unique(alignment$date)), "2050-06-01")
  expect_equal(alignment$n_days, 1)

  # and the stranded day is reported rather than omitted
  expect_equal(as.character(alignment$removed$date), "2050-06-02")
  expect_match(alignment$removed$errors, "no upstream observation at any row's travel-time offset")
})

test_that("a clean fit reports no invalid days", {
  skip_on_cran()
  skip_if_not_installed('rstan')

  sp <- fast_2station_specs()
  dat <- suppressMessages(mm_align_data_2s(subset_2station_days(two_station_example, 3)))

  mm <- suppressMessages(metab_bayes_2s(specs=sp, data=dat))

  expect_equal(nrow(mm@fit$daily), 3)
  expect_true(all(mm@fit$daily$valid_day))
  expect_true(all(mm@fit$daily$errors == ''))
  expect_true(all(is.finite(predict_metab(mm)$GPP)))
})

test_that("a day dropped by day_tests comes back as a valid_day=FALSE row (joint fit)", {
  skip_on_cran()
  skip_if_not_installed('rstan')

  fx <- corrupt_middle_day()
  sp <- fast_2station_specs()

  mm <- expect_message(
    metab_bayes_2s(specs=sp, data=fx$data),
    "dropping 1 day\\(s\\) that fail day_tests")

  daily <- mm@fit$daily

  # every date supplied is accounted for, in date order
  expect_equal(nrow(daily), 3)
  expect_equal(daily$date, fx$dates)
  expect_false(is.unsorted(daily$date))

  bad <- daily[daily$date == fx$bad_date, ]
  good <- daily[daily$date != fx$bad_date, ]
  expect_false(bad$valid_day)
  expect_true(all(good$valid_day))
  expect_equal(bad$errors, 'depth <= 0')
  expect_true(all(good$errors == ''))
  expect_true(is.na(bad$GPP_daily_50pct))

  # the run as a whole succeeded, with no run-level error to blank out the
  # two good days' estimates
  expect_length(mm@fit$errors, 0)
  expect_true(all(!is.na(good$GPP_daily_50pct)))

  # predict_metab() shows the dropped date too, with NA estimates
  pm <- predict_metab(mm)
  expect_equal(nrow(pm), 3)
  expect_true(is.na(pm$GPP[pm$date == fx$bad_date]))
  expect_true(all(is.finite(pm$GPP[pm$date != fx$bad_date])))

  # instantaneous output covers only the modeled days
  expect_equal(nrow(predict_DO(mm)), 2 * 96)
  expect_false(fx$bad_date %in% mm_date_2s(predict_DO(mm)$solar.time))
})

test_that("a day dropped by day_tests comes back as a valid_day=FALSE row (per-day fit)", {
  skip_on_cran()
  skip_if_not_installed('rstan')

  fx <- corrupt_middle_day()
  sp <- fast_2station_specs()
  sp$split_dates <- TRUE

  mm <- suppressMessages(metab_bayes_2s(specs=sp, data=fx$data))
  daily <- mm@fit$daily

  expect_equal(nrow(daily), 3)
  expect_equal(daily$date, fx$dates)

  bad <- daily[daily$date == fx$bad_date, ]
  expect_false(bad$valid_day)
  expect_equal(bad$errors, 'depth <= 0')
  expect_true(is.na(bad$GPP_daily_50pct))
  expect_true(all(daily$valid_day[daily$date != fx$bad_date]))
  expect_true(all(!is.na(daily$GPP_daily_50pct[daily$date != fx$bad_date])))

  # reported identically in both modes: valid_day=FALSE, not a failed fit
  expect_equal(nrow(predict_DO(mm)), 2 * 96)
})

test_that("a dataset with no usable days errors clearly rather than reaching Stan", {
  dat <- suppressMessages(mm_align_data_2s(subset_2station_days(two_station_example, 2)))
  for(dt in as.list(unique(dat$date))) {
    dat$depth[which(dat$date == dt)[1]] <- 0
  }

  expect_error(
    suppressMessages(metab_bayes_2s(specs=fast_2station_specs(), data=dat)),
    "no days remain after day_tests")
})

test_that("specs$day_tests is honored, not hard-coded", {
  dat <- corrupt_middle_day()$data
  sp <- fast_2station_specs()
  sp$day_tests <- 'complete_data' # drop pos_depth

  # the depth<=0 day is no longer rejected before fitting
  res <- suppressMessages(mm_filter_valid_days_2s(dat, day_tests=sp$day_tests))
  expect_equal(nrow(res$removed), 0)
})


# alignment record carried into the fit -------------------------------------

test_that("a fit on data with unknown alignment drops says so, and reports only day_tests drops", {
  skip_if_not_installed('rstan')
  rec <- new.env()
  testthat::local_mocked_bindings(runstan_bayes=mock_runstan_recorder(rec))

  # a hand-marked frame has no record of which days its alignment dropped
  fx <- corrupt_middle_day()
  unknown <- mm_as_aligned_2s(as.data.frame(fx$data))
  expect_null(attr(unknown, 'removed'))

  expect_message(
    mm <- suppressWarnings(metab_bayes_2s(specs=fast_2station_specs(), data=unknown)),
    "days dropped during alignment are unknown")
  daily <- get_fit(mm)$daily
  expect_equal(daily$date, fx$dates)
  expect_equal(daily$valid_day, c(TRUE, FALSE, TRUE))
  expect_equal(daily$errors[2], 'depth <= 0')

  # with the record present, its dropped days come back too
  ceiling_aligned <- suppressMessages(mm_align_data_2s(make_2day_ceiling_data()))
  mm2 <- suppressWarnings(suppressMessages(metab_bayes_2s(specs=fast_2station_specs(), data=ceiling_aligned)))
  expect_equal(as.character(get_fit(mm2)$daily$date), c("2050-06-01", "2050-06-02"))
  expect_equal(get_fit(mm2)$daily$valid_day, c(TRUE, FALSE))
})

test_that("a fit reports days dropped by the travel-time ceiling and by day_tests together, in date order", {
  skip_if_not_installed('rstan')
  rec <- new.env()
  testthat::local_mocked_bindings(runstan_bayes=mock_runstan_recorder(rec))

  # three 06:00-06:00 days: day 1 fails day_tests (depth 0), day 2 is clean,
  # day 3's travel.time exceeds the default ceiling
  dat <- make_2station_data(n=3 + 3*288)
  dates <- as.Date(c("2050-06-01", "2050-06-02", "2050-06-03"))
  day_of <- mm_date_2s(dat$solar.time)
  dat$depth[which(day_of == dates[1])[10]] <- 0
  dat$travel.time[day_of == dates[3]] <- 15/24

  aligned <- suppressMessages(mm_align_data_2s(dat))
  expect_equal(attr(aligned, 'removed')$date, dates[3])

  mm <- suppressWarnings(suppressMessages(metab_bayes_2s(specs=fast_2station_specs(), data=aligned)))
  daily <- get_fit(mm)$daily
  expect_equal(daily$date, dates)
  expect_equal(daily$valid_day, c(FALSE, TRUE, FALSE))
  expect_equal(daily$errors[1], 'depth <= 0')
  expect_match(daily$errors[3], "travel.time exceeds the 0\\.42 days \\(10 hours\\) ceiling")

  # only the clean day reached Stan
  expect_length(rec$calls, 1)
  expect_equal(rec$calls[[1]]$n_days, 1)
})

test_that("the fit records the travel-time ceiling the data were aligned with", {
  skip_if_not_installed('rstan')
  rec <- new.env()
  testthat::local_mocked_bindings(runstan_bayes=mock_runstan_recorder(rec))
  dat <- make_2station_data()

  mm <- suppressWarnings(suppressMessages(metab_bayes_2s(
    specs=fast_2station_specs(), data=suppressMessages(mm_align_data_2s(dat, max_travel_time_days=0.3)))))
  expect_equal(get_specs(mm)$aligned_max_travel_time_days, 0.3)

  # unknown for a hand-marked frame: present, and NULL
  hand <- mm_as_aligned_2s(as.data.frame(suppressMessages(mm_align_data_2s(dat))))
  mm2 <- suppressWarnings(suppressMessages(metab_bayes_2s(specs=fast_2station_specs(), data=hand)))
  expect_true('aligned_max_travel_time_days' %in% names(get_specs(mm2)))
  expect_null(get_specs(mm2)$aligned_max_travel_time_days)

  # taken from the data every time, not carried over in reused specs
  mm3 <- suppressWarnings(suppressMessages(metab_bayes_2s(specs=get_specs(mm), data=hand)))
  expect_null(get_specs(mm3)$aligned_max_travel_time_days)

  # no gap-filling setting is recorded
  expect_false('aligned_max_gap_hours' %in% names(get_specs(mm)))
  expect_false('aligned_max_gap_hours' %in% names(get_specs(mm2)))
})

test_that("aligned data missing its date column fails validation clearly", {
  aligned <- suppressMessages(mm_align_data_2s(make_2station_data()))
  no_date <- aligned[, setdiff(names(aligned), 'date')]
  expect_s3_class(no_date, 'aligned_2s')
  expect_error(metab_bayes_2s(data=no_date), "aligned data is missing these columns: date")
})

test_that("metab_bayes_2s() rejects aligned data whose solar.time is not UTC", {
  aligned <- suppressMessages(mm_align_data_2s(make_2station_data()))
  attr(aligned$solar.time, 'tzone') <- 'America/Denver'
  expect_error(mm_validate_data_2station(aligned), "timezone 'UTC'")
  expect_error(metab_bayes_2s(data=aligned), "timezone 'UTC'")
})

test_that("specs() no longer takes the travel-time ceiling or the gap-filling tolerance", {
  sp <- specs(mm_name('bayes_2s'))
  expect_false(any(c('max_travel_time_days', 'max_gap_hours') %in% names(sp)))
  expect_error(specs(mm_name('bayes_2s'), max_travel_time_days=0.3), "unused argument")
  expect_error(specs(mm_name('bayes_2s'), max_gap_hours=2), "unused argument")
})
