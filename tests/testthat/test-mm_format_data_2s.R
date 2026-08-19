
# Small, hourly, traceable upstream/downstream/light data.frames -- easier to
# reason about exact drop counts and NA positions than the real (15-min,
# gappy) two_station_raw_example.
make_hourly <- function(start, n, ...) {
  data.frame(
    timestamp = as.POSIXct(start, tz="UTC") + as.difftime(seq_len(n) - 1, units="hours"),
    ...
  )
}

test_that("mm_format_data_2s errors clearly on missing required columns", {
  upstream <- make_hourly("2050-06-01 00:00:00", 3, DO.obs=1:3, DO.sat=4:6)
  downstream <- make_hourly("2050-06-01 00:00:00", 3, DO.obs=1:3, DO.sat=4:6,
                             temp.water=20, depth=0.5, travel.time=0.01)
  light <- make_hourly("2050-06-01 00:00:00", 3, light=100)

  expect_error(
    mm_format_data_2s(dplyr::select(upstream, -DO.sat), downstream, light),
    "upstream is missing these columns: DO.sat")
  expect_error(
    mm_format_data_2s(upstream, dplyr::select(downstream, -travel.time), light),
    "downstream is missing these columns: travel.time")
  expect_error(
    mm_format_data_2s(upstream, downstream, dplyr::select(light, -light)),
    "light is missing these columns: light")
})

test_that("a non-numeric column, including character DO.sat with literal \"NA\" strings, is a clear error rather than a silent fix", {
  upstream <- make_hourly("2050-06-01 00:00:00", 3, DO.obs=c(9.1, 9.2, 9.3),
                           DO.sat=c("9.9", "NA", "10.1"))
  downstream <- make_hourly("2050-06-01 00:00:00", 3, DO.obs=c(8.1, 8.2, 8.3), DO.sat=c(8.9, 9.0, 9.1),
                             temp.water=20, depth=0.5, travel.time=0.01)
  light <- make_hourly("2050-06-01 00:00:00", 3, light=c(100, 200, 300))

  expect_error(
    mm_format_data_2s(upstream, downstream, light),
    "upstream\\$DO.sat must be numeric, but is character")
})

test_that("a non-POSIXct timestamp, or a POSIXct timestamp with the wrong timezone, is a clear error", {
  downstream <- make_hourly("2050-06-01 00:00:00", 3, DO.obs=1:3, DO.sat=1:3 + 10,
                             temp.water=20, depth=0.5, travel.time=0.01)
  light <- make_hourly("2050-06-01 00:00:00", 3, light=1:3)

  upstream_char_ts <- make_hourly("2050-06-01 00:00:00", 3, DO.obs=1:3, DO.sat=1:3 + 10)
  upstream_char_ts$timestamp <- as.character(upstream_char_ts$timestamp)
  expect_error(
    mm_format_data_2s(upstream_char_ts, downstream, light),
    "upstream\\$timestamp must be of class POSIXct")

  upstream_wrong_tz <- make_hourly("2050-06-01 00:00:00", 3, DO.obs=1:3, DO.sat=1:3 + 10)
  upstream_wrong_tz$timestamp <- lubridate::with_tz(upstream_wrong_tz$timestamp, "America/Denver")
  expect_error(
    mm_format_data_2s(upstream_wrong_tz, downstream, light),
    "upstream\\$timestamp must have timezone 'UTC'")
})

test_that("upstream/downstream deployment windows are trimmed to their overlap, with one message naming the counts", {
  # downstream: 00:00-09:00 (10 rows). upstream & light both extend 2 hours
  # before and 2 hours after that window (14 rows each); only the 4
  # out-of-window rows in each should be dropped, and downstream (already
  # inside the overlap) should lose none.
  downstream <- make_hourly("2050-06-01 00:00:00", 10, DO.obs=1:10, DO.sat=1:10 + 10,
                             temp.water=20, depth=0.5, travel.time=0.01)
  upstream <- make_hourly("2050-05-31 22:00:00", 14, DO.obs=1:14, DO.sat=1:14 + 10)
  light <- make_hourly("2050-05-31 22:00:00", 14, light=1:14 * 10)

  expect_message(
    out <- mm_format_data_2s(upstream, downstream, light),
    "dropping 4 upstream row\\(s\\), 0 downstream row\\(s\\), 4 light row\\(s\\)")

  expect_equal(nrow(out), 10)
  # every downstream row still has a real upstream match after trimming (the
  # trimmed-in upstream rows fully cover downstream's window); light itself
  # is separately NA here because these 10 hourly rows don't fill either of
  # the two partial 06:00-06:00 days they span -- see mm_lag_light_2s()'s
  # day-completeness rule, exercised elsewhere, not what this test is about
  expect_false(any(is.na(v(out$DO.obs.up))))
})

test_that("no trimming and no dropped rows when inputs already share one deployment window (two_station_raw_example)", {
  data("two_station_raw_example", envir=environment())
  up <- two_station_raw_example$upstream
  down <- two_station_raw_example$downstream
  lt <- two_station_raw_example$light

  # gap filling reports what it bridged, so this call is no longer silent;
  # what this test is about is that nothing gets *trimmed* -- the three
  # inputs already share one deployment window
  msgs <- capture_messages(out <- mm_format_data_2s(up, down, lt))
  expect_false(any(grepl('trimming to the shared', msgs)))

  # no downstream row is dropped. Rows are now *added*, though: filling
  # bridges this record's short gaps (see test-mm_fill_gaps_2s.R), so the
  # output is a superset of downstream's timestamps rather than a match
  expect_true(all(v(mm_snap_to_bin_2s(down$timestamp)) %in% v(out$solar.time)))
  expect_gte(nrow(out), nrow(down))
})

test_that("upstream and downstream with no overlapping window is an error", {
  downstream <- make_hourly("2050-06-01 00:00:00", 3, DO.obs=1:3, DO.sat=1:3 + 10,
                             temp.water=20, depth=0.5, travel.time=0.01)
  upstream <- make_hourly("2050-06-05 00:00:00", 3, DO.obs=1:3, DO.sat=1:3 + 10)
  light <- make_hourly("2050-06-01 00:00:00", 3, light=1:3)

  expect_error(mm_format_data_2s(upstream, downstream, light), "no overlapping deployment window")
})

test_that("a downstream timestamp with no matching upstream/light bin is interpolated when short, NA when too long", {
  # upstream/light are missing bins that downstream has -- a real mid-series
  # gap. At an hourly timestep the 1-hour tolerance spans exactly one bin,
  # so one missing bin is bridged and two consecutive ones are not
  downstream <- make_hourly("2050-06-01 00:00:00", 8, DO.obs=1:8, DO.sat=1:8 + 10,
                             temp.water=20, depth=0.5, travel.time=0.01)
  drop_bins <- function(df, drop) df[!(v(df$timestamp) %in% drop), ]
  make_up <- function(drop) drop_bins(
    make_hourly("2050-06-01 00:00:00", 8, DO.obs=1:8, DO.sat=1:8 + 10), drop)
  make_lt <- function(drop) drop_bins(
    make_hourly("2050-06-01 00:00:00", 8, light=1:8 * 10), drop)

  one <- as.POSIXct("2050-06-01 02:00:00", tz="UTC")
  out <- suppressMessages(mm_format_data_2s(make_up(one), downstream, make_lt(one)))

  # the row is still there and no longer NA: it is interpolated from its
  # neighbours (DO.obs.up 2 and 4 bracket the missing 3)
  expect_equal(nrow(out), 8)
  gap_row <- which(v(out$solar.time) == one)
  expect_false(any(is.na(v(out$DO.obs.up))))
  expect_equal(v(out$DO.obs.up)[gap_row], 3)

  # two consecutive missing bins exceed the tolerance, so the original
  # contract still holds there: an unmatched bin is an NA, never a lost row
  two <- one + as.difftime(c(0, 1), units="hours")
  out2 <- suppressMessages(mm_format_data_2s(make_up(two), downstream, make_lt(two)))

  expect_equal(nrow(out2), 8)
  gap_rows <- which(v(out2$solar.time) %in% two)
  expect_true(all(is.na(v(out2$DO.obs.up)[gap_rows])))
  expect_true(all(is.na(v(out2$DO.sat.up)[gap_rows])))
  expect_false(any(is.na(v(out2$DO.obs.up)[-gap_rows])))
})

test_that("the timestep is taken from downstream alone, so a coarser upstream/light is binned onto downstream's grid", {
  # downstream is hourly, upstream/light are 2-hourly: the three series have
  # different modal timesteps. Only downstream's is used, for all three, so
  # upstream's rows land on every other hourly bin. Were the timestep instead
  # re-detected per series, upstream's 2-hour bins would be numbered half as
  # far apart and would match downstream's hourly bins at the wrong rows.
  make_2hourly <- function(start, n, ...) {
    data.frame(
      timestamp = as.POSIXct(start, tz="UTC") + as.difftime((seq_len(n) - 1) * 2, units="hours"),
      ...
    )
  }
  downstream <- make_hourly("2050-06-01 00:00:00", 9, DO.obs=1:9, DO.sat=1:9 + 10,
                             temp.water=20, depth=0.5, travel.time=0.01)
  upstream <- make_2hourly("2050-06-01 00:00:00", 5, DO.obs=c(100,102,104,106,108), DO.sat=200)
  light <- make_2hourly("2050-06-01 00:00:00", 5, light=c(10,12,14,16,18))

  out <- suppressMessages(mm_format_data_2s(upstream, downstream, light))

  # downstream's row structure is preserved exactly
  expect_equal(nrow(out), 9)
  expect_equal(v(out$solar.time), downstream$timestamp)

  # upstream landed on bins 0,2,4,6,8 of downstream's hourly grid; the odd
  # bins in between are single-bin gaps, bridged by interpolation
  expect_equal(v(out$DO.obs.up), c(100,101,102,103,104,105,106,107,108))

  # (light is NA throughout here only because mm_lag_light_2s() reports a
  # within-day proportion and this 9-row fixture spans no 06:00-06:00 day)
})

test_that("a duplicate timestep bin in downstream is a clear error", {
  # 6 regularly-spaced rows (enough for mm_get_timestep()'s modal-timestep
  # detection to work cleanly) plus one extra row duplicating an existing bin
  downstream <- make_hourly("2050-06-01 00:00:00", 6, DO.obs=1:6, DO.sat=1:6 + 10,
                             temp.water=20, depth=0.5, travel.time=0.01)
  downstream <- rbind(downstream, downstream[3, ])
  upstream <- make_hourly("2050-06-01 00:00:00", 6, DO.obs=1:6, DO.sat=1:6 + 10)
  light <- make_hourly("2050-06-01 00:00:00", 6, light=1:6)

  expect_error(
    mm_format_data_2s(upstream, downstream, light),
    "downstream has more than one row in the same timestep bin")
})

test_that("mm_format_data_2s()'s output passes mm_validate_data()/mm_validate_data_2station() and runs end-to-end through mm_align_2s()/prepdata_bayes_2s() on the real, gappy two_station_raw_example", {
  data("two_station_raw_example", envir=environment())
  out <- mm_format_data_2s(
    two_station_raw_example$upstream, two_station_raw_example$downstream, two_station_raw_example$light)

  dat_list <- mm_validate_data(out, metab_class="metab_bayes_2s")
  expect_no_error(mm_validate_data_2station(dat_list$data))

  # the same route metab_bayes_2s() takes: align, then drop days whose modeled
  # values fail day_tests, then prep. Skipping the middle step leaves NAs in
  # the matrices -- 86 of this record's 891 aligned days carry one, mostly
  # from upstream outages -- and Stan rejects those outright
  sp <- specs(mm_name('bayes_2s'))
  aln <- suppressMessages(mm_align_2s(v(dat_list$data), max_travel_time_hours=sp$max_travel_time_hours))
  filtered <- suppressMessages(mm_filter_valid_days_2s(dat_list$data, aln, day_tests=sp$day_tests))
  expect_gt(nrow(filtered$removed), 0) # this record does exercise the filter

  prep <- suppressMessages(prepdata_bayes_2s(dat_list$data, specs=sp, aln=filtered$aln))
  expect_equal(prep$n_obs, 96) # 15-min timestep -> 96 obs per 06:00-06:00 day
  expect_gt(prep$n_days, 0)
  expect_equal(prep$n_days, filtered$aln$n_days)
  for(varname in c('DO_obs_up','DO_sat_up','DO_obs_down','DO_sat_down','light','depth','temp_water','travel_time')) {
    expect_equal(dim(prep[[varname]]), c(prep$n_obs, prep$n_days), info=varname)
    expect_false(anyNA(prep[[varname]]), info=varname)
  }
})

test_that("a user's own already-formatted data (not built via mm_format_data_2s) passes validation independently", {
  # matches two_station_example's shape directly, with no provenance marker
  dat <- data.frame(
    solar.time = as.POSIXct("2050-06-01 00:00:00", tz="UTC") + as.difftime(0:9 * 5, units="mins"),
    DO.obs.up = rep(9, 10), DO.sat.up = rep(10, 10),
    DO.obs.down = rep(8.8, 10), DO.sat.down = rep(9.9, 10),
    light = rep(0.1, 10), depth = rep(0.5, 10), temp.water = rep(20, 10),
    travel.time = rep(0.01, 10))

  dat_list <- mm_validate_data(dat, metab_class="metab_bayes_2s")
  expect_no_error(mm_validate_data_2station(dat_list$data))
})
