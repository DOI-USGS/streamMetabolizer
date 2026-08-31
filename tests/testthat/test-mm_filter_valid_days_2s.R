# Two 06:00-06:00 days at a 5-minute timestep, preceded by 3 lead-in rows.
# travel.time=0.01 d is 2.88 timesteps, so lag = 3: each modeled row draws its
# upstream values from 3 rows earlier. Rows 1-3 are lead-in (no upstream match
# of their own), rows 4-291 are 2050-06-01, rows 292-579 are 2050-06-02. Day
# 2's first three upstream values therefore come from rows 289-291, which sit
# inside day 1 -- the overlap the modeled-frame indexing exists to get right.
make_2day_2station_data <- function() {
  n <- 3 + 2*288
  data.frame(
    solar.time = as.POSIXct("2050-06-01 05:45:00", tz="UTC") +
      as.difftime((seq_len(n) - 1) * 5, units="mins"),
    DO.obs.up = rep(9, n),
    DO.sat.up = rep(10, n),
    DO.obs.down = rep(8.8, n),
    DO.sat.down = rep(9.9, n),
    light = rep(300, n),
    depth = rep(0.5, n),
    temp.water = rep(20, n),
    travel.time = rep(0.01, n)
  )
}

filter_2day <- function(dat, ...) {
  aln <- suppressMessages(mm_align_2s(v(dat)))
  mm_filter_valid_days_2s(dat, aln, ...)
}

test_that("the two-day fixture is what the rest of this file assumes", {
  dat <- make_2day_2station_data()
  aln <- suppressMessages(mm_align_2s(v(dat)))

  expect_equal(aln$n_days, 2)
  expect_equal(aln$n_obs, 288)
  expect_equal(as.character(unique(aln$date)), c("2050-06-01", "2050-06-02"))
  expect_equal(range(aln$keep), c(4, 579))
  # the lag is 3 rows, so day 2 reaches back into day 1's rows for its first
  # few upstream values
  expect_equal(aln$shift_idx, aln$keep - 3)
  expect_equal(min(aln$shift_idx[aln$date == as.Date("2050-06-02")]), 289)
  expect_true(289 %in% aln$keep[aln$date == as.Date("2050-06-01")])
})

test_that("mm_modeled_rows_2s draws upstream from shift_idx and everything else from keep", {
  dat <- make_2day_2station_data()
  # distinct values per row make the indexing visible
  dat$DO.obs.up <- seq_len(nrow(dat))
  dat$DO.obs.down <- seq_len(nrow(dat)) + 10000
  aln <- suppressMessages(mm_align_2s(v(dat)))

  modeled <- mm_modeled_rows_2s(dat, aln)

  expect_equal(nrow(modeled), length(aln$keep))
  expect_named(modeled, c('solar.time','DO.obs.up','DO.sat.up','DO.obs.down',
                          'DO.sat.down','light','depth','temp.water','travel.time'))
  expect_equal(modeled$DO.obs.up, dat$DO.obs.up[aln$shift_idx])
  expect_equal(modeled$DO.obs.down, dat$DO.obs.down[aln$keep])
  expect_equal(modeled$solar.time, dat$solar.time[aln$keep])
})

test_that("mm_modeled_rows_2s strips units", {
  dat <- make_2day_2station_data()
  template <- mm_data(solar.time, DO.obs.up, DO.sat.up, DO.obs.down, DO.sat.down,
                      light, depth, temp.water, travel.time)
  for(col in names(template)) dat[[col]] <- u(dat[[col]], get_units(template[[col]]))
  dat$light <- u(v(dat$light), NA)
  aln <- suppressMessages(mm_align_2s(v(dat)))

  expect_false(is.unitted(mm_modeled_rows_2s(dat, aln)))
})

test_that("a clean dataset passes untouched", {
  dat <- make_2day_2station_data()
  aln <- suppressMessages(mm_align_2s(v(dat)))

  res <- expect_silent(mm_filter_valid_days_2s(dat, aln))

  expect_equal(nrow(res$removed), 0)
  expect_named(res$removed, c('date','errors'))
  expect_equal(res$aln, aln)
})

test_that("complete_data drops a day with an NA in a downstream column, leaving the other intact", {
  dat <- make_2day_2station_data()
  dat$temp.water[300] <- NA # a modeled row of 2050-06-02

  res <- expect_message(filter_2day(dat), "dropping 1 day\\(s\\) that fail day_tests")

  expect_equal(as.character(res$removed$date), "2050-06-02")
  expect_equal(res$removed$errors, "NAs in temp.water")
  expect_equal(as.character(unique(res$aln$date)), "2050-06-01")
  expect_equal(res$aln$n_days, 1)
  expect_length(res$aln$keep, 288)
  expect_length(res$aln$shift_idx, 288)
})

test_that("an upstream NA is charged to the day that depends on it, not the day it sits in", {
  # This is the whole reason the tests run on the modeled frame. Row 289 is
  # one of 2050-06-01's own rows, but nothing in 2050-06-01 is modeled from
  # it -- 2050-06-02's first upstream value is. Checking each day's own rows
  # would blame the wrong date.
  dat <- make_2day_2station_data()
  dat$DO.obs.up[289] <- NA

  res <- suppressMessages(filter_2day(dat))

  expect_equal(as.character(res$removed$date), "2050-06-02")
  expect_equal(res$removed$errors, "NAs in DO.obs.up")
  expect_equal(as.character(unique(res$aln$date)), "2050-06-01")
})

test_that("an NA in a lead-in row that nothing is modeled from drops nothing", {
  # rows 1-3 supply no modeled row's upstream values (day 1's first modeled
  # row is row 4, whose upstream is row 1) -- row 1 IS used, rows 2-3 are used
  # by rows 5-6. Downstream columns of rows 1-3, though, are never modeled.
  dat <- make_2day_2station_data()
  dat$DO.obs.down[2] <- NA
  dat$depth[2] <- NA

  res <- expect_silent(filter_2day(dat))
  expect_equal(nrow(res$removed), 0)
  expect_equal(res$aln$n_days, 2)
})

test_that("pos_depth drops a day whose depth is not positive", {
  # zero depth: the model divides by it
  dat <- make_2day_2station_data()
  dat$depth[300] <- 0

  res <- suppressMessages(filter_2day(dat))
  expect_equal(as.character(res$removed$date), "2050-06-02")
  expect_equal(res$removed$errors, "depth <= 0")
  expect_equal(as.character(unique(res$aln$date)), "2050-06-01")

  # negative depth, same treatment
  dat2 <- make_2day_2station_data()
  dat2$depth[10] <- -0.2 # a modeled row of 2050-06-01
  res2 <- suppressMessages(filter_2day(dat2))
  expect_equal(as.character(res2$removed$date), "2050-06-01")
  expect_equal(res2$removed$errors, "depth <= 0")

  # and it is genuinely pos_depth doing the work, not complete_data
  res3 <- suppressMessages(filter_2day(dat, day_tests='complete_data'))
  expect_equal(nrow(res3$removed), 0)
})

test_that("a day failing several tests reports every reason", {
  dat <- make_2day_2station_data()
  dat$depth[300] <- 0
  dat$light[301] <- NA

  res <- suppressMessages(filter_2day(dat))
  expect_equal(nrow(res$removed), 1)
  expect_equal(res$removed$errors, "NAs in light; depth <= 0")
})

test_that("a single-day dataset can be dropped entirely", {
  # the n_days == 0 boundary, reached with nothing left to fall back on
  dat <- make_2day_2station_data()[1:291, ] # 3 lead-in rows + 2050-06-01 only
  dat$depth[10] <- 0

  res <- suppressMessages(filter_2day(dat))
  expect_equal(nrow(res$removed), 1)
  expect_equal(res$aln$n_days, 0)
  expect_length(res$aln$keep, 0)
  expect_length(res$aln$date, 0)
})

test_that("both days can be dropped, leaving an empty alignment for the caller to handle", {
  dat <- make_2day_2station_data()
  dat$depth[10] <- 0
  dat$depth[300] <- 0

  res <- suppressMessages(filter_2day(dat))
  expect_equal(nrow(res$removed), 2)
  expect_equal(res$aln$n_days, 0)
  expect_length(res$aln$keep, 0)
})

test_that("an empty day_tests skips testing entirely, spelled either way", {
  dat <- make_2day_2station_data()
  dat$depth[300] <- 0
  dat$temp.water[301] <- NA

  for(none in list(character(0), c())) {
    res <- expect_silent(filter_2day(dat, day_tests=none))
    expect_equal(nrow(res$removed), 0)
    expect_equal(res$aln$n_days, 2)
  }

  # c() evaluates to NULL and is what the rest of the package uses to mean
  # "no tests", so specs() must accept it too rather than making two-station
  # the one model type that insists on character(0)
  expect_null(specs(mm_name('bayes_2s'), day_tests=c())$day_tests)
  expect_equal(specs(mm_name('bayes_2s'), day_tests=character(0))$day_tests, character(0))
})

test_that("pos_discharge is accepted but does nothing without a discharge column", {
  dat <- make_2day_2station_data()
  res <- expect_silent(filter_2day(dat, day_tests=c('complete_data','pos_discharge')))
  expect_equal(nrow(res$removed), 0)
})

test_that("one-station-only day_tests are rejected rather than silently applied", {
  # full_day tests the 4 AM / 28-hour diel window and would reject every
  # two-station day; even_timesteps would undo the gap tolerance mm_lag_2s()
  # provides on purpose
  dat <- make_2day_2station_data()
  aln <- suppressMessages(mm_align_2s(v(dat)))

  expect_error(
    mm_filter_valid_days_2s(dat, aln, day_tests=c('complete_data','full_day')),
    "may only include.*got full_day")
  expect_error(
    mm_filter_valid_days_2s(dat, aln, day_tests='even_timesteps'),
    "may only include.*got even_timesteps")
  expect_error(
    mm_filter_valid_days_2s(dat, aln, day_tests=c('nonsense')),
    "may only include.*got nonsense")
  expect_error(
    mm_filter_valid_days_2s(dat, aln, day_tests=42),
    "day_tests must be a character vector")
})

test_that("specs() rejects the same day_tests, at specs-creation time", {
  expect_error(specs(mm_name('bayes_2s'), day_tests=c('complete_data','full_day')),
               "may only include.*got full_day")
  expect_silent(specs(mm_name('bayes_2s'), day_tests='complete_data'))
})

test_that("specs(mm_name('bayes_2s')) defaults day_tests to complete_data + pos_depth", {
  sp <- specs(mm_name('bayes_2s'))
  expect_equal(sp$day_tests, c('complete_data','pos_depth'))

  # one-station's five-test default is deliberately not inherited
  expect_false('full_day' %in% sp$day_tests)
  expect_false('even_timesteps' %in% sp$day_tests)

  # day_start/day_end/required_timestep configure machinery two-station
  # doesn't use, and stay absent
  expect_null(sp$day_start)
  expect_null(sp$day_end)
  expect_null(sp$required_timestep)
})

test_that("mm_is_valid_day() accepts a day_tests subset that omits full_day and even_timesteps", {
  # regression test: those two tests' guards evaluate is.finite(timestep.days)
  # with `&` rather than `&&`, so timestep.days must be defined even when
  # neither test was requested (two-station's default requests neither)
  ply <- data.frame(
    solar.time = as.POSIXct("2050-06-01 06:00:00", tz="UTC") +
      as.difftime((0:23) * 1, units="hours"),
    depth = rep(0.5, 24), temp.water = rep(20, 24))

  expect_true(mm_is_valid_day(ply, day_tests=c('complete_data','pos_depth')))
  expect_true(mm_is_valid_day(ply, day_tests='complete_data'))
  expect_true(mm_is_valid_day(ply, day_tests='pos_depth'))

  ply$depth[3] <- NA
  expect_equal(mm_is_valid_day(ply, day_tests=c('complete_data','pos_depth')), "NAs in depth")
})
