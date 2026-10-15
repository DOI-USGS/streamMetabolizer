# a plain two-day aligned frame, date column included, for the validation tests
aligned_2day_plain <- function() {
  as.data.frame(suppressMessages(mm_align_data_2s(make_2day_2station_data())))
}


# mm_align_data_2s() ------------------------------------------------------

test_that("mm_align_data_2s() returns a unitless aligned frame with date first", {
  aligned <- suppressMessages(mm_align_data_2s(make_2day_2station_data()))
  expect_s3_class(aligned, 'aligned_2s')
  expect_identical(names(aligned), c(
    'date', 'solar.time', 'DO.obs.up', 'DO.sat.up', 'DO.obs.down', 'DO.sat.down',
    'light', 'depth', 'temp.water', 'travel.time'))

  aligned_from_unitted <- suppressMessages(mm_align_data_2s(formatted_raw_example_2s()$data))
  expect_false(any(vapply(unclass(aligned_from_unitted), unitted::is.unitted, logical(1))))
})

test_that("a dataset with nothing dropped records an empty, not a missing, removed", {
  aligned <- suppressMessages(mm_align_data_2s(make_2day_2station_data()))
  removed <- attr(aligned, 'removed')
  expect_true(is.data.frame(removed))
  expect_equal(nrow(removed), 0)
})

test_that("insufficient lead-in is reported by mm_align_data_2s()", {
  expect_error(mm_align_data_2s(make_2station_data(n=2)), "insufficient lead-in data")
})

test_that("mm_align_data_2s() records whether the light-window check was applied, and no fill setting", {
  # formatted data that pass the check carry an empty record, not a missing one
  raw <- as_raw_2s(make_2day_2station_data())
  formatted <- suppressMessages(mm_format_data_2s(raw$upstream, raw$downstream, raw$light))
  expect_equal(nrow(attr(formatted, 'light_check_failed')), 0)
  checked <- suppressMessages(mm_align_data_2s(formatted))
  expect_identical(attr(checked, 'light_check'), TRUE)
  expect_null(attr(checked, 'max_gap_hours'))

  # data carrying no light-window record are aligned without the check
  unchecked <- suppressMessages(mm_align_data_2s(make_2day_2station_data()))
  expect_identical(attr(unchecked, 'light_check'), FALSE)
})

test_that("mm_align_data_2s() output passes aligned-data validation", {
  expect_silent(mm_validate_data_2station(suppressMessages(mm_align_data_2s(subset_2station_days(two_station_example, 3)))))
})

test_that("each aligned row pairs downstream values with the upstream values one travel time earlier", {
  # checked against the raw input directly, not against the alignment engine:
  # the upstream partner of a row at solar.time t is the raw row at
  # t - travel.time, rounded to a whole number of timesteps
  expect_upstream_paired <- function(raw) {
    raw <- v(raw)
    aligned <- suppressMessages(mm_align_data_2s(raw))
    raw_t <- as.numeric(raw$solar.time)
    timestep_s <- stats::median(diff(raw_t))
    target_t <- as.numeric(aligned$solar.time) - round(aligned$travel.time * 86400 / timestep_s) * timestep_s
    up <- match(round(target_t), round(raw_t))
    down <- match(round(as.numeric(aligned$solar.time)), round(raw_t))
    expect_false(anyNA(up))
    expect_equal(aligned$DO.obs.up, raw$DO.obs.up[up])
    expect_equal(aligned$DO.sat.up, raw$DO.sat.up[up])
    expect_equal(aligned$DO.obs.down, raw$DO.obs.down[down])
  }

  # traceable values, so a same-time pairing can't pass by coincidence
  dat <- make_2station_data()
  dat$DO.obs.up <- seq_len(nrow(dat))
  dat$DO.sat.up <- seq_len(nrow(dat)) + 1000
  expect_upstream_paired(dat)

  expect_upstream_paired(formatted_raw_example_2s()$data)
})

test_that("NA travel.time is rejected with the affected dates", {
  dat <- make_2day_2station_data()
  dat$travel.time[c(100, 400)] <- NA
  expect_error(mm_align_data_2s(dat), "travel.time is NA on 2050-06-01, 2050-06-02; fill it")
})

test_that("errors from mm_align_data_2s() name no internal function", {
  offgrid <- make_2day_2station_data()
  offgrid$solar.time[10] <- offgrid$solar.time[10] + 60
  err <- tryCatch(mm_align_data_2s(offgrid), error=identity)
  expect_match(conditionMessage(err), "not on a regular timestep grid.*mm_format_data_2s")
  expect_null(conditionCall(err))

  err <- tryCatch(mm_align_data_2s(make_2station_data(n=2)), error=identity)
  expect_null(conditionCall(err))
})



# mm_modeled_rows_2s() ------------------------------------------------------

test_that("mm_modeled_rows_2s draws upstream from shift_idx and everything else from keep", {
  dat <- make_2day_2station_data()
  # distinct values per row make the indexing visible
  dat$DO.obs.up <- seq_len(nrow(dat))
  dat$DO.obs.down <- seq_len(nrow(dat)) + 10000
  alignment <- suppressMessages(mm_align_2s(v(dat), day_start_hour=mm_day_start_2s))

  modeled <- mm_modeled_rows_2s(dat, alignment)

  expect_equal(nrow(modeled), length(alignment$keep))
  expect_named(modeled, c('solar.time','DO.obs.up','DO.sat.up','DO.obs.down',
                          'DO.sat.down','light','depth','temp.water','travel.time'))
  expect_equal(modeled$DO.obs.up, dat$DO.obs.up[alignment$shift_idx])
  expect_equal(modeled$DO.obs.down, dat$DO.obs.down[alignment$keep])
  expect_equal(modeled$solar.time, dat$solar.time[alignment$keep])
})

test_that("mm_modeled_rows_2s strips units", {
  dat <- make_2day_2station_data()
  template <- mm_data(solar.time, DO.obs.up, DO.sat.up, DO.obs.down, DO.sat.down,
                      light, depth, temp.water, travel.time)
  for(col in names(template)) dat[[col]] <- u(dat[[col]], get_units(template[[col]]))
  dat$light <- u(v(dat$light), NA)
  alignment <- suppressMessages(mm_align_2s(v(dat), day_start_hour=mm_day_start_2s))

  expect_false(is.unitted(mm_modeled_rows_2s(dat, alignment)))
})

# [ method ------------------------------------------------------------------

test_that("row subsetting keeps the class, removed, and light-check flag", {
  aligned <- suppressMessages(mm_align_data_2s(make_2day_lit_boundary_data()))
  day1 <- aligned[aligned$date == aligned$date[1], , drop=FALSE]

  expect_s3_class(day1, 'aligned_2s')
  expect_identical(attr(day1, 'removed'), attr(aligned, 'removed'))
  expect_identical(attr(day1, 'light_check'), TRUE)
})

test_that("extracting a single column returns a vector, not an aligned frame", {
  aligned <- suppressMessages(mm_align_data_2s(make_2day_2station_data()))
  expect_type(aligned[, 'depth'], 'double')
  expect_s3_class(aligned['depth'], 'aligned_2s')
})

test_that("an absent attribute stays absent (unknown) through subsetting", {
  aligned <- mm_as_aligned_2s(aligned_2day_plain())
  sub <- aligned[1:10, ]
  expect_null(attr(sub, 'removed'))
  expect_null(attr(sub, 'light_check'))
})


# mm_as_aligned_2s() ------------------------------------------------------

test_that("mm_as_aligned_2s() accepts a hand-built frame with or without a date column", {
  plain <- aligned_2day_plain()

  with_date <- mm_as_aligned_2s(plain)
  expect_s3_class(with_date, 'aligned_2s')
  expect_identical(with_date$date, plain$date)

  without_date <- mm_as_aligned_2s(plain[setdiff(names(plain), 'date')])
  expect_identical(names(without_date), names(plain))
  expect_identical(without_date$date, plain$date)
})

test_that("mm_as_aligned_2s() validates its input", {
  plain <- aligned_2day_plain()
  expect_error(mm_as_aligned_2s(plain[-nrow(plain), ]), "same number of rows")
})

test_that("mm_as_aligned_2s() keeps a user-supplied date, and so rejects a wrong one", {
  plain <- aligned_2day_plain()
  plain$date[1] <- plain$date[1] + 1
  expect_error(mm_as_aligned_2s(plain), "pass that hour as day_start_hour")
})

test_that("mm_as_aligned_2s() records dropped days and the light-check flag as unknown", {
  aligned <- suppressMessages(mm_align_data_2s(with_light_check_2s(make_2day_2station_data())))
  rebuilt <- mm_as_aligned_2s(aligned)
  expect_null(attr(rebuilt, 'removed'))
  expect_null(attr(rebuilt, 'light_check'))
})

test_that("mm_as_aligned_2s() validates the frame with its date column", {
  seen <- NULL
  real_validate <- mm_validate_data
  testthat::local_mocked_bindings(mm_validate_data=function(data, ...) {
    seen <<- names(data)
    real_validate(data, ...)
  })
  mm_as_aligned_2s(aligned_2day_plain())
  expect_true('date' %in% seen)
})

test_that("mm_as_aligned_2s() rejects NA travel.time with the affected dates", {
  plain <- aligned_2day_plain()
  plain$travel.time[which(plain$date == plain$date[1])[5]] <- NA
  err <- expect_error(
    mm_as_aligned_2s(plain),
    paste0("travel.time is NA on ", plain$date[1], "; fill it or remove those days"))
  expect_false(grepl("mm_[a-z_0-9]+\\(|missing value", conditionMessage(err)))
})

test_that("mm_as_aligned_2s() rejects a solar.time that is not UTC", {
  plain <- aligned_2day_plain()
  attr(plain$solar.time, 'tzone') <- 'America/Denver'
  expect_error(mm_as_aligned_2s(plain), "timezone 'UTC'")
  expect_error(mm_as_aligned_2s(plain[setdiff(names(plain), 'date')]), "timezone 'UTC'")
})

# aligned-data validation ----------------------------------------------------

test_that("aligned-data validation rejects each malformed frame", {
  plain <- aligned_2day_plain()
  aligned <- suppressMessages(mm_align_data_2s(make_2day_2station_data()))
  with_flag <- function(df, light_check) {
    new_aligned_2s(df, light_check=light_check, day_start_hour=mm_day_start_2s)
  }

  cases <- list(
    'missing required column' = list(
      data=with_flag(plain[setdiff(names(plain), 'travel.time')], NULL),
      error="missing these columns: travel.time"),
    'zero rows' = list(
      data=with_flag(plain[0, ], NULL),
      error="no rows"),
    'non-POSIXct solar.time' = list(
      data={d <- plain; d$solar.time <- format(d$solar.time); with_flag(d, NULL)},
      error="non-NA 'solar.time' column of class POSIXct"),
    'NA solar.time' = list(
      data={d <- plain; d$solar.time[1] <- NA; with_flag(d, NULL)},
      error="non-NA 'solar.time' column of class POSIXct"),
    'date not matching solar.time' = list(
      data={d <- plain; d$date[1] <- d$date[1] + 1; with_flag(d, NULL)},
      error="24-hour day, starting at hour 6"),
    'NA date' = list(
      data={d <- plain; d$date[1] <- NA; with_flag(d, NULL)},
      error="non-NA 'date' column"),
    'non-Date date' = list(
      data={d <- plain; d$date <- as.character(d$date); with_flag(d, NULL)},
      error="class Date"),
    'unsorted rows' = list(
      data=with_flag(plain[c(2, 1, 3:nrow(plain)), ], NULL),
      error="sorted"),
    'rows out of order after rbind' = list(
      data=rbind(aligned, aligned[aligned$date == aligned$date[1], ]),
      error="sorted"),
    'sorted but duplicated rows' = list(
      data=with_flag(plain[c(1, 1:nrow(plain)), ], NULL),
      error="no duplicate timestamps"),
    'a day with fewer rows' = list(
      data=with_flag(plain[-nrow(plain), ], NULL),
      error="same number of rows"),
    'irregular timestep within a day' = list(
      data={d <- plain; d$solar.time[100] <- d$solar.time[100] + as.difftime(2, units='mins'); with_flag(d, NULL)},
      error="regular timestep"),
    'a 5-minute day next to a 15-minute day' = list(
      data={
        steps <- function(start, minutes) as.POSIXct(start, tz='UTC') + as.difftime((0:95) * minutes, units='mins')
        d <- plain[1:192, ]
        d$solar.time <- c(steps('2050-06-01 06:00', 5), steps('2050-06-02 06:00', 15))
        d$date <- mm_date_2s(d$solar.time, mm_day_start_2s)
        with_flag(d, NULL)},
      error="single regular timestep; found steps from 5 to 15 minutes"),
    'non-positive travel.time' = list(
      data={d <- plain; d$travel.time[5] <- 0; with_flag(d, NULL)},
      error="travel.time must be > 0"),
    'malformed light-check flag' = list(
      data=with_flag(plain, 'yes'),
      error="malformed 'light_check' attribute"))

  for(case in names(cases)) {
    expect_error(mm_validate_data_2station(cases[[case]]$data), cases[[case]]$error, info=case)
  }
})


# rbind method -------------------------------------------------------------

test_that("rbind of aligned pieces gives an aligned frame with unknown removed days", {
  aligned <- suppressMessages(mm_align_data_2s(make_2day_2station_data()))
  dates <- unique(aligned$date)
  combined <- rbind(aligned[aligned$date == dates[1], ], aligned[aligned$date == dates[2], ])

  expect_s3_class(combined, 'aligned_2s')
  expect_null(attr(combined, 'removed'))
  expect_identical(attr(combined, 'light_check'), attr(aligned, 'light_check'))
  expect_identical(unclass(combined)[names(aligned)], unclass(aligned)[names(aligned)])
  expect_silent(mm_validate_data_2station(combined))
})

test_that("rbind keeps the light-check flag only when every piece has the same one", {
  aligned <- suppressMessages(mm_align_data_2s(with_light_check_2s(make_2day_2station_data())))
  dates <- unique(aligned$date)
  day1 <- aligned[aligned$date == dates[1], ]
  day2 <- aligned[aligned$date == dates[2], ]
  expect_identical(attr(rbind(day1, day2), 'light_check'), TRUE)

  attr(day2, 'light_check') <- FALSE
  expect_null(attr(rbind(day1, day2), 'light_check'))

  attr(day2, 'light_check') <- NULL
  expect_null(attr(rbind(day1, day2), 'light_check'))
})

test_that("rbind with a plain data.frame depends on which comes first", {
  # base R dispatches rbind to this method whenever the first data.frame
  # argument is aligned, and to rbind.data.frame otherwise
  aligned <- suppressMessages(mm_align_data_2s(make_2day_2station_data()))
  dates <- unique(aligned$date)
  day1 <- aligned[aligned$date == dates[1], ]
  day2_plain <- as.data.frame(aligned[aligned$date == dates[2], ])

  aligned_first <- rbind(day1, day2_plain)
  expect_s3_class(aligned_first, 'aligned_2s')
  expect_null(attr(aligned_first, 'removed'))
  expect_null(attr(aligned_first, 'light_check'))

  plain_first <- rbind(day2_plain, day1)
  expect_false(inherits(plain_first, 'aligned_2s'))
})

test_that("rbind passes rbind.data.frame's named arguments through", {
  aligned <- suppressMessages(mm_align_data_2s(make_2day_2station_data()))
  combined <- rbind(aligned[1:5, ], aligned[6:10, ], make.row.names=FALSE)
  expect_s3_class(combined, 'aligned_2s')
  expect_equal(nrow(combined), 10)
})


# as.data.frame and dplyr -----------------------------------------------------

test_that("as.data.frame() returns a plain data.frame with no alignment record", {
  aligned <- suppressMessages(mm_align_data_2s(make_2day_lit_boundary_data()))
  plain <- as.data.frame(aligned)
  expect_identical(class(plain), 'data.frame')
  expect_null(attr(plain, 'removed'))
  expect_null(attr(plain, 'light_check'))
  expect_identical(lapply(plain, identity), lapply(unclass(aligned), identity))
})

test_that("dplyr verbs keep the alignment record, but bind_rows() of separate alignments doesn't", {
  a1 <- suppressMessages(mm_align_data_2s(make_2day_lit_boundary_data()))
  later <- make_2day_2station_data()
  later$solar.time <- later$solar.time + as.difftime(10, units='days')
  a2 <- suppressMessages(mm_align_data_2s(later))

  filtered <- dplyr::filter(a1, date == min(date))
  expect_s3_class(filtered, 'aligned_2s')
  expect_identical(attr(filtered, 'removed'), attr(a1, 'removed'))
  expect_identical(attr(filtered, 'light_check'), TRUE)

  bound <- dplyr::bind_rows(a1, a2)
  expect_null(attr(bound, 'removed'))
  expect_null(attr(bound, 'light_check'))
})



# day_start_hour -------------------------------------------------------------

test_that("a non-default day_start_hour moves the day boundary", {
  # make_ts_data(): hourly rows from 2050-06-01 03:00, 3-row lag. With days
  # starting at midnight, 2050-06-02 is the one complete day
  aligned <- suppressMessages(mm_align_data_2s(make_ts_data(), day_start_hour=0))
  expect_identical(unique(aligned$date), as.Date('2050-06-02'))
  expect_identical(aligned$solar.time[1], as.POSIXct('2050-06-02 00:00:00', tz='UTC'))
  expect_identical(attr(aligned, 'day_start_hour'), 0)
  expect_silent(mm_validate_data_2station(aligned))

  default <- suppressMessages(mm_align_data_2s(make_ts_data()))
  expect_identical(attr(default, 'day_start_hour'), 6)

  # the hour is validated: 0 is allowed, 24 and malformed values are not
  expect_error(mm_align_data_2s(make_ts_data(), day_start_hour=24), "day_start_hour must be >= 0 and < 24")
  expect_error(mm_check_day_start_hour_2s(NA_real_), "single non-NA number")
  expect_error(mm_check_day_start_hour_2s(-1), ">= 0 and < 24")
  expect_silent(mm_check_day_start_hour_2s(0))
})

test_that("mm_as_aligned_2s() labels and checks dates with its own day_start_hour", {
  plain <- as.data.frame(suppressMessages(mm_align_data_2s(make_ts_data(), day_start_hour=0)))
  expect_identical(attr(mm_as_aligned_2s(plain, day_start_hour=0), 'day_start_hour'), 0)
  expect_identical(
    as.data.frame(mm_as_aligned_2s(plain[names(plain) != 'date'], day_start_hour=0)), plain)
})

test_that("rbind carries the day start hour unless pieces disagree", {
  a6 <- suppressMessages(mm_align_data_2s(make_ts_data()))
  later <- make_ts_data()
  later$solar.time <- later$solar.time + as.difftime(10, units='days')
  a0 <- suppressMessages(mm_align_data_2s(later, day_start_hour=0))
  expect_error(mm_validate_data_2station(rbind(a6, a0)), "no recorded day_start_hour")

  # a plain piece has no hour, so it doesn't conflict
  a6_later <- suppressMessages(mm_align_data_2s(later))
  expect_silent(mm_validate_data_2station(rbind(a6, as.data.frame(a6_later))))
})


# light-window check -----------------------------------------------------------

# n_days whole calendar days from 2050-06-01 00:00 at timestep_min, with light
# from daylight_2s() (lit 07:00-19:00, exactly 0 otherwise)
light_series_2s <- function(timestep_min=60, n_days=4, travel_time=3/24,
                            start='2050-06-01 00:00') {
  solar.time <- as.POSIXct(start, tz='UTC') +
    as.difftime(seq(0, n_days * 1440 - timestep_min, by=timestep_min), units='mins')
  data.frame(
    solar.time=solar.time, DO.obs.up=9, DO.sat.up=10, DO.obs.down=8.8, DO.sat.down=9.9,
    light=daylight_2s(solar.time), depth=0.5, temp.water=20, travel.time=travel_time)
}
at <- function(dat, time) which(dat$solar.time == as.POSIXct(time, tz='UTC'))
light_check <- function(dat, day_start_hour=6) {
  mm_light_check_2s(dat$solar.time, dat$light, dat$travel.time, day_start_hour)
}

test_that("a day start at night passes; a daytime start fails every day, and the error points to travel.time and day_start_hour", {
  dat <- light_series_2s()
  # a 3-hour window from a midnight start reaches back only to 21:00
  expect_equal(nrow(light_check(dat, 0)), 0)

  # with a lag of 0 the window is the row itself, so even a noon start passes
  no_lag <- light_series_2s(travel_time=1/1440)
  expect_true(all(mm_lag_2s(no_lag$solar.time, no_lag$travel.time)$lag == 0))
  expect_equal(nrow(light_check(no_lag, 12)), 0)

  # from a noon start it reaches the previous day's lit morning, so every day
  # after the first fails and none is left
  failed <- light_check(dat, 12)
  expect_equal(failed$date, unique(mm_date_2s(dat$solar.time, 12))[-1])
  expect_error(
    suppressMessages(mm_align_data_2s(with_light_check_2s(dat, 12), day_start_hour=12)),
    paste0(
      "^Every day was dropped: each day's travel-time light window reaches light from an earlier day, ",
      "which would contaminate the day's light proportions\\. Check that `travel\\.time` is in days and ",
      "not too long, and move `day_start_hour` to a time with no light in the hours before it\\.$"))

  # the general message when no day failed the light check
  expect_error(
    suppressMessages(mm_align_data_2s(make_2station_data(n=100))),
    "^no complete days remain: every day was dropped \\(see the messages above\\)$")
})

test_that("NA light in an earlier day fails the day with its own message", {
  dat <- light_series_2s()
  dat$light[at(dat, '2050-06-02 04:00')] <- NA
  expect_message(
    alignment <- mm_align_2s(dat, day_start_hour=6, light_check_failed=light_check(dat)),
    "whose travel-time light windows reach light unknown in an earlier day: 2050-06-02 \\(2 row\\(s\\)\\)")
  expect_equal(
    alignment$removed$errors[alignment$removed$date == as.Date('2050-06-02')],
    "light unknown in an earlier day of the light window (2 row(s))")
})

test_that("a lit lead-in row fails the first modeled day, and a dark one passes", {
  dat <- light_series_2s(start='2050-06-01 03:00', n_days=2)
  expect_equal(nrow(light_check(dat)), 0)
  dat$light[at(dat, '2050-06-01 04:00')] <- 300
  expect_equal(light_check(dat)$date, as.Date('2050-06-01'))
})

test_that("a window spanning more than one earlier day is checked across all of them", {
  dat <- light_series_2s(n_days=5, travel_time=1.5)
  dat$light <- 0
  # 36 hours before 06:00 on 06-04 is 18:00 on 06-02, two days earlier
  dat$light[at(dat, '2050-06-02 18:00')] <- 1
  expect_true(as.Date('2050-06-04') %in% light_check(dat)$date)
})

test_that("mm_light_check_2s() matches a brute-force check of every window", {
  set.seed(4475)
  # each class runs at random timesteps, with rows dropped to leave gaps, a
  # random day start hour, and per-day travel times that reach no row before
  # the day, rows of the night before it, or rows of the afternoon before it
  classes <- list(
    ordinary=function(st) sample(c(0, 0, 0, 5, 100, NA), length(st), replace=TRUE),
    zero_or_na=function(st) sample(c(0, 0, 0, NA), length(st), replace=TRUE),
    tiny=function(st) sample(c(0, 0, 0, 1e-12), length(st), replace=TRUE),
    exact_zero_at_boundary=function(st) daylight_2s(st))
  brute_force <- function(solar.time, light, travel.time, hour) {
    lagged <- mm_lag_2s(solar.time, travel.time)
    day <- mm_date_2s(solar.time, hour)
    lit <- unknown <- logical(length(solar.time))
    for(i in which(lagged$has_leadin)) {
      j <- lagged$shift_idx[i]:i
      j <- j[day[j] < day[i]]
      lit[i] <- any(light[j] > 0, na.rm=TRUE)
      unknown[i] <- anyNA(light[j])
    }
    dates <- sort(unique(day[lit | unknown]))
    reason <- c('unknown', 'lit')[(dates %in% day[lit]) + 1]
    n_rows <- vapply(seq_along(dates), function(k) {
      sum(day == dates[k] & if(reason[k] == 'lit') lit else unknown)
    }, integer(1))
    data.frame(date=dates, reason=reason, n_rows=n_rows, stringsAsFactors=FALSE)
  }

  covered <- c()
  for(class in names(classes)) {
    got <- expected <- list()
    any_fail <- any_pass <- FALSE
    for(trial in 1:3) {
      timestep_min <- sample(c(5, 15, 30, 60), 1)
      grid <- as.POSIXct('2050-06-01 00:00', tz='UTC') +
        as.difftime(seq(0, 6 * 1440 - timestep_min, by=timestep_min), units='mins')
      solar.time <- sort(sample(grid, round(0.85 * length(grid))))
      hour <- sample(c(0, 2, 4, 6), 1)
      day <- mm_date_2s(solar.time, hour)
      days <- unique(day)
      reach <- sample(c('none', 'night', 'afternoon'), length(days), replace=TRUE)
      travel_by_day <- vapply(reach, function(r) switch(r,
        none=1e-4, night=runif(1, 0.01, 0.2), afternoon=runif(1, 0.5, 0.7)), numeric(1))
      travel.time <- unname(travel_by_day[match(day, days)])
      light <- classes[[class]](solar.time)

      got[[trial]] <- mm_light_check_2s(solar.time, light, travel.time, hour)
      expected[[trial]] <- brute_force(solar.time, light, travel.time, hour)
      any_fail <- any_fail || nrow(expected[[trial]]) > 0
      any_pass <- any_pass || any(!(days[-1] %in% expected[[trial]]$date))
    }
    expect_equal(got, expected, info=class)
    covered[class] <- any_fail && any_pass
  }
  # every class has at least one failing day and one passing day after the first
  expect_true(all(covered), info=paste(names(covered)[!covered], collapse=', '))
})

test_that("the light-window drop message names at most 10 days, and removed records them all", {
  dat <- light_series_2s(n_days=14)
  # lit at 04:00 on 06-01 through 06-12, failing the days that start two hours later
  dat$light[format(dat$solar.time, '%H') == '04' & dat$solar.time < as.POSIXct('2050-06-13', tz='UTC')] <- 50
  msgs <- testthat::capture_messages(aligned <- mm_align_data_2s(with_light_check_2s(dat)))
  light_msg <- grep('light above zero', msgs, value=TRUE)
  expect_length(light_msg, 1)
  expect_equal(lengths(regmatches(light_msg, gregexpr('[0-9]{4}-[0-9]{2}-[0-9]{2}', light_msg))), 10)
  expect_match(light_msg, '^dropping 12 day\\(s\\) .*, and 2 more\n$')
  expect_equal(sum(grepl('light window', attr(aligned, 'removed')$errors)), 12)
})
