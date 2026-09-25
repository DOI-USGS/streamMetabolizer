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

  raw <- two_station_raw_example
  formatted <- suppressMessages(mm_format_data_2s(raw$upstream, raw$downstream, raw$light))
  aligned_from_unitted <- suppressMessages(mm_align_data_2s(formatted))
  expect_false(any(vapply(unclass(aligned_from_unitted), unitted::is.unitted, logical(1))))
})

test_that("mm_align_data_2s() reports the days it drops", {
  expect_message(mm_align_data_2s(make_2day_ceiling_data()), "exceeds the .* ceiling")
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

test_that("mm_align_data_2s() records its travel-time ceiling and no fill setting", {
  aligned <- suppressMessages(mm_align_data_2s(make_2day_2station_data(), max_travel_time_days=0.3))
  expect_identical(attr(aligned, 'max_travel_time_days'), 0.3)
  expect_null(attr(aligned, 'max_gap_hours'))
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

  raw <- two_station_raw_example
  expect_upstream_paired(suppressMessages(mm_format_data_2s(raw$upstream, raw$downstream, raw$light)))
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


# [ method ------------------------------------------------------------------

test_that("row subsetting keeps the class, removed, and travel-time ceiling", {
  aligned <- suppressMessages(mm_align_data_2s(make_2day_ceiling_data()))
  day1 <- aligned[aligned$date == aligned$date[1], , drop=FALSE]

  expect_s3_class(day1, 'aligned_2s')
  expect_identical(attr(day1, 'removed'), attr(aligned, 'removed'))
  expect_identical(attr(day1, 'max_travel_time_days'), attr(aligned, 'max_travel_time_days'))
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
  expect_null(attr(sub, 'max_travel_time_days'))
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
  expect_error(mm_as_aligned_2s(plain), "06:00-06:00 day")
})

test_that("mm_as_aligned_2s() records dropped days and the ceiling as unknown, and so skips the ceiling check", {
  aligned <- suppressMessages(mm_align_data_2s(make_2day_2station_data()))
  rebuilt <- mm_as_aligned_2s(aligned)
  expect_null(attr(rebuilt, 'removed'))
  expect_null(attr(rebuilt, 'max_travel_time_days'))

  capped <- suppressMessages(mm_align_data_2s(make_2day_2station_data(), max_travel_time_days=0.2))
  capped$travel.time[5] <- 0.3
  expect_silent(mm_as_aligned_2s(capped))
})


# aligned-data validation ----------------------------------------------------

test_that("aligned-data validation rejects each malformed frame", {
  plain <- aligned_2day_plain()
  aligned <- suppressMessages(mm_align_data_2s(make_2day_2station_data()))
  with_ceiling <- function(df, ceiling_days) {
    new_aligned_2s(df, max_travel_time_days=ceiling_days)
  }

  cases <- list(
    'missing required column' = list(
      data=with_ceiling(plain[setdiff(names(plain), 'travel.time')], NULL),
      error="missing these columns: travel.time"),
    'zero rows' = list(
      data=with_ceiling(plain[0, ], NULL),
      error="no rows"),
    'non-POSIXct solar.time' = list(
      data={d <- plain; d$solar.time <- format(d$solar.time); with_ceiling(d, NULL)},
      error="non-NA 'solar.time' column of class POSIXct"),
    'NA solar.time' = list(
      data={d <- plain; d$solar.time[1] <- NA; with_ceiling(d, NULL)},
      error="non-NA 'solar.time' column of class POSIXct"),
    'date not matching solar.time' = list(
      data={d <- plain; d$date[1] <- d$date[1] + 1; with_ceiling(d, NULL)},
      error="06:00-06:00 day"),
    'NA date' = list(
      data={d <- plain; d$date[1] <- NA; with_ceiling(d, NULL)},
      error="non-NA 'date' column"),
    'non-Date date' = list(
      data={d <- plain; d$date <- as.character(d$date); with_ceiling(d, NULL)},
      error="class Date"),
    'unsorted rows' = list(
      data=with_ceiling(plain[c(2, 1, 3:nrow(plain)), ], NULL),
      error="sorted"),
    'rows out of order after rbind' = list(
      data=rbind(aligned, aligned[aligned$date == aligned$date[1], ]),
      error="sorted"),
    'sorted but duplicated rows' = list(
      data=with_ceiling(plain[c(1, 1:nrow(plain)), ], NULL),
      error="no duplicate timestamps"),
    'a day with fewer rows' = list(
      data=with_ceiling(plain[-nrow(plain), ], NULL),
      error="same number of rows"),
    'irregular timestep within a day' = list(
      data={d <- plain; d$solar.time[100] <- d$solar.time[100] + as.difftime(2, units='mins'); with_ceiling(d, NULL)},
      error="regular timestep"),
    'a 5-minute day next to a 15-minute day' = list(
      data={
        steps <- function(start, minutes) as.POSIXct(start, tz='UTC') + as.difftime((0:95) * minutes, units='mins')
        d <- plain[1:192, ]
        d$solar.time <- c(steps('2050-06-01 06:00', 5), steps('2050-06-02 06:00', 15))
        d$date <- mm_date_2s(d$solar.time)
        with_ceiling(d, NULL)},
      error="single regular timestep; found steps from 5 to 15 minutes"),
    'non-positive travel.time' = list(
      data={d <- plain; d$travel.time[5] <- 0; with_ceiling(d, NULL)},
      error="travel.time must be > 0"),
    'travel.time above the stored ceiling' = list(
      data={d <- plain; d$travel.time[5] <- 0.3; with_ceiling(d, 0.2)},
      error="ceiling"))

  for(case in names(cases)) {
    expect_error(mm_validate_data_2station(cases[[case]]$data), cases[[case]]$error, info=case)
  }
})

test_that("travel.time exactly at the stored ceiling passes", {
  d <- aligned_2day_plain()
  d$travel.time[5] <- 0.3
  expect_silent(mm_validate_data_2station(new_aligned_2s(d, max_travel_time_days=0.3)))
})


# rbind method -------------------------------------------------------------

test_that("rbind of aligned pieces gives an aligned frame with unknown removed days", {
  aligned <- suppressMessages(mm_align_data_2s(make_2day_2station_data()))
  dates <- unique(aligned$date)
  combined <- rbind(aligned[aligned$date == dates[1], ], aligned[aligned$date == dates[2], ])

  expect_s3_class(combined, 'aligned_2s')
  expect_null(attr(combined, 'removed'))
  expect_identical(attr(combined, 'max_travel_time_days'), attr(aligned, 'max_travel_time_days'))
  expect_identical(unclass(combined)[names(aligned)], unclass(aligned)[names(aligned)])
  expect_silent(mm_validate_data_2station(combined))
})

test_that("rbind keeps the ceiling only when every piece has the same one", {
  aligned <- suppressMessages(mm_align_data_2s(make_2day_2station_data()))
  dates <- unique(aligned$date)
  day1 <- aligned[aligned$date == dates[1], ]
  day2 <- aligned[aligned$date == dates[2], ]

  attr(day2, 'max_travel_time_days') <- 0.3
  expect_null(attr(rbind(day1, day2), 'max_travel_time_days'))

  attr(day2, 'max_travel_time_days') <- NULL
  expect_null(attr(rbind(day1, day2), 'max_travel_time_days'))
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
  expect_null(attr(aligned_first, 'max_travel_time_days'))

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
  aligned <- suppressMessages(mm_align_data_2s(make_2day_ceiling_data()))
  plain <- as.data.frame(aligned)
  expect_identical(class(plain), 'data.frame')
  expect_null(attr(plain, 'removed'))
  expect_null(attr(plain, 'max_travel_time_days'))
  expect_identical(lapply(plain, identity), lapply(unclass(aligned), identity))
})

test_that("dplyr verbs keep the alignment record, but bind_rows() of separate alignments doesn't", {
  a1 <- suppressMessages(mm_align_data_2s(make_2day_ceiling_data(), max_travel_time_days=0.42))
  later <- make_2day_2station_data()
  later$solar.time <- later$solar.time + as.difftime(10, units='days')
  a2 <- suppressMessages(mm_align_data_2s(later, max_travel_time_days=0.3))

  filtered <- dplyr::filter(a1, date == min(date))
  expect_s3_class(filtered, 'aligned_2s')
  expect_identical(attr(filtered, 'removed'), attr(a1, 'removed'))
  expect_identical(attr(filtered, 'max_travel_time_days'), 0.42)

  bound <- dplyr::bind_rows(a1, a2)
  expect_null(attr(bound, 'removed'))
  expect_null(attr(bound, 'max_travel_time_days'))
})
