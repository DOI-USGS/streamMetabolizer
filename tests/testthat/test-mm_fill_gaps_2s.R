# the data without the fill_share table, for comparing data values alone
drop_fill_share <- function(x) {
  attr(x, 'fill_share') <- NULL
  x
}

# a fill-share table's columns as a plain data.frame, for comparing values
plain_fill_share <- function(fill_share) {
  data.frame(date=fill_share$date, share=fill_share$share)
}

# gap sizing ----------------------------------------------------------------

test_that("missing-row gaps are filled up to the tolerance and left alone beyond it", {
  dat <- make_gappy_data(n=40)

  for(gap in c(1, 4)) {
    gappy <- dat[-(10:(9 + gap)), ]
    filled <- suppressMessages(mm_fill_gaps_2s(gappy))
    expect_equal(nrow(filled), nrow(dat))
    # the inserted rows reproduce what was removed: every column here is
    # linear in row number, so linear interpolation is exact
    expect_equal(filled$solar.time, dat$solar.time)
    expect_equal(filled$DO.obs.down, dat$DO.obs.down)
  }

  # 5 bins is one past the tolerance: nothing is inserted, and the hole is
  # left for the day-window check to act on
  gappy <- dat[-(10:14), ]
  filled <- suppressMessages(mm_fill_gaps_2s(gappy))
  expect_equal(nrow(filled), nrow(gappy))
  rownames(gappy) <- NULL # filling renumbers rows; the values are the point
  expect_equal(drop_fill_share(filled), gappy)
})

test_that("missing-value gaps use the same tolerance as missing-row gaps", {
  dat <- make_gappy_data(n=40)

  gappy <- dat
  gappy$DO.obs.up[10:13] <- NA # 4 bins: fillable
  gappy$temp.water[20:24] <- NA # 5 bins: not
  filled <- suppressMessages(mm_fill_gaps_2s(gappy))

  expect_equal(nrow(filled), nrow(dat))
  expect_false(any(is.na(filled$DO.obs.up)))
  expect_equal(filled$DO.obs.up, dat$DO.obs.up)
  expect_equal(sum(is.na(filled$temp.water)), 5)
})

test_that("gap size is measured in time, not in rows", {
  # two NA rows that are adjacent in row order but sit on either side of a
  # multi-day outage. Measured in rows the run looks like 2; measured on the
  # bin grid its bracketing observations are hundreds of bins apart, so
  # filling it would interpolate straight across the outage
  early <- make_gappy_data(n=10, start="2050-06-01 00:00:00")
  late <- make_gappy_data(n=10, start="2050-06-05 00:00:00")
  dat <- rbind(early, late)
  dat$DO.obs.down[10] <- NA
  dat$DO.obs.down[11] <- NA

  filled <- suppressMessages(mm_fill_gaps_2s(dat))

  expect_equal(sum(is.na(filled$DO.obs.down)), 2)
  expect_equal(nrow(filled), nrow(dat))
})

test_that("a tolerance finer than the timestep fills nothing, and says so", {
  dat <- make_gappy_data(n=20)[-(10:11), ]
  expect_message(
    filled <- mm_fill_gaps_2s(dat, max_gap_hours=0.1),
    'shorter than one timestep')
  expect_equal(drop_fill_share(filled), dat)

  # the fill-share table is still set, with nothing filled
  expect_equal(plain_fill_share(attr(filled, 'fill_share')),
               data.frame(date=as.Date('2050-05-31'), share=0))
})


# edges ---------------------------------------------------------------------

test_that("gaps at the start and end of the record are never filled", {
  dat <- make_gappy_data(n=40)

  # missing values with an observation on only one side
  gappy <- dat
  gappy$DO.obs.up[1:2] <- NA
  gappy$DO.obs.up[39:40] <- NA
  filled <- suppressMessages(mm_fill_gaps_2s(gappy))
  expect_true(all(is.na(filled$DO.obs.up[1:2])))
  expect_true(all(is.na(filled$DO.obs.up[39:40])))

  # missing rows off either end aren't gaps at all -- there is no bracketing
  # observation to notice them by, so the record simply starts/ends later
  trimmed <- dat[3:38, ]
  filled <- suppressMessages(mm_fill_gaps_2s(trimmed))
  expect_equal(nrow(filled), nrow(trimmed))
})


# the no-NA-row invariant ---------------------------------------------------

test_that("inserted rows that cannot be fully filled are dropped again", {
  # a 2-bin missing-row gap sitting inside a 12-bin missing-value run in one
  # column: the rows are insertable, but DO.sat.up can't be interpolated for
  # them, so keeping them would satisfy the day row count while carrying NA
  dat <- make_gappy_data(n=40)
  gappy <- dat
  gappy$DO.sat.up[15:26] <- NA
  gappy <- gappy[-(20:21), ]

  filled <- suppressMessages(mm_fill_gaps_2s(gappy))

  expect_equal(nrow(filled), nrow(gappy))
  expect_equal(sum(is.na(filled$DO.sat.up)), 10)
  # the invariant itself: no row anywhere holds a partially filled record
  interpolated_cols <- setdiff(names(filled), 'solar.time')
  na_rows <- rowSums(is.na(filled[interpolated_cols])) > 0
  expect_equal(sum(na_rows), 10)
  expect_true(all(rowSums(is.na(filled[na_rows, interpolated_cols])) == 1))
})

test_that("filling never introduces an NA-bearing row on real data", {
  # the raw downstream series is genuinely gappy in both forms at once: 200-odd
  # runs of missing timesteps and 30-odd runs of NA depth, at every length from
  # one bin to several days
  data(two_station_raw_example)
  down <- unitted::v(two_station_raw_example$downstream)
  names(down)[names(down) == 'timestamp'] <- 'solar.time'
  down$solar.time <- mm_snap_to_bin_2s(down$solar.time)

  filled <- suppressMessages(mm_fill_gaps_2s(down, warn_fill_share=1))

  # original rows pass through untouched, so they are exactly the timestamps
  # the input already had; anything else was inserted
  inserted <- !(as.numeric(filled$solar.time) %in% as.numeric(down$solar.time))
  expect_gt(sum(inserted), 0)
  expect_false(any(is.na(filled[inserted, ])))

  # and the two-tier policy is intact: over-tolerance gaps are still NA or
  # still absent, rather than having been quietly bridged
  expect_true(any(is.na(filled)))
  expect_lt(nrow(filled), diff(range(as.numeric(filled$solar.time))) / (15 * 60) + 1)
})


# idempotence and units -----------------------------------------------------

test_that("filling is idempotent", {
  dat <- make_gappy_data(n=40)[-(10:12), ]
  once <- suppressMessages(mm_fill_gaps_2s(dat))
  twice <- suppressMessages(mm_fill_gaps_2s(once))
  expect_equal(once, twice)
  # the second pass has nothing to report
  expect_silent(mm_fill_gaps_2s(once))
})

test_that("units survive filling, in whichever shape they arrived", {
  data(two_station_example)
  units_in <- unitted::get_units(two_station_example)

  # (a) the common shape: a plain data.frame whose columns are unitted
  gappy <- two_station_example[-(10:12), ]
  filled <- suppressMessages(mm_fill_gaps_2s(gappy))
  expect_equal(unitted::get_units(filled), units_in)
  expect_false(unitted::is.unitted(filled))
  expect_equal(nrow(filled), nrow(two_station_example))

  # (b) a true unitted data.frame stays one
  unitted_in <- unitted::u(unitted::v(two_station_example), units_in)[-(10:12), ]
  filled <- suppressMessages(mm_fill_gaps_2s(unitted_in))
  expect_true(unitted::is.unitted(filled))
  expect_equal(unitted::get_units(filled), units_in)

  # (c) a plain data.frame stays plain
  plain <- unitted::v(two_station_example)[-(10:12), ]
  filled <- suppressMessages(mm_fill_gaps_2s(plain))
  expect_false(any(vapply(filled, unitted::is.unitted, logical(1))))
})


# preconditions -------------------------------------------------------------

test_that("mm_fill_gaps_2s enforces the same grid preconditions as mm_lag_2s", {
  dat <- make_gappy_data(n=20)

  offgrid <- dat
  offgrid$solar.time <- offgrid$solar.time + as.difftime(3, units="mins")
  expect_error(mm_fill_gaps_2s(offgrid), 'not on a regular timestep grid')

  duped <- rbind(dat, dat[10, ])
  duped <- duped[order(duped$solar.time), ]
  expect_error(mm_fill_gaps_2s(duped), 'more than one row in the same timestep bin')

  unsorted <- dat[c(2, 1, 3:20), ]
  expect_error(mm_fill_gaps_2s(unsorted), 'sorted ascending')
})

test_that("max_gap_hours and warn_fill_share are validated", {
  dat <- make_gappy_data(n=20)
  expect_error(mm_fill_gaps_2s(dat, max_gap_hours=0), 'must be > 0')
  expect_error(mm_fill_gaps_2s(dat, max_gap_hours=c(1, 2)), 'single non-NA number')
  expect_error(mm_fill_gaps_2s(dat, max_gap_hours=NA), 'single non-NA number')
  expect_error(
    mm_fill_gaps_2s(dat, max_gap_hours=mm_max_gap_hours_cap + 0.5),
    'must be <= 2 hours')

  for(bad in list(0, 1.5, NA, c(0.1, 0.2))) {
    expect_error(mm_fill_gaps_2s(dat, warn_fill_share=bad), 'above 0 and at most 1')
  }
})


# heavily filled days ------------------------------------------------------

test_that("a day just over the fill share warns, naming it; just under is silent", {
  # 08:00 onward, so every row is on the 2050-06-01 06:00 day. Two inserted
  # rows fill every column, and two more values are filled inside existing
  # rows of DO.obs.up: 4 of a full day's 96 values in that column, 4.2%
  dat <- make_gappy_data(n=40, start="2050-06-01 08:00:00")
  dat$DO.obs.up[20:21] <- NA
  gappy <- dat[-(10:11), ]

  expect_warning(
    suppressMessages(mm_fill_gaps_2s(gappy, warn_fill_share=0.04)),
    '2050-06-01')
  expect_warning(suppressMessages(mm_fill_gaps_2s(gappy, warn_fill_share=0.042)), NA)
})

test_that("filled values count toward the day of their own timestamp", {
  # a 4-bin gap from 05:30 to 06:15: split 2 and 2 across a 06:00 day start,
  # all 4 on one day with a midnight start
  gappy <- make_gappy_data(n=8, start="2050-06-01 05:00:00")[-(3:6), ]
  expect_warning(suppressMessages(mm_fill_gaps_2s(gappy, warn_fill_share=0.03)), NA)
  expect_warning(
    suppressMessages(mm_fill_gaps_2s(gappy, warn_fill_share=0.03, day_start_hour=0)),
    '2050-06-01')

  # the fill-share table groups days the same way
  split <- attr(suppressMessages(mm_fill_gaps_2s(gappy)), 'fill_share')
  expect_equal(plain_fill_share(split),
               data.frame(date=as.Date(c('2050-05-31', '2050-06-01')), share=c(2, 2)/96))
  whole <- attr(suppressMessages(mm_fill_gaps_2s(gappy, day_start_hour=0)), 'fill_share')
  expect_equal(plain_fill_share(whole), data.frame(date=as.Date('2050-06-01'), share=4/96))
  expect_identical(attr(whole, 'day_start_hour'), 0)
})

test_that("with no day_start_hour given, the warning groups days by the hour stored on the data", {
  # the same gap as above: it splits across a 06:00 day start, but not across
  # the stored midnight start, so only grouping by the stored hour warns
  gappy <- make_gappy_data(n=8, start="2050-06-01 05:00:00")[-(3:6), ]
  attr(gappy, 'day_start_hour') <- 0
  expect_warning(suppressMessages(mm_fill_gaps_2s(gappy, warn_fill_share=0.03)), '2050-06-01')
})

test_that("inserted rows dropped again do not count toward the filled share", {
  # the dropped-row fixture from the invariant test above, plus one genuine
  # fill: 1 of 96 values (1.0%) is under the threshold, but counting the 2
  # dropped rows too would make it 3 of 96 (3.1%)
  dat <- make_gappy_data(n=40, start="2050-06-01 08:00:00")
  dat$DO.sat.up[15:26] <- NA
  dat$DO.obs.down[35] <- NA
  gappy <- dat[-(20:21), ]
  expect_warning(suppressMessages(mm_fill_gaps_2s(gappy, warn_fill_share=0.015)), NA)
  fill_share <- attr(suppressMessages(mm_fill_gaps_2s(gappy)), 'fill_share')
  expect_equal(fill_share$share, 1/96)
})

test_that("the fill-share table has a row for every day, with the share from the most-filled column", {
  # the fixture above, extended past 06:00 into an ungapped second day.
  # DO.obs.up has 4 of 96 values filled on the first day, every other column 2
  dat <- make_gappy_data(n=120, start="2050-06-01 08:00:00")
  dat$DO.obs.up[20:21] <- NA
  gappy <- dat[-(10:11), ]

  fill_share <- attr(suppressMessages(mm_fill_gaps_2s(gappy, warn_fill_share=1)), 'fill_share')
  expect_s3_class(fill_share, 'tbl_df')
  expect_equal(names(fill_share), c('date', 'share'))
  expect_s3_class(fill_share$date, 'Date')
  expect_equal(fill_share$date, as.Date(c('2050-06-01', '2050-06-02')))
  expect_equal(fill_share$share, c(4/96, 0))
  expect_identical(attr(fill_share, 'day_start_hour'), mm_day_start_2s)
})

test_that("the warning names exactly the days the fill-share table puts above the threshold", {
  # three whole days: 4 rows filled on the first (4.2%), 2 on the second
  # (2.1%), none on the third
  dat <- make_gappy_data(n=3 * 96, start="2050-06-01 06:00:00")
  gappy <- dat[-c(10:13, 110:111), ]

  warned <- character(0)
  filled <- withCallingHandlers(
    suppressMessages(mm_fill_gaps_2s(gappy, warn_fill_share=0.03)),
    warning=function(w) {
      warned <<- c(warned, conditionMessage(w))
      invokeRestart('muffleWarning')
    })
  fill_share <- attr(filled, 'fill_share')
  above <- fill_share[fill_share$share > 0.03, ]

  expect_length(warned, 1)
  expect_equal(regmatches(warned, gregexpr('[0-9]{4}-[0-9]{2}-[0-9]{2}', warned))[[1]],
               as.character(above$date))
  expect_match(warned, sprintf('%s (%.1f%%)', above$date, 100 * above$share), fixed=TRUE)
})

test_that("a second call adds its shares to the carried ones", {
  # a 2-bin gap fillable at a half-hour tolerance and a 4-bin gap fillable
  # only at one hour, both on the same day
  dat <- make_gappy_data(n=60, start="2050-06-01 08:00:00")
  gappy <- dat[-c(10:11, 30:33), ]

  first <- suppressMessages(mm_fill_gaps_2s(gappy, max_gap_hours=0.5))
  expect_equal(attr(first, 'fill_share')$share, 2/96)
  second <- suppressMessages(mm_fill_gaps_2s(first, max_gap_hours=1))
  expect_equal(attr(second, 'fill_share')$share, 6/96)
})

test_that("a carried fill-share table made with a different day_start_hour is an error", {
  gappy <- make_gappy_data(n=40)[-(10:11), ]
  filled <- suppressMessages(mm_fill_gaps_2s(gappy, day_start_hour=0))
  expect_error(
    mm_fill_gaps_2s(filled),
    'fill_share made with day_start_hour=0 but are being filled with day_start_hour=6')
})

test_that("a malformed carried fill-share table is an error", {
  filled <- suppressMessages(mm_fill_gaps_2s(make_gappy_data(n=40)[-(10:11), ]))
  fill_share <- attr(filled, 'fill_share')
  malformed <- function(value) {
    attr(filled, 'fill_share') <- value
    filled
  }
  msg <- "malformed 'fill_share'"
  expect_error(mm_fill_gaps_2s(malformed(fill_share$share)), msg)
  expect_error(mm_fill_gaps_2s(malformed(fill_share['date'])), msg)
  no_date <- fill_share
  no_date$date <- as.character(no_date$date)
  expect_error(mm_fill_gaps_2s(malformed(no_date)), msg)
  no_hour <- fill_share
  attr(no_hour, 'day_start_hour') <- NULL
  expect_error(mm_fill_gaps_2s(malformed(no_hour)), msg)

  with_share <- function(share) {
    value <- fill_share
    value$share <- share
    value
  }
  expect_error(mm_fill_gaps_2s(malformed(with_share(NA_real_))), msg)
  expect_error(mm_fill_gaps_2s(malformed(with_share(-0.1))), msg)
  expect_error(mm_fill_gaps_2s(malformed(with_share(1.1))), msg)
  repeated <- fill_share[c(1, 1), ]
  attr(repeated, 'day_start_hour') <- attr(fill_share, 'day_start_hour')
  expect_error(mm_fill_gaps_2s(malformed(repeated)), msg)

  with_hour <- function(hour) {
    value <- fill_share
    attr(value, 'day_start_hour') <- hour
    value
  }
  expect_error(mm_fill_gaps_2s(malformed(with_hour('6'))), 'day_start_hour must be a single non-NA number')
  expect_error(mm_fill_gaps_2s(malformed(with_hour(c(6, 6)))), 'day_start_hour must be a single non-NA number')
  expect_error(mm_fill_gaps_2s(malformed(with_hour(24))), 'day_start_hour must be >= 0 and < 24')
})

test_that("shares added across calls are capped at 1", {
  # one whole day, with every row but every ninth missing in one column: each
  # of ten 8-bin gaps is fillable at the 2-hour tolerance, 80 of 96 values.
  # the first call fills DO.obs.up that way, the second DO.obs.down
  dat <- make_gappy_data(n=96, start="2050-06-01 06:00:00")
  sparse <- seq_len(96) %% 9 != 1 & seq_len(96) < 91
  dat$DO.obs.up[sparse] <- NA
  first <- suppressWarnings(suppressMessages(mm_fill_gaps_2s(dat, max_gap_hours=2)))
  first$DO.obs.down[sparse] <- NA
  second <- suppressWarnings(suppressMessages(mm_fill_gaps_2s(first, max_gap_hours=2)))

  alone <- suppressWarnings(suppressMessages(mm_fill_gaps_2s(drop_fill_share(first), max_gap_hours=2)))
  expect_gt(attr(first, 'fill_share')$share + attr(alone, 'fill_share')$share, 1)
  expect_identical(attr(second, 'fill_share')$share, 1)
})

test_that("carried days no longer in the data drop out of the fill-share table", {
  dat <- make_gappy_data(n=120, start="2050-06-01 08:00:00")
  filled <- suppressMessages(mm_fill_gaps_2s(dat[-(10:11), ]))
  expect_equal(nrow(attr(filled, 'fill_share')), 2)

  first_day <- filled[filled$solar.time < as.POSIXct('2050-06-02 06:00:00', tz='UTC'), ]
  refilled <- suppressMessages(mm_fill_gaps_2s(first_day))
  expect_equal(plain_fill_share(attr(refilled, 'fill_share')),
               data.frame(date=as.Date('2050-06-01'), share=2/96))
})

test_that("the fill-share table is not carried past alignment", {
  dat <- make_full_days(n_days=2)
  filled <- suppressMessages(mm_fill_gaps_2s(dat[-(30:31), ]))
  expect_false(is.null(attr(filled, 'fill_share')))
  expect_null(attr(suppressMessages(mm_align_data_2s(filled)), 'fill_share'))
})

test_that("a day_start_hour that differs from the recorded one is an error", {
  dat <- make_gappy_data(n=20)
  attr(dat, 'day_start_hour') <- 6
  expect_error(
    mm_fill_gaps_2s(dat, day_start_hour=0),
    'day_start_hour=6.*day_start_hour=0')
})


# integration with the day-window check -------------------------------------

test_that("a day dropped for a short gap is recovered by filling", {
  dat <- make_full_days(n_days=2)
  dates <- mm_date_2s(dat$solar.time, mm_day_start_2s)
  target <- unique(dates)[2] # the first full 06:00-06:00 day

  # knock a 3-bin hole in the middle of that day
  hole <- which(dates == target)[40:42]
  gappy <- dat[-hole, ]

  alignment_before <- suppressMessages(mm_align_2s(gappy, day_start_hour=mm_day_start_2s))
  expect_false(target %in% alignment_before$date)

  filled <- suppressMessages(mm_fill_gaps_2s(gappy))
  alignment_after <- suppressMessages(mm_align_2s(filled, day_start_hour=mm_day_start_2s))
  expect_true(target %in% alignment_after$date)
  expect_equal(sum(alignment_after$date == target), 96)
})

test_that("a gap straddling 06:00 recovers both of the days it touches", {
  dat <- make_full_days(n_days=3)
  dates <- mm_date_2s(dat$solar.time, mm_day_start_2s)
  boundary <- min(which(dates == unique(dates)[3])) # first row of the 2nd full day

  # two rows before the 06:00 boundary and two after it
  gappy <- dat[-((boundary - 2):(boundary + 1)), ]

  alignment_before <- suppressMessages(mm_align_2s(gappy, day_start_hour=mm_day_start_2s))
  expect_false(any(unique(dates)[2:3] %in% alignment_before$date))

  filled <- suppressMessages(mm_fill_gaps_2s(gappy))
  alignment_after <- suppressMessages(mm_align_2s(filled, day_start_hour=mm_day_start_2s))
  expect_true(all(unique(dates)[2:3] %in% alignment_after$date))
})

test_that("an over-tolerance gap still costs its day", {
  dat <- make_full_days(n_days=2)
  dates <- mm_date_2s(dat$solar.time, mm_day_start_2s)
  target <- unique(dates)[2]

  hole <- which(dates == target)[40:45] # 6 bins, past the tolerance
  gappy <- dat[-hole, ]

  filled <- suppressMessages(mm_fill_gaps_2s(gappy))
  expect_equal(nrow(filled), nrow(gappy))
  alignment <- suppressMessages(mm_align_2s(filled, day_start_hour=mm_day_start_2s))
  expect_false(target %in% alignment$date)
})


# light: the na.rm bias this fixes ------------------------------------------

test_that("filling raw light before normalizing removes mm_lag_light_2s's na.rm bias", {
  # mm_lag_light_2s() sums with na.rm=TRUE in both the travel-time window and
  # the day total, which silently treats a missing light reading as darkness
  # rather than as missing. Filling the raw series first makes na.rm a no-op.
  # This is an intentional behavior change: the proportions below differ
  # before and after, and the "after" values are the correct ones.
  n <- 26
  solar.time <- as.POSIXct("2050-06-01 04:00:00", tz="UTC") +
    as.difftime((0:(n - 1)) * 1, units="hours")
  travel.time <- rep(2/24, n) # 2 hours -> lag=2, window width 3 rows
  light <- rep(12, n)

  complete <- mm_lag_light_2s(solar.time, light, travel.time, day_start_hour=mm_day_start_2s)

  # one missing reading inside the complete 06:00-06:00 day
  gappy <- light
  gappy[10] <- NA
  biased <- mm_lag_light_2s(solar.time, gappy, travel.time, day_start_hour=mm_day_start_2s)

  # na.rm=TRUE drops the missing term from the day total (24 rows x 12 = 288
  # becomes 276), so every row of the day is divided by too small a number
  expect_equal(unname(biased[3]), 36 / 276)
  expect_equal(unname(complete[3]), 36 / 288)
  expect_gt(biased[3], complete[3]) # inflated, not merely different

  # filling the raw series restores the correct denominator everywhere
  dat <- data.frame(
    solar.time=solar.time, light=gappy, travel.time=travel.time)
  filled <- suppressMessages(mm_fill_gaps_2s(dat, max_gap_hours=2))
  repaired <- mm_lag_light_2s(filled$solar.time, filled$light, filled$travel.time, day_start_hour=mm_day_start_2s)
  expect_equal(repaired, complete)
})

test_that("mm_format_data_2s fills gaps before converting light to a proportion", {
  n <- 200
  stamps <- as.POSIXct("2050-06-01 00:00:00", tz="UTC") +
    as.difftime((seq_len(n) - 1) * 15, units="mins")
  up <- data.frame(timestamp=stamps, DO.obs=9 + seq_len(n)/100, DO.sat=rep(10, n))
  down <- data.frame(
    timestamp=stamps, DO.obs=8.8 + seq_len(n)/100, DO.sat=rep(9.9, n),
    temp.water=rep(20, n), depth=rep(0.5, n), travel.time=rep(0.01, n))
  lit <- data.frame(timestamp=stamps, light=100 + seq_len(n))

  full <- suppressMessages(mm_format_data_2s(up, down, lit))

  # remove 3 timesteps from every input, as a real sensor outage would
  hole <- 100:102
  gappy <- suppressMessages(mm_format_data_2s(
    up[-hole, ], down[-hole, ], lit[-hole, ]))

  expect_equal(nrow(gappy), nrow(full))
  expect_equal(unitted::v(gappy$solar.time), unitted::v(full$solar.time))
  # light was filled while still raw, so the proportion matches the
  # never-gappy result rather than being computed over a short day total
  expect_equal(unitted::v(gappy$light), unitted::v(full$light))
})


# real data -----------------------------------------------------------------

test_that("filling behaves on the real two_station_raw_example", {
  data(two_station_raw_example)
  raw <- two_station_raw_example
  expect_warning(
    dat <- suppressMessages(mm_format_data_2s(raw$upstream, raw$downstream, raw$light)),
    '2011-01-12')
  # the threshold reaches the fill step through mm_format_data_2s()
  expect_warning(
    suppressMessages(mm_format_data_2s(raw$upstream, raw$downstream, raw$light, warn_fill_share=1)),
    NA)

  # every timestamp is still on the bin grid, still sorted, still unique.
  # "on the grid" allows the same 1-second slack mm_lag_2s() allows, since
  # timestamp arithmetic leaves sub-second noise on real data
  stamps <- unitted::v(dat$solar.time)
  expect_false(is.unsorted(stamps))
  expect_false(any(duplicated(stamps)))
  expect_lt(max(abs(as.numeric(stamps) - as.numeric(mm_snap_to_bin_2s(stamps)))), 1)

  # filling is what lets more days through the completeness check
  unfilled <- suppressMessages(mm_format_data_2s(
    raw$upstream, raw$downstream, raw$light, max_gap_hours=15/60, warn_fill_share=1))
  expect_gt(nrow(dat), nrow(unfilled))

  n_days <- function(d) suppressMessages(mm_align_2s(unitted::v(d), day_start_hour=mm_day_start_2s))$n_days
  expect_gt(n_days(dat), n_days(unfilled))
})

test_that("filling formatted data again drops the light-window record, so alignment skips the check", {
  raw <- as_raw_2s(make_2day_lit_boundary_data())
  out <- suppressMessages(mm_format_data_2s(raw$upstream, raw$downstream, raw$light))
  expect_false(is.null(attr(out, 'light_check_failed')))

  refilled <- suppressMessages(mm_fill_gaps_2s(out))
  expect_null(attr(refilled, 'light_check_failed'))
  expect_identical(attr(suppressMessages(mm_align_data_2s(refilled)), 'light_check'), FALSE)

  # also when the tolerance is finer than the timestep and nothing can be filled
  unfillable <- suppressMessages(mm_fill_gaps_2s(out, max_gap_hours=1/60))
  expect_null(attr(unfillable, 'light_check_failed'))
})
