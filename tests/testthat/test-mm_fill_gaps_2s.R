# Build a two-station data.frame on a 15-minute grid. At that timestep the
# default 1-hour tolerance is exactly 4 bins, so a 4-bin gap is the largest
# fillable one and a 5-bin gap is the smallest unfillable one -- every
# boundary test below is written against those two numbers.
make_gappy_data <- function(n=40, start="2050-06-01 00:00:00", travel_time=0.01) {
  data.frame(
    solar.time = as.POSIXct(start, tz="UTC") +
      as.difftime((seq_len(n) - 1) * 15, units="mins"),
    DO.obs.up = 9 + seq_len(n) / 100,
    DO.sat.up = rep(10, n),
    DO.obs.down = 8.8 + seq_len(n) / 100,
    DO.sat.down = rep(9.9, n),
    light = 100 + seq_len(n),
    depth = rep(0.5, n),
    temp.water = rep(20, n),
    travel.time = rep(travel_time, n)
  )
}

# Whole 06:00-06:00 days at a 15-minute timestep: 96 rows per day, plus
# n_leadin rows before the first 06:00 so the modeled rows have upstream data
# to lag from.
make_full_days <- function(n_days=2, n_leadin=1, start_date="2050-06-01", travel_time=0.01) {
  day_start <- as.POSIXct(paste0(start_date, " 06:00:00"), tz="UTC")
  first <- day_start - as.difftime(n_leadin * 15, units="mins")
  n <- n_leadin + 96 * n_days
  dat <- make_gappy_data(n=n, travel_time=travel_time)
  dat$solar.time <- first + as.difftime((seq_len(n) - 1) * 15, units="mins")
  dat
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
  expect_equal(filled, gappy)
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
  expect_equal(filled, dat)
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

  filled <- suppressMessages(mm_fill_gaps_2s(down))

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
  expect_error(mm_fill_gaps_2s(offgrid), 'not on a snap-to-bin grid')

  duped <- rbind(dat, dat[10, ])
  duped <- duped[order(duped$solar.time), ]
  expect_error(mm_fill_gaps_2s(duped), 'more than one row in the same timestep bin')

  unsorted <- dat[c(2, 1, 3:20), ]
  expect_error(mm_fill_gaps_2s(unsorted), 'sorted ascending')
})

test_that("max_gap_hours is validated against the cap", {
  dat <- make_gappy_data(n=20)
  expect_error(mm_fill_gaps_2s(dat, max_gap_hours=0), 'must be > 0')
  expect_error(mm_fill_gaps_2s(dat, max_gap_hours=c(1, 2)), 'single non-NA number')
  expect_error(mm_fill_gaps_2s(dat, max_gap_hours=NA), 'single non-NA number')
  expect_error(
    mm_fill_gaps_2s(dat, max_gap_hours=mm_max_gap_hours_cap + 0.5),
    'must be <= 2 hours')
  expect_error(specs(mm_name('bayes_2s'), max_gap_hours=6), 'must be <= 2 hours')
})


# integration with the day-window check -------------------------------------

test_that("a day dropped for a short gap is recovered by filling", {
  dat <- make_full_days(n_days=2)
  dates <- mm_date_2s(dat$solar.time)
  target <- unique(dates)[2] # the first full 06:00-06:00 day

  # knock a 3-bin hole in the middle of that day
  hole <- which(dates == target)[40:42]
  gappy <- dat[-hole, ]

  aln_before <- suppressMessages(mm_align_2s(gappy))
  expect_false(target %in% aln_before$date)

  filled <- suppressMessages(mm_fill_gaps_2s(gappy))
  aln_after <- suppressMessages(mm_align_2s(filled))
  expect_true(target %in% aln_after$date)
  expect_equal(sum(aln_after$date == target), 96)
})

test_that("a gap straddling 06:00 recovers both of the days it touches", {
  dat <- make_full_days(n_days=3)
  dates <- mm_date_2s(dat$solar.time)
  boundary <- min(which(dates == unique(dates)[3])) # first row of the 2nd full day

  # two rows before the 06:00 boundary and two after it
  gappy <- dat[-((boundary - 2):(boundary + 1)), ]

  aln_before <- suppressMessages(mm_align_2s(gappy))
  expect_false(any(unique(dates)[2:3] %in% aln_before$date))

  filled <- suppressMessages(mm_fill_gaps_2s(gappy))
  aln_after <- suppressMessages(mm_align_2s(filled))
  expect_true(all(unique(dates)[2:3] %in% aln_after$date))
})

test_that("an over-tolerance gap still costs its day", {
  dat <- make_full_days(n_days=2)
  dates <- mm_date_2s(dat$solar.time)
  target <- unique(dates)[2]

  hole <- which(dates == target)[40:45] # 6 bins, past the tolerance
  gappy <- dat[-hole, ]

  filled <- suppressMessages(mm_fill_gaps_2s(gappy))
  expect_equal(nrow(filled), nrow(gappy))
  aln <- suppressMessages(mm_align_2s(filled))
  expect_false(target %in% aln$date)
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

  complete <- mm_lag_light_2s(solar.time, light, travel.time)

  # one missing reading inside the complete 06:00-06:00 day
  gappy <- light
  gappy[10] <- NA
  biased <- mm_lag_light_2s(solar.time, gappy, travel.time)

  # na.rm=TRUE drops the missing term from the day total (24 rows x 12 = 288
  # becomes 276), so every row of the day is divided by too small a number
  expect_equal(unname(biased[3]), 36 / 276)
  expect_equal(unname(complete[3]), 36 / 288)
  expect_gt(biased[3], complete[3]) # inflated, not merely different

  # filling the raw series restores the correct denominator everywhere
  dat <- data.frame(
    solar.time=solar.time, light=gappy, travel.time=travel.time)
  filled <- suppressMessages(mm_fill_gaps_2s(dat, max_gap_hours=2))
  repaired <- mm_lag_light_2s(filled$solar.time, filled$light, filled$travel.time)
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
  dat <- suppressMessages(mm_format_data_2s(raw$upstream, raw$downstream, raw$light))

  # every timestamp is still on the bin grid, still sorted, still unique.
  # "on the grid" allows the same 1-second slack mm_lag_2s() allows, since
  # timestamp arithmetic leaves sub-second noise on real data
  stamps <- unitted::v(dat$solar.time)
  expect_false(is.unsorted(stamps))
  expect_false(any(duplicated(stamps)))
  expect_lt(max(abs(as.numeric(stamps) - as.numeric(mm_snap_to_bin_2s(stamps)))), 1)

  # filling is what lets more days through the completeness check
  unfilled <- suppressMessages(mm_format_data_2s(
    raw$upstream, raw$downstream, raw$light, max_gap_hours=15/60))
  expect_gt(nrow(dat), nrow(unfilled))

  n_days <- function(d) suppressMessages(mm_align_2s(unitted::v(d)))$n_days
  expect_gt(n_days(dat), n_days(unfilled))
})
