# Two-station test fixtures, shared across test files.

# Build a minimal, valid two-station data.frame. Defaults give a 5-minute
# timestep (0.0034722 days) and a 0.01-day travel time, so
# max_lag = round(0.01 / 0.0034722) = 3 timesteps of required upstream lead-in.
# n=291 by default: 3 lead-in rows before 2050-06-01's 06:00 start, plus one
# genuinely complete 06:00-06:00 day (288 rows at 5-min resolution), so the
# default fixture is a real day under mm_align_2s()'s completeness check
# rather than a partial-day fragment; n can still be overridden down to a
# toy size (e.g. n=2) for tests that specifically want too little data.
make_2station_data <- function(n=291, timestep_min=5, travel_time=0.01) {
  data.frame(
    solar.time = as.POSIXct("2050-06-01 05:45:00", tz="UTC") +
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

# A day whose travel.time exceeds specs$max_travel_time_days (10/24-day
# default, 0.5-day cap) is not a dataset-wide error: mm_align_2s() drops
# just that day, with a message naming the date and travel time (see
# mm_lag_2s.R). Built from two make_2station_data() day-shaped blocks: day 1
# at the default travel_time=0.01 (well under the ceiling), day 2 reusing
# day 1's complete-day rows as a template, re-dated to follow immediately
# after, with travel.time raised to 15 hours, 0.625 days (above the ceiling). Day 2
# still has real upstream lead-in -- it draws on day 1's data -- so the
# ceiling, not lead-in availability, is what drops it.
make_2day_ceiling_data <- function() {
  day1 <- make_2station_data()
  day2 <- day1[-(1:3), ] # drop day 1's own lead-in rows; keep the 288-row complete-day block as a template
  day2$solar.time <- max(day1$solar.time) + as.difftime(seq_len(nrow(day2)) * 5, units="mins")
  day2$travel.time <- 15/24
  rbind(day1, day2)
}

# Build a two-day, unit-labeled data.frame with a known, traceable
# DO.obs.up/DO.sat.up series (sequential integers) so the shift can be
# checked by exact value, plus a leading lead-in block. Hourly timestep
# (0.0416667 days) and 3-hour (0.125-day) travel.time give
# max_lag = round(0.125 / 0.0416667) = 3 lead-in timesteps. n_leadin=3 rows
# precede day 1's 06:00 start (an incomplete, and thus dropped, partial day
# of their own); day 1 and day 2 are each a genuinely complete 06:00-06:00
# window (24 hourly rows), so both modeled days end up with n_obs = 24 rows.
make_ts_data <- function(n_leadin=3, n_day1=24, n_day2=24, travel_time=0.125, unitted=FALSE) {
  n_total <- n_leadin + n_day1 + n_day2
  solar.time <- as.POSIXct("2050-06-01 03:00:00", tz="UTC") +
    as.difftime((seq_len(n_total) - 1), units="hours")
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

# The first run of n_days consecutive two-station days that survive both
# alignment and the day-validity tests, together with the alignment they came
# from. two_station_example spans six years of real sensor record and carries
# 90 multi-day gaps, so its valid days are not all adjacent: choosing the
# first n_days valid days by position would drop a real gap inside a window
# the tests below describe as consecutive, and days bordering a gap don't
# behave like interior ones. Fails loudly rather than quietly returning a
# gap-spanning run.
first_valid_2station_days <- function(full_data, n_days) {
  aln <- suppressMessages(mm_align_2s(v(full_data)))
  aln <- suppressMessages(mm_filter_valid_days_2s(v(full_data), aln))$aln
  dates <- unique(aln$date)
  runs <- vapply(
    seq_len(max(0, length(dates) - n_days + 1)),
    function(i) all(diff(dates[i:(i + n_days - 1)]) == 1),
    logical(1))
  if(!any(runs)) {
    stop('no run of ', n_days, ' consecutive valid two-station days available')
  }
  start <- which(runs)[1]
  list(aln=aln, dates=dates[start:(start + n_days - 1)])
}

# Subset two_station_example to just a few modeled days for a faster test
# fit. Naively slicing rows doesn't work: max_lag (the number of upstream
# lead-in rows required) is recomputed from whatever travel.time values are
# present in the slice, so an arbitrary row range can leave a partial first
# date once prepdata_bayes_2s() trims max_lag rows off the front -- the same
# lead-in-sizing logic used in data-raw/two_station_example.R is needed here
# too. Unlike subset_2station_days(), the lead-in block sized here is the
# window-wide worst case rather than each row's own requirement, so part of
# it is itself modelable and surfaces as one extra valid_day=FALSE day.
subset_2station_data <- function(full_data, n_modeled_days) {
  solar_time <- v(full_data$solar.time)
  timestep_days <- stats::median(as.numeric(diff(solar_time), units='days'))

  modeled_dates <- first_valid_2station_days(full_data, n_modeled_days)$dates
  # two-station days run 06:00 -> 06:00 the next day (not calendar midnight
  # to midnight), so the last requested day's window isn't complete until one
  # timestep before the following day's 06:00
  modeled_start <- as.POSIXct(paste0(modeled_dates[1], ' 06:00:00'), tz='UTC')
  modeled_end <- as.POSIXct(paste0(modeled_dates[length(modeled_dates)] + 1, ' 06:00:00'), tz='UTC') -
    as.difftime(timestep_days, units='days')

  candidate_start <- modeled_start - as.difftime(1, units='days')
  candidate <- full_data[solar_time >= candidate_start & solar_time <= modeled_end, ]
  max_lag <- max(round(v(candidate$travel.time) / timestep_days))
  lead_in_start <- modeled_start - as.difftime(max_lag * timestep_days, units='days')

  full_data[solar_time >= lead_in_start & solar_time <= modeled_end, ]
}

# Subset two_station_example to its first n_days consecutive complete
# 06:00-06:00 days, keeping the leading rows those days need for upstream
# lead-in. Driven by the alignment itself rather than by calendar dates, so
# the subset can't disagree with the day partition the fitting code will
# recompute from it. Because the selected days are consecutive, the row range
# below spans only the lead-in block and those days -- with a gap-spanning
# selection it would also sweep in the partial days sitting in the gap.
subset_2station_days <- function(full_data, n_days) {
  sel <- first_valid_2station_days(full_data, n_days)
  rows <- which(sel$aln$date %in% sel$dates)
  full_data[min(sel$aln$shift_idx[rows]):max(sel$aln$keep[rows]), ]
}

fast_2station_specs <- function() {
  specs(mm_name('bayes_2s'), n_chains=1, n_cores=1,
        burnin_steps=100, saved_steps=100, verbose=FALSE)
}

# Corrupt one modeled row of the middle day of an n-day slice, so that
# day_tests drops exactly that day and the days on either side are untouched.
corrupt_middle_day <- function(n_days=3, col='depth', value=0) {
  dat <- subset_2station_days(two_station_example, n_days)
  aln <- suppressMessages(mm_align_2s(v(dat)))
  bad_date <- unique(aln$date)[2]
  dat[[col]][aln$keep[aln$date == bad_date][10]] <- u(value, get_units(dat[[col]]))
  list(data=dat, bad_date=bad_date, dates=unique(aln$date))
}

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
