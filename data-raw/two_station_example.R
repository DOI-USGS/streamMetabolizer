# Builds data/two_station_example.rda from the VFTS (Variable Flow
# Two-Station) paper's published input data. Not run automatically as part
# of the package build/check; run
# manually (with the package root as the working directory) whenever the
# example dataset needs to be regenerated.
#

# Download from ScienceBase: https://www.sciencebase.gov/catalog/item/6887d457d4be024722b4aae2

# Paper URL: https://aslopubs.onlinelibrary.wiley.com/doi/pdf/10.1002/lom3.70066

library(dplyr)
library(unitted)

# --- read & filter -----------------------------------------------------

raw <- read.csv(
  file.path('..', '2_station', 'Data', '2_VFTS_and_One-station_model_input.csv'),
  stringsAsFactors=FALSE)

# VFTS-2 is the variable-travel-time two-station run (as opposed to the
# fixed-travel-time VFTS-1/VFTS-3 runs or the one-station-only 'OS' run) and
# covers an intermediate date range (2008-03-11 to 2014-02-27)
vfts2 <- raw %>%
  filter(model_run == 'VFTS-2') %>%
  mutate(datetime = as.POSIXct(datetime, format='%Y-%m-%dT%H:%M:%SZ', tz='UTC')) %>%
  arrange(datetime)

# 30 consecutive days from the middle of the dataset (2011, roughly halfway
# between 2008 and 2014), avoiding the start/end edges of both the overall
# dataset and of this particular gap-free stretch of observations (verified
# separately to run gap-free from 2011-07-28 to 2011-08-30 at the source
# data's native 15-minute timestep).
modeled_start <- as.POSIXct('2011-07-31 00:00:00', tz='UTC')
modeled_end <- as.POSIXct('2011-08-29 23:45:00', tz='UTC')
timestep_days <- 15/(24*60)

# metab_bayes_2s() applies the upstream-DO lag shift at fit time: each modeled
# row reaches back round(travel.time / timestep_days) timesteps for its
# upstream observation. Rows whose shift would reach before the start of the
# data have no upstream match and are dropped per-row, so the dataset carries
# a lead-in block of max_lag rows immediately before modeled_start for the
# earliest modeled rows to reach into. max_lag is computed from a generous
# 2-day candidate lead-in window and then trimmed to size.
candidate_start <- modeled_start - as.difftime(2, units='days')
candidate <- vfts2 %>% filter(datetime >= candidate_start, datetime <= modeled_end)
max_lag <- max(round(candidate$travel_time / timestep_days))
lead_in_start <- modeled_start - as.difftime(max_lag * 15, units='mins')

vfts2_window <- vfts2 %>%
  filter(datetime >= lead_in_start, datetime <= modeled_end)

# confirm the window is gap-free at the native 15-min timestep and that the
# post-lead-in rows split evenly across calendar dates.
#
# NOTE: the 30 below counts the CALENDAR-DAY span of this raw input slice, not
# the number of two-station days that survive fitting. Two-station days run
# 06:00-06:00, so the calendar-date count here and the modeled-day count after
# alignment are different quantities: this slice spans 30 calendar dates but
# yields 29 modeled days, for the calendar-boundary reason tracked separately
# in the design doc. Don't read this 30 as "30 modeled days."
stopifnot(all(abs(diff(as.numeric(vfts2_window$datetime)) - 15*60) < 1e-6))
modeled_dates <- as.Date(vfts2_window$datetime[seq.int(max_lag+1, nrow(vfts2_window))])
stopifnot(length(unique(table(modeled_dates))) == 1)
stopifnot(length(unique(modeled_dates)) == 30)

# --- recover the genuinely contemporaneous upstream series ----------------
#
# The VFTS-2 run's upstream_DO/upstream_DO_sat columns are NOT raw: they are
# already shifted by that run's own per-row lag, so each row holds the
# upstream value from lag timesteps EARLIER. Using them directly would let
# metab_bayes_2s() apply its own lag on top, shifting upstream data twice.
#
# The VFTS-3 run supplies the fix. It is the same published file, the same
# reach and the same physical sensors -- its downstream_DO/downstream_DO_sat
# columns are bit-identical to VFTS-2's -- but it was run with a single fixed
# travel time, so its lag is a constant 25 timesteps across all rows rather
# than varying per row. A constant shift is exactly invertible: re-labelling
# each VFTS-3 row with the time its upstream value was actually measured
# (datetime - lag) reconstructs the raw upstream series on the full grid, with
# no gaps. (VFTS-2's varying lag is NOT invertible this way -- where its lag
# decreases between rows, an upstream instant is referenced by no row at all
# and is simply absent from that run's column.)
#
# Verified: this reconstruction reproduces VFTS-2's own pre-lagged column
# exactly (158,678 of 158,678 rows, zero deviation) once VFTS-2's per-row lag
# is re-applied to it, confirming both runs draw on one identical raw upstream
# series; and it matches the raw sonde export (upstream_inputs.xlsx) on
# 99.989% of rows over the earlier period that file covers, the residual being
# this CSV's 2-decimal rounding.
#
# So: upstream DO comes from the VFTS-3 rows, everything else from VFTS-2.
vfts3 <- raw %>%
  filter(model_run == 'VFTS-3') %>%
  mutate(datetime = as.POSIXct(datetime, format='%Y-%m-%dT%H:%M:%SZ', tz='UTC')) %>%
  arrange(datetime)

# the constant lag is the precondition for the whole approach -- assert it
# rather than assuming it, so a revised source file fails loudly here instead
# of silently mis-shifting the upstream series
stopifnot(!anyNA(vfts3$lag), length(unique(vfts3$lag)) == 1)
vfts3_lag <- unique(vfts3$lag)

upstream_raw <- data.frame(
  datetime = vfts3$datetime - as.difftime(vfts3_lag * 15, units='mins'),
  DO.obs.up = vfts3$upstream_DO,
  DO.sat.up = vfts3$upstream_DO_sat)

# --- rename to package conventions & attach units -----------------------

renamed <- vfts2_window %>%
  transmute(
    solar.time = datetime,
    DO.obs.down = downstream_DO,
    DO.sat.down = downstream_DO_sat,
    light = light,
    depth = reach_depth,
    temp.water = downstream_temp,
    travel.time = travel_time) %>%
  left_join(upstream_raw, by=c('solar.time'='datetime'))
# model_run and lag are intentionally dropped (not package columns)

# every modeled row must have a real upstream observation; a gap here would
# mean the VFTS-3 grid didn't cover this window
stopifnot(!anyNA(renamed$DO.obs.up), !anyNA(renamed$DO.sat.up))

# The travel.time shipped here is VFTS-2's (the variable-travel-time run this
# dataset is built from), NOT the fixed 0.26 d that VFTS-3 used to build the
# upstream column above. That is deliberate: travel.time drives the model's
# own lookback, and VFTS-2's varying values are what this example is meant to
# exercise. A consequence worth stating plainly, so it isn't mistaken later
# for a residual bug: the lag metab_bayes_2s() recomputes at fit time from
# these travel.time values disagrees with the CSV's own `lag` column on ~34%
# of window rows. That is expected. The model's lookback is its own modeling
# choice; the CSV's `lag` column only records how the published file was
# built. The two need not agree -- what matters is that the upstream series
# these rows carry is now raw, so the model's lag is applied exactly once.

template <- mm_data(solar.time, DO.obs.up, DO.sat.up, DO.obs.down, DO.sat.down,
                     light, depth, temp.water, travel.time)
two_station_example <- renamed
for(col in names(template)) {
  two_station_example[[col]] <- u(renamed[[col]], get_units(template[[col]]))
}
two_station_example <- two_station_example[names(template)]

# mm_data()'s shared 'light' template describes one-station's raw PAR
# (umol m^-2 s^-1), but VFTS-2's 'light' column is already the day-normalized,
# travel-time-weighted proportion that mm_lag_light_2s() computes (Bishop et
# al. 2026 Eq. 2) -- unitless, not PAR. Override to NA units, matching the
# package convention for dimensionless columns (e.g. err.obs.phi).
two_station_example$light <- u(v(two_station_example$light), NA)

# --- save -----------------------------------------------------------------

usethis::use_data(two_station_example, overwrite=TRUE)
