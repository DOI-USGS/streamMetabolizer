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

# The '_workshoprevision' file is the authors' corrected export. It carries
# the same VFTS/one-station payload as the original, column for column, but
# names its timestamp column downstream_datetime (making explicit that the
# timestamp labels the DOWNSTREAM measurement) and stamps each row on a clock
# advanced 24 timesteps relative to the original file.
raw <- read.csv(
  file.path('..', '2_station', 'Data',
            '2_VFTS_and_One-station_model_input_workshoprevision.csv'),
  stringsAsFactors=FALSE)

# VFTS-2 is the variable-travel-time two-station run (as opposed to the
# fixed-travel-time VFTS-1/VFTS-3 runs or the one-station-only 'OS' run) and
# covers an intermediate date range (2008-03-11 to 2014-02-28)
vfts2 <- raw %>%
  filter(model_run == 'VFTS-2') %>%
  mutate(datetime = as.POSIXct(downstream_datetime, format='%Y-%m-%dT%H:%M:%SZ', tz='UTC')) %>%
  arrange(datetime)

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
# (datetime - lag) reconstructs the raw upstream series on the source grid.
# (VFTS-2's varying lag is NOT invertible this way -- where its lag decreases
# between rows, an upstream instant is referenced by no row at all and is
# simply absent from that run's column.)
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
  mutate(datetime = as.POSIXct(downstream_datetime, format='%Y-%m-%dT%H:%M:%SZ', tz='UTC')) %>%
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

# The full VFTS-2 range is staged deliberately -- no filtering down to a
# gap-free stretch. The series carries 90 interruptions, every one of them a
# whole day or longer (there are no sub-day gaps at all), and the days
# bordering them are dropped downstream by mm_align_2s()'s 06:00-06:00
# completeness rule rather than being papered over here. Real sensor records
# look like this; an example dataset that didn't would misrepresent what the
# model has to cope with.
renamed <- vfts2 %>%
  transmute(
    solar.time = datetime,
    DO.obs.down = downstream_DO,
    DO.sat.down = downstream_DO_sat,
    light = light,
    depth = reach_depth,
    temp.water = downstream_temp,
    travel.time = travel_time) %>%
  left_join(upstream_raw, by=c('solar.time'='datetime'))
# model_run, lag and downstream_discharge are intentionally dropped (not
# package columns)

# Rows whose upstream counterpart falls outside the reconstructed grid keep
# their NA rather than being dropped here. Removing them would punch extra
# holes in an otherwise regular 15-minute series; left in place, the day they
# belong to is caught by mm_filter_valid_days_2s()'s complete_data test and
# reported as an invalid day, which is the behavior this example should show.
stopifnot(sum(is.na(renamed$DO.obs.up)) == sum(is.na(renamed$DO.sat.up)))

# The travel.time shipped here is VFTS-2's (the variable-travel-time run this
# dataset is built from), NOT the fixed travel time VFTS-3 used to build the
# upstream column above. That is deliberate: travel.time drives the model's
# own lookback, and VFTS-2's varying values are what this example is meant to
# exercise. A consequence worth stating plainly, so it isn't mistaken later
# for a residual bug: the lag metab_bayes_2s() recomputes at fit time from
# these travel.time values disagrees with the CSV's own `lag` column on many
# rows. That is expected. The model's lookback is its own modeling choice;
# the CSV's `lag` column only records how the published file was built. The
# two need not agree -- what matters is that the upstream series these rows
# carry is now raw, so the model's lag is applied exactly once.

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
