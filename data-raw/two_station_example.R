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

# metab_bayes_2s()'s upstream-DO lag shift (see prepdata_bayes_2s() and
# metab_bayes_2s()'s "Two-station data requirements" section) needs
# max_lag = max(round(travel.time / timestep_days)) rows of lead-in
# immediately before modeled_start -- and because prepdata_bayes_2s() trims
# max_lag rows off the *start of the whole array*, not off each calendar
# day, that lead-in window must be exactly max_lag rows (not e.g. a whole
# extra day) or the first modeled date ends up with a different row count
# than the rest, which prepdata_bayes_2s() rejects. max_lag is computed from a
# generous 2-day candidate lead-in window and then trimmed to size.
candidate_start <- modeled_start - as.difftime(2, units='days')
candidate <- vfts2 %>% filter(datetime >= candidate_start, datetime <= modeled_end)
max_lag <- max(round(candidate$travel_time / timestep_days))
lead_in_start <- modeled_start - as.difftime(max_lag * 15, units='mins')

vfts2_window <- vfts2 %>%
  filter(datetime >= lead_in_start, datetime <= modeled_end)

# confirm the window is gap-free at the native 15-min timestep, and that
# trimming the lead-in rows (as prepdata_bayes_2s() does) leaves exactly 30
# modeled dates with equal row counts
stopifnot(all(abs(diff(as.numeric(vfts2_window$datetime)) - 15*60) < 1e-6))
modeled_dates <- as.Date(vfts2_window$datetime[seq.int(max_lag+1, nrow(vfts2_window))])
stopifnot(length(unique(table(modeled_dates))) == 1)
stopifnot(length(unique(modeled_dates)) == 30)

# --- rename to package conventions & attach units -----------------------

renamed <- vfts2_window %>%
  transmute(
    solar.time = datetime,
    DO.obs.up = upstream_DO,
    DO.sat.up = upstream_DO_sat,
    DO.obs.down = downstream_DO,
    DO.sat.down = downstream_DO_sat,
    light = light,
    depth = reach_depth,
    temp.water = downstream_temp,
    travel.time = travel_time)
# model_run and lag are intentionally dropped (not package columns)

template <- mm_data(solar.time, DO.obs.up, DO.sat.up, DO.obs.down, DO.sat.down,
                     light, depth, temp.water, travel.time)
two_station_example <- renamed
for(col in names(template)) {
  two_station_example[[col]] <- u(renamed[[col]], get_units(template[[col]]))
}
two_station_example <- two_station_example[names(template)]

# --- save -----------------------------------------------------------------

usethis::use_data(two_station_example, overwrite=TRUE)
