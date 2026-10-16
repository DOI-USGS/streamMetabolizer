# Builds data/two_station_raw_example.rda from the raw, pre-formatting
# VFTS two-station source files: the upstream and downstream sonde exports
# and a separately-sourced modeled light series. Not run automatically as
# part of the package build/check; run manually (with the package root as
# the working directory) whenever the example dataset needs to be
# regenerated.
#
# Unlike two_station_example (already joined, lagged, and formatted for
# metab_bayes_2s()), this dataset stages the three source files as-is: columns
# renamed to package convention, artifact columns dropped, mixed types coerced,
# and the upstream/downstream deployment windows truncated to their shared
# overlap -- but NOT joined, lagged, or converted to mean solar time. That
# processing belongs to a later formatting step; see
# mm_ts_prep_data()/prepdata_bayes_2s() for the equivalent
# already-formatted-data logic).
#
# Source files (not included in the package; obtained separately):
#   2_station/Data/upstream_inputs.xlsx
#   2_station/Data/downstream_inputs.xlsx
#   2_station/Data/LF_twostation_yard_light_2008-2024.csv

library(dplyr)
library(readxl)
library(unitted)

data_dir <- file.path('..', '2_station', 'Data')

# --- read ------------------------------------------------------------------

upstream_xlsx <- read_excel(file.path(data_dir, 'upstream_inputs.xlsx'))
downstream_xlsx <- read_excel(file.path(data_dir, 'downstream_inputs.xlsx'))
light_csv <- read.csv(
  file.path(data_dir, 'LF_twostation_yard_light_2008-2024.csv'),
  stringsAsFactors=FALSE)

# --- rename to package conventions, drop export artifacts, coerce types ----

# oSat is exported as character in the upstream file (including two
# "NA" strings) but numeric in the downstream file; as.numeric() correctly
# maps both real numbers and the literal string "NA" to NA. ...1 (a row-index
# column) and decimal_days_since_deployment_start (which resets at each
# concatenated deployment) are export artifacts, not measurements.
# barometric_pressure is dropped too -- DO.sat is supplied directly and
# doesn't need to be recalculated from it.
upstream <- upstream_xlsx %>%
  transmute(
    timestamp = datetime,
    DO.obs = dissolved_oxygen,
    DO.sat = as.numeric(oSat))

downstream <- downstream_xlsx %>%
  transmute(
    timestamp = datetime,
    DO.obs = dissolved_oxygen,
    DO.sat = as.numeric(oSat),
    temp.water = temp,
    depth = depth,
    travel.time = tt)

# Direct + Indirect are summed into a single combined irradiance value per
# timestep
light <- light_csv %>%
  transmute(
    timestamp = as.POSIXct(DateTime, format='%Y-%m-%dT%H:%M:%SZ', tz='UTC'),
    light = Direct + Indirect)

# --- truncate to the upstream/downstream overlap window --------------------

# The upstream and downstream sondes were not deployed over the same window
# (upstream stops in Jan 2011; downstream continues to Feb 2014), and the
# light series is a modeled, gap-free set running all the way to 2025.
# Compare ranges up front and truncate all three data.frames to the shared
# upstream/downstream overlap in one step, with a single message, rather
# than joining/filtering row-by-row.
overlap_start <- max(min(upstream$timestamp), min(downstream$timestamp))
overlap_end <- min(max(upstream$timestamp), max(downstream$timestamp))

upstream <- filter(upstream, timestamp >= overlap_start, timestamp <= overlap_end)
downstream <- filter(downstream, timestamp >= overlap_start, timestamp <= overlap_end)
light <- filter(light, timestamp >= overlap_start, timestamp <= overlap_end)

# --- attach units ------------------------------------------------------

# timestamp is raw clock time, not yet converted to mean solar time (that
# conversion is a later formatting step), so it's tagged NA rather than
# reusing mm_data()'s 'solar.time' column.
attach_mm_units <- function(df, ...) {
  template <- mm_data(...)
  for(col in names(template)) {
    df[[col]] <- u(df[[col]], get_units(template[[col]]))
  }
  df$timestamp <- u(df$timestamp, NA)
  as.data.frame(df[c('timestamp', names(template))])
}

upstream <- attach_mm_units(upstream, DO.obs, DO.sat)
downstream <- attach_mm_units(downstream, DO.obs, DO.sat, temp.water, depth, travel.time)

# light here is raw per-timestep irradiance (Direct + Indirect summed), not
# yet converted to the day-normalized, travel-time-weighted proportion that
# two_station_example's 'light' column holds (see mm_lag_light_2s()) -- that
# conversion happens in a later formatting step, not here. The source file
# doesn't state units explicitly; tagged to match mm_data()'s 'light' units
# (umol m^-2 s^-1, photosynthetic photon flux density) based on the
# two-station metadata's description of the Yard et al. (2005) incident-light
# model this series derives from, and because its value range (up to ~2025)
# matches the metadata's documented range for that same PPFD quantity.
light <- attach_mm_units(light, light)

# --- assemble & save ---------------------------------------------------

two_station_raw_example <- list(
  upstream = upstream,
  downstream = downstream,
  light = light)

usethis::use_data(two_station_raw_example, overwrite=TRUE)
