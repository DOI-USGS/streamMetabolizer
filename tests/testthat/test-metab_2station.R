
# Build a minimal, valid two-station data.frame. Defaults give a 5-minute
# timestep (0.0034722 days) and a 0.01-day travel time, so
# max_lag = round(0.01 / 0.0034722) = 3 timesteps of required upstream lead-in.
make_2station_data <- function(n=10, timestep_min=5, travel_time=0.01) {
  data.frame(
    solar.time = as.POSIXct("2050-06-01 00:00:00", tz="UTC") +
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

test_that("mm_validate_data catches missing required columns", {
  dat <- dplyr::select(make_2station_data(), -DO.obs.up)
  expect_error(metab_2station(data=dat), "missing these columns")
})

test_that("travel.time <= 0 triggers an error", {
  dat <- make_2station_data(travel_time=0)
  expect_error(metab_2station(data=dat), "travel.time must be > 0")

  dat <- make_2station_data(travel_time=-0.01)
  expect_error(metab_2station(data=dat), "travel.time must be > 0")
})

test_that("travel.time >= 1 triggers an error with a units hint", {
  dat <- make_2station_data(travel_time=1)
  expect_error(metab_2station(data=dat), "travel.time must be < 1 day.*incorrect units")

  dat <- make_2station_data(travel_time=1.5)
  expect_error(metab_2station(data=dat), "travel.time must be < 1 day.*incorrect units")
})

test_that("insufficient lead-in data triggers an error", {
  # 2 rows but max_lag=3 timesteps of upstream lead-in are needed
  dat <- make_2station_data(n=2)
  expect_error(metab_2station(data=dat), "insufficient lead-in data")
})

test_that("a minimal valid data.frame passes all data-format checks", {
  dat <- make_2station_data(n=10)
  # metab_2station is a stub, so a fully valid data.frame should make it all
  # the way through the column/units/travel.time/lead-in checks and fail only
  # on the final "not yet implemented" stop -- not on any validation check.
  expect_error(metab_2station(data=dat), "metab_2station not yet implemented")
})
