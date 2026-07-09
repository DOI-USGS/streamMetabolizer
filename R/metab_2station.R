#' Two-station Bayesian metabolism model fitting function (stub)
#'
#' Fits a two-station (upstream/downstream, a.k.a. VFTS) Bayesian model to
#' estimate GPP and ER from paired upstream and downstream DO, temperature,
#' light, and travel-time data. This function is currently a stub: it
#' validates the \code{data} argument and enforces two-station-specific data
#' requirements, but does not yet fit a model. See \code{\link{mm_name}} to
#' choose a Bayesian model and \code{\link{specs}} for relevant options for
#' the \code{specs} argument.
#'
#' @inheritParams metab
#' @return Not yet implemented; currently always errors after data validation.
#'
#' @section Two-station data requirements: In addition to the checks
#'   performed by \code{\link{mm_validate_data}}, \code{data$travel.time} (the
#'   reach travel time between the upstream and downstream stations, in days)
#'   must be strictly positive and less than 1 day (values >= 1 usually
#'   indicate travel time was supplied in the wrong units). There must also be
#'   enough lead-in observations of upstream DO before the first row of
#'   \code{data} to cover the longest travel time in the dataset, given the
#'   (median) timestep of \code{data$solar.time}.
#'
#' @export
#' @family metab_model
metab_2station <- function(
  specs=specs(mm_name('bayes_2station')),
  data=mm_data(solar.time, DO.obs.up, DO.sat.up, DO.obs.down, DO.sat.down,
               light, depth, temp.water, travel.time),
  data_daily=mm_data(date, optional='all'),
  info=NULL
) {

  # Check data for correct column names & units
  dat_list <- mm_validate_data(data, data_daily, 'metab_2station')

  travel_time <- v(dat_list$data$travel.time)
  solar_time <- v(dat_list$data$solar.time)

  # a. travel.time must be strictly positive
  if(any(travel_time <= 0)) {
    stop('travel.time must be > 0')
  }

  # b. travel.time is expected in days, so any value >= 1 almost certainly
  # reflects a units mistake (e.g., minutes or hours rather than days)
  if(any(travel_time >= 1)) {
    stop('travel.time must be < 1 day; values >= 1 suggest incorrect units (expected days)')
  }

  # c. there must be enough lead-in rows of upstream DO before the first
  # modeled row to cover the longest travel time in the dataset. timestep_days
  # is the median observation interval, in days; max_lag is the number of
  # timesteps by which upstream data must lead downstream predictions. The
  # first max_lag rows of data serve only as lead-in and cannot themselves be
  # modeled, so at least max_lag + 1 rows are required overall.
  timestep_days <- stats::median(as.numeric(diff(solar_time), units='days'))
  max_lag <- max(round(travel_time / timestep_days))
  if(nrow(dat_list$data) <= max_lag) {
    lead_in_needed <- max_lag - nrow(dat_list$data) + 1
    stop(paste0(
      'insufficient lead-in data for upstream DO: the longest travel.time implies a lag of ',
      max_lag, ' timestep(s), but only ', nrow(dat_list$data), ' row(s) were supplied; ',
      'need ', lead_in_needed, ' more lead-in timestep(s) of upstream data before the first modeled row'))
  }

  stop('metab_2station not yet implemented')
}
