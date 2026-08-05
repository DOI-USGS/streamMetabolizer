#' Default and maximum travel-time ceiling for two-station models
#'
#' The ceiling bounds how far back in time a modeled day may reach for its
#' upstream observations. Past the ceiling, the upstream parcel (and the light
#' it experienced) originates before the modeled day's own 06:00 start, so the
#' day-normalized light proportion no longer has a well-defined day total to be
#' a proportion of; such days are dropped rather than modeled.
#'
#' \code{mm_max_travel_time_cap} is deliberately well below the 24-hour day
#' window: a ceiling at the full window length would let a single day's travel
#' time consume the entire lead-in period, which defeats the purpose of having
#' a ceiling.
#'
#' \code{mm_max_travel_time_default} is the canonical default (currently 10).
#' \code{\link{specs}} cannot reference it directly: its \code{prefer_missing}
#' check flags any formal whose default is a bare symbol (via \code{is.symbol()}
#' on the unevaluated default) as "no real default," which would trigger a
#' false warning telling users to set \code{max_travel_time_hours} in
#' \code{\link{revise}()} instead -- exactly backwards, since overriding the
#' ceiling at \code{specs()}-creation time is the intended use. A
#' \code{NULL}-default workaround was considered and rejected: it avoids the
#' warning, but requires an extra \code{is.null()} resolution step inside
#' \code{specs()} plus a second guard at the \code{mm_align_2s()} call site in
#' \code{metab_bayes_2s.R} (duplicating one that already exists in
#' \code{prepdata_bayes_2s()}) -- trading one hand-sync problem for a second,
#' more fragile one. \code{specs()}'s literal default is kept in sync with
#' this constant by hand; a guard test in \code{test-metab_bayes_2s.R}
#' (search for "max_travel_time_hours default") fails if the two diverge.
#'
#' @keywords internal
#' @name mm_max_travel_time
NULL

#' @rdname mm_max_travel_time
mm_max_travel_time_default <- 10

#' @rdname mm_max_travel_time
mm_max_travel_time_cap <- 12

#' Hour at which a two-station day begins
#'
#' Two-station days run from 06:00 to 06:00 the next day, a full 24-hour
#' window, following Bishop et al. (2026). This is unrelated to one-station's
#' \code{day_start}/\code{day_end} (4 AM, 28 hours), which describe an
#' overlapping diel window rather than a partition of the time series; the two
#' conventions must not be conflated or shared.
#'
#' @keywords internal
mm_day_start_2s <- 6

#' Assign two-station day labels to timestamps
#'
#' Labels each timestamp with the two-station day it belongs to: the day
#' beginning at \code{mm_day_start_2s} (06:00) and running 24 hours. A
#' timestamp at 05:45 therefore belongs to the *previous* calendar date, and
#' one at 06:00 to the current calendar date.
#'
#' @param solar.time POSIXct vector of timestamps, in UTC (as enforced by
#'   \code{\link{mm_validate_data}}).
#' @return a Date vector of two-station day labels, the same length as
#'   \code{solar.time}
#' @keywords internal
mm_date_2s <- function(solar.time) {
  as.Date(solar.time - as.difftime(mm_day_start_2s, units='hours'), tz='UTC')
}

#' Compute the per-row upstream lag for two-station data
#'
#' The upstream observation that "matches" a given downstream observation at
#' row \code{i} was recorded \code{lag[i] <- round(travel.time[i] /
#' timestep_days)} timesteps earlier, where \code{timestep_days} is the median
#' timestep of \code{solar.time}. This is the single definition of that
#' calculation: \code{\link{mm_validate_data_2station}} (existence check) and
#' \code{mm_align_2s} (used in turn by \code{\link{prepdata_bayes_2s}} and
#' \code{\link{metab_bayes_2s}}) all route through it, so the three cannot
#' drift apart.
#'
#' No rows are dropped here; \code{has_leadin} reports which rows *could* be
#' modeled, and callers decide what to do about the rest.
#'
#' @param solar.time POSIXct vector of timestamps, sorted ascending.
#' @param travel.time numeric vector of reach travel times, in days, the same
#'   length as \code{solar.time}.
#' @return a list with \code{timestep_days} (median timestep, in days),
#'   \code{lag} (per-row lag, in timesteps), \code{shift_idx} (per-row index
#'   into the input rows supplying that row's upstream values; may be < 1 where
#'   the shift reaches before the start of the data), and \code{has_leadin}
#'   (logical, \code{shift_idx >= 1})
#' @keywords internal
mm_lag_2s <- function(solar.time, travel.time) {

  n_total <- length(solar.time)
  if(n_total < 2) {
    stop('need at least 2 rows of data to determine a timestep', call.=FALSE)
  }

  timestep_days <- stats::median(as.numeric(diff(solar.time), units='days'))
  lag <- round(travel.time / timestep_days)
  shift_idx <- seq_len(n_total) - lag

  list(
    timestep_days = timestep_days,
    lag = lag,
    shift_idx = shift_idx,
    has_leadin = shift_idx >= 1)
}

#' Align two-station data onto its 06:00 day window
#'
#' Applies, in order: the per-row lead-in test, two-station day labeling, the
#' travel-time ceiling, and the whole-day completeness requirement. The result
#' describes which input rows are modeled, which input rows supply their
#' upstream values, and how those modeled rows partition into days.
#'
#' Rows are dropped for lack of lead-in on a *per-row* basis: row \code{i} is
#' modelable whenever \code{i - lag[i] >= 1}. (An earlier implementation
#' trimmed \code{max(lag)} rows off the front of the dataset, discarding rows
#' that needed less lead-in than the dataset-wide worst case.)
#'
#' Whole days are dropped, with a message, in two cases:
#' \itemize{
#'   \item the day's longest travel time exceeds
#'     \code{max_travel_time_hours} (see \code{\link{mm_max_travel_time}});
#'   \item the day does not fill its 24-hour window, i.e. it holds fewer than
#'     \code{round(1 / timestep_days)} modeled rows. Partial days arise at the
#'     edges of any dataset whose bounds don't fall on 06:00, and from missing
#'     sensor readings. Dropping them is required by the fixed-shape
#'     \code{n_obs x n_days} matrices the Stan model expects; tolerating small
#'     within-day gaps by interpolation is separate, later work.
#' }
#'
#' @param data data.frame (units already stripped) containing at least
#'   \code{solar.time} and \code{travel.time}, sorted ascending by
#'   \code{solar.time}.
#' @param max_travel_time_hours the travel-time ceiling, in hours. Days whose
#'   longest travel time exceeds this are dropped.
#' @return a list with \code{keep} (integer indices of the modeled rows, into
#'   the rows of \code{data}), \code{shift_idx} (integer indices, also into
#'   \code{data}, of the rows supplying each modeled row's upstream values),
#'   \code{date} (Date vector of two-station day labels, parallel to
#'   \code{keep}), \code{n_obs} (modeled rows per day), \code{n_days}, and
#'   \code{timestep_days}
#' @keywords internal
mm_align_2s <- function(data, max_travel_time_hours=mm_max_travel_time_default) {

  solar_time <- data$solar.time
  travel_time <- data$travel.time

  lagged <- mm_lag_2s(solar_time, travel_time)

  # per-row lead-in: keep every row whose upstream shift lands within the data
  keep <- which(lagged$has_leadin)
  if(length(keep) == 0) {
    # mm_validate_data_2station() raises a more informative version of this
    # error before fitting begins; this guards direct calls
    stop('no rows have enough upstream lead-in data to be modeled', call.=FALSE)
  }
  date <- mm_date_2s(solar_time[keep])

  # travel-time ceiling: drop the offending day, not the whole dataset
  max_travel_time_days <- max_travel_time_hours / 24
  worst_by_day <- tapply(travel_time[keep], date, max)
  over_ceiling <- worst_by_day[worst_by_day > max_travel_time_days]
  if(length(over_ceiling) > 0) {
    message(paste0(
      'dropping ', length(over_ceiling), ' day(s) whose travel.time exceeds the ',
      max_travel_time_hours, '-hour ceiling: ',
      paste(sprintf('%s (%.2f hours)', names(over_ceiling), unname(over_ceiling) * 24), collapse=', ')))
    in_bounds <- !(as.character(date) %in% names(over_ceiling))
    keep <- keep[in_bounds]
    date <- date[in_bounds]
  }

  # completeness: every remaining day must fill its 24-hour window, because
  # the Stan model's n_obs x n_days matrices have no ragged/masking mechanism
  expected_n_obs <- round(1 / lagged$timestep_days)
  n_by_day <- table(date)
  partial <- n_by_day[n_by_day != expected_n_obs]
  if(length(partial) > 0) {
    message(paste0(
      'dropping ', length(partial), ' day(s) that do not fill the 06:00-06:00 window: ',
      paste(sprintf('%s (%d of %d expected observations)', names(partial), unname(partial), expected_n_obs),
            collapse=', ')))
    complete <- !(as.character(date) %in% names(partial))
    keep <- keep[complete]
    date <- date[complete]
  }

  if(length(keep) == 0) {
    stop(paste0(
      'no complete days remain after applying the ', max_travel_time_hours,
      '-hour travel.time ceiling and the 06:00-06:00 day-window requirement'), call.=FALSE)
  }

  list(
    keep = keep,
    shift_idx = lagged$shift_idx[keep],
    date = date,
    n_obs = expected_n_obs,
    n_days = length(unique(date)),
    timestep_days = lagged$timestep_days)
}
