#' Reshape long-format two-station data into the list expected by the
#' two-station Stan model
#'
#' Time-shifts the upstream DO series to match the travel time between
#' stations, then pivots the result into the \code{n_obs x n_days} matrices
#' expected by the \code{data} block of \code{inst/models/b2_np_oi_tr_plrckm.stan}
#' (see \code{\link{metab_2station}}).
#'
#' The upstream observation that "matches" a given downstream observation at
#' row \code{i} was recorded \code{lag[i] <- round(travel.time[i] /
#' timestep_days)} timesteps earlier, where \code{timestep_days} is the
#' median timestep of \code{data$solar.time}. This must be computed the same
#' way as in \code{\link{metab_2station}}'s lead-in check, so that the
#' \code{max(lag)} computed here always agrees with the lead-in requirement
#' already validated there. The first \code{max(lag)} rows of \code{data} are
#' lead-in rows: they supply upstream DO for the shift but are never
#' themselves treated as modeled (downstream) observations.
#'
#' @param data data.frame as validated by \code{\link{mm_validate_data}} for
#'   \code{\link{metab_2station}}: must contain \code{solar.time},
#'   \code{DO.obs.up}, \code{DO.sat.up}, \code{DO.obs.down},
#'   \code{DO.sat.down}, \code{light}, \code{depth}, \code{temp.water},
#'   \code{travel.time}, sorted ascending by \code{solar.time}, and must
#'   include the lead-in rows required to cover the longest travel time (see
#'   \code{\link{metab_2station}}).
#' @param specs optional list of model specs. If it contains
#'   \code{K600_lnorm_meanlog} and/or \code{K600_lnorm_sdlog}, those values
#'   are used; otherwise placeholder defaults of \code{log(3.48)} and
#'   \code{0.5} are used.
#' @return a named list with all variables in the Stan model's data block:
#'   \code{n_obs}, \code{n_days}, \code{DO_obs_up}, \code{DO_sat_up},
#'   \code{DO_obs_down}, \code{DO_sat_down}, \code{light}, \code{depth},
#'   \code{temp_water}, \code{travel_time} (each an \code{n_obs x n_days}
#'   matrix, unitless), and \code{K600_lnorm_meanlog}/\code{K600_lnorm_sdlog}
#' @importFrom unitted v
#' @export
mm_ts_prep_data <- function(data, specs=NULL) {

  # strip units; Stan cannot handle unitted vectors/matrices
  data <- v(data)

  # timestep_days must match metab_2station()'s lead-in check exactly (median
  # timestep, in days), so that max_lag here agrees with what was validated
  # there
  timestep_days <- stats::median(as.numeric(diff(data$solar.time), units='days'))

  # lag, in timesteps, that the upstream series must be shifted by to line up
  # with each row's downstream observation
  lag <- round(data$travel.time / timestep_days)
  max_lag <- max(lag)
  n_total <- nrow(data)
  if(n_total <= max_lag) {
    stop(
      'not enough lead-in rows to cover the longest travel.time (', max_lag, ' timesteps ',
      'implied, but only ', n_total, ' rows supplied); this should already have been caught ',
      'by metab_2station()')
  }

  # the first max_lag rows are lead-in only (upstream DO used for the shift,
  # but never modeled themselves). every row i > max_lag is guaranteed to
  # have a valid shift target (i - lag[i] >= 1) because lag[i] <= max_lag
  keep <- seq.int(max_lag + 1, n_total)
  shift_idx <- keep - lag[keep]
  if(any(shift_idx < 1)) {
    # should be unreachable given max_lag's definition; guards against
    # programming errors rather than expected user input
    stop('internal error: upstream shift index falls before the first row of data')
  }

  modeled <- data.frame(
    solar.time = data$solar.time[keep],
    DO_obs_up = data$DO.obs.up[shift_idx],
    DO_sat_up = data$DO.sat.up[shift_idx],
    DO_obs_down = data$DO.obs.down[keep],
    DO_sat_down = data$DO.sat.down[keep],
    light = data$light[keep],
    depth = data$depth[keep],
    temp_water = data$temp.water[keep],
    travel_time = data$travel.time[keep]
  )

  # pivot into n_obs x n_days matrices, one column per unique date, following
  # the same time_by_date_matrix approach used for the 1-station Stan models
  # (see prepdata_bayes() in metab_bayes.R)
  date_vec <- as.character(as.Date(modeled$solar.time))
  date_table <- table(date_vec)
  n_days <- length(date_table)
  n_obs_per_day <- unique(unname(date_table))
  if(length(n_obs_per_day) > 1) {
    stop(
      'dates have differing numbers of modeled rows after lead-in removal; ',
      'observations cannot be combined into a matrix: ',
      paste(sprintf('%s (%d rows)', names(date_table), date_table), collapse=', '))
  }
  n_obs <- n_obs_per_day

  to_matrix <- function(vec) matrix(vec, nrow=n_obs, ncol=n_days, byrow=FALSE)

  # confirm each date occupies a contiguous block of rows, i.e., that data
  # was sorted by solar.time; otherwise the matrix pivot below would silently
  # scramble which rows belong to which date
  date_mat <- to_matrix(date_vec)
  unique_dates_per_col <- apply(date_mat, MARGIN=2, FUN=unique)
  if(is.list(unique_dates_per_col) || !isTRUE(all.equal(unname(unique_dates_per_col), names(date_table)))) {
    stop('data must be sorted by solar.time so that each date occupies a contiguous block of rows')
  }

  K600_lnorm_meanlog <- if(!is.null(specs$K600_lnorm_meanlog)) specs$K600_lnorm_meanlog else log(3.48)
  K600_lnorm_sdlog <- if(!is.null(specs$K600_lnorm_sdlog)) specs$K600_lnorm_sdlog else 0.5

  list(
    n_obs = n_obs,
    n_days = n_days,
    DO_obs_up = to_matrix(modeled$DO_obs_up),
    DO_sat_up = to_matrix(modeled$DO_sat_up),
    DO_obs_down = to_matrix(modeled$DO_obs_down),
    DO_sat_down = to_matrix(modeled$DO_sat_down),
    light = to_matrix(modeled$light),
    depth = to_matrix(modeled$depth),
    temp_water = to_matrix(modeled$temp_water),
    travel_time = to_matrix(modeled$travel_time),
    K600_lnorm_meanlog = K600_lnorm_meanlog,
    K600_lnorm_sdlog = K600_lnorm_sdlog
  )
}
