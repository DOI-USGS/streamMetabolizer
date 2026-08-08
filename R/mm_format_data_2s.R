#' @include mm_lag_2s.R mm_fill_gaps_2s.R
NULL

#' Combine raw two-station upstream/downstream/light data into
#' \code{metab_bayes_2s}'s expected shape
#'
#' Takes the three raw, unaligned two-station data.frames -- upstream DO,
#' downstream DO/temp/depth/travel time, and a single light series
#' (see \code{\link{two_station_raw_example}} for the expected input shape)
#' -- and produces the single data.frame that \code{\link{mm_validate_data}}
#' (including its two-station-specific lead-in check) and
#' \code{\link{metab_bayes_2s}} expect (see \code{\link{two_station_example}}).
#' Formatting data by hand to that same shape works equally well; this
#' function is a convenience for the common case of separate
#' upstream/downstream/light, not a requirement.
#'
#' "Raw" here means unaligned: inputs must already use
#' package-convention column names (\code{timestamp}, \code{DO.obs},
#' \code{DO.sat}, etc.).
#'
#' \code{downstream} is the central piece: every downstream row survives into the
#' output, with \code{NA} wherever no upstream or light observation exists at
#' that timestep. \code{upstream} and \code{light} are matched to it by
#' timestep bin.
#'
#' \code{upstream} and \code{downstream} are first trimmed to their shared
#' deployment window: rows outside it could never be matched to anything, and
#' comparing the two ranges up front gives one clear message about what was
#' dropped and why, rather than a scattered NA in every unmatched row.
#' \code{light} is trimmed to the same window (it is expected to be a
#' modeled, gap-free series that outlasts both sondes, per
#' \code{\link{two_station_raw_example}}).
#'
#' All three inputs' timestamps are snapped onto a common bin grid before
#' matching, so that sensor deployments concatenated at slightly different
#' clock-minute offsets don't read as unmatched.
#'
#' \code{light} is converted from a raw per-timestep value into the
#' within-day, travel-time-weighted proportion \code{\link{metab_bayes_2s}}
#' requires, by summing each row's upstream travel-time window and dividing
#' by that day's total light sum, and its units are set to \code{NA} to
#' match -- mirroring \code{\link{metab_bayes_2s}}'s own \code{data} default.
#'
#' @section Gap filling: brief interruptions in the record -- a sonde that
#'   logged nothing for a few timesteps, or logged \code{NA} -- are bridged by
#'   linear interpolation, so that a day marred by a short dropout can be
#'   modeled rather than discarded. Gaps longer than \code{max_gap_hours} are
#'   left in place, and \code{\link{metab_bayes_2s}} drops the days holding
#'   them. Filling happens here while \code{light} is still a raw
#'   per-timestep value, before the within-day proportion is computed, which
#'   is the accurate order: the proportion is then formed from a complete
#'   series divided by a complete day total. This is a reason to supply raw
#'   light and let this function convert it, rather than formatting data by
#'   hand -- \code{\link{metab_bayes_2s}} can still fill gaps in
#'   hand-formatted data, but by then \code{light} has already been
#'   normalized, and interpolating it is only an approximation.
#'
#' @param upstream data.frame or unitted data.frame with columns
#'   \code{timestamp} (POSIXct, UTC), \code{DO.obs}, \code{DO.sat}
#' @param downstream data.frame or unitted data.frame with columns
#'   \code{timestamp} (POSIXct, UTC), \code{DO.obs}, \code{DO.sat},
#'   \code{temp.water}, \code{depth}, \code{travel.time}
#' @param light data.frame or unitted data.frame with columns
#'   \code{timestamp} (POSIXct, UTC), \code{light} (a single combined value
#'   per timestep)
#' @param max_gap_hours the gap-filling tolerance, in hours: runs of missing
#'   data spanning no more than this are interpolated, longer runs are left in
#'   place. Defaults to 1 hour and may not be set above 2. To have gaps
#'   treated identically here and during fitting, pass the same value to
#'   \code{\link{specs}}'s \code{max_gap_hours}.
#' @return a unitted data.frame with columns \code{solar.time},
#'   \code{DO.obs.up}, \code{DO.sat.up}, \code{DO.obs.down},
#'   \code{DO.sat.down}, \code{light}, \code{depth}, \code{temp.water},
#'   \code{travel.time}, sorted ascending by \code{solar.time}, ready for
#'   \code{\link{mm_validate_data}}'s two-station-specific lead-in check
#' @importFrom unitted v u get_units
#' @export
mm_format_data_2s <- function(upstream, downstream, light,
                              max_gap_hours=mm_max_gap_hours_default) {

  mm_check_required_cols(upstream, c('timestamp','DO.obs','DO.sat'), 'upstream')
  mm_check_required_cols(downstream, c('timestamp','DO.obs','DO.sat','temp.water','depth','travel.time'), 'downstream')
  mm_check_required_cols(light, c('timestamp','light'), 'light')

  # inputs must already be the expected types -- malformed types are an
  # error here, not something this function fixes
  check_types <- function(df, arg_name) {
    df <- v(df)
    if(!lubridate::is.POSIXct(df$timestamp)) {
      stop(paste0(arg_name, '$timestamp must be of class POSIXct'), call.=FALSE)
    }
    if(!(lubridate::tz(df$timestamp) %in% c('UTC','GMT'))) {
      stop(paste0(arg_name, "$timestamp must have timezone 'UTC'"), call.=FALSE)
    }
    for(col in setdiff(names(df), 'timestamp')) {
      if(!is.numeric(df[[col]])) {
        stop(paste0(
          arg_name, '$', col, ' must be numeric, but is ', class(df[[col]])[1]), call.=FALSE)
      }
    }
    df
  }
  upstream <- check_types(upstream, 'upstream')
  downstream <- check_types(downstream, 'downstream')
  light <- check_types(light, 'light')

  # --- trim to the upstream/downstream shared deployment window -----------

  overlap_start <- max(min(upstream$timestamp), min(downstream$timestamp))
  overlap_end <- min(max(upstream$timestamp), max(downstream$timestamp))
  if(overlap_start > overlap_end) {
    stop(
      'upstream and downstream have no overlapping deployment window: upstream spans ',
      min(upstream$timestamp), ' to ', max(upstream$timestamp), ', downstream spans ',
      min(downstream$timestamp), ' to ', max(downstream$timestamp), call.=FALSE)
  }

  in_window <- function(df) df$timestamp >= overlap_start & df$timestamp <= overlap_end
  keep_upstream <- in_window(upstream)
  keep_downstream <- in_window(downstream)
  keep_light <- in_window(light)
  n_dropped <- sum(!keep_upstream) + sum(!keep_downstream) + sum(!keep_light)
  if(n_dropped > 0) {
    message(paste0(
      'trimming to the shared upstream/downstream deployment window (', overlap_start, ' to ',
      overlap_end, '): dropping ', sum(!keep_upstream), ' upstream row(s), ',
      sum(!keep_downstream), ' downstream row(s), ', sum(!keep_light), ' light row(s)'))
  }
  upstream <- upstream[keep_upstream, ]
  downstream <- downstream[keep_downstream, ]
  light <- light[keep_light, ]

  # --- snap all three onto a common bin grid, then match by timestep bin --

  downstream <- downstream[order(downstream$timestamp), ]
  upstream <- upstream[order(upstream$timestamp), ]
  light <- light[order(light$timestamp), ]
  downstream$timestamp <- mm_snap_to_bin_2s(downstream$timestamp)
  upstream$timestamp <- mm_snap_to_bin_2s(upstream$timestamp)
  light$timestamp <- mm_snap_to_bin_2s(light$timestamp)

  timestep_days <- unlist(mm_get_timestep(downstream$timestamp, format='modal'))
  if(length(timestep_days) != 1 || is.na(timestep_days)) {
    stop('could not detect a single nominal timestep for downstream$timestamp', call.=FALSE)
  }

  bin_down <- mm_bin_index_2s(downstream$timestamp, timestep_days)
  bin_up <- mm_bin_index_2s(upstream$timestamp, timestep_days)
  bin_light <- mm_bin_index_2s(light$timestamp, timestep_days)
  if(any(duplicated(bin_down))) {
    stop('downstream has more than one row in the same timestep bin; check for duplicate or near-duplicate timestamps', call.=FALSE)
  }
  if(any(duplicated(bin_up))) {
    stop('upstream has more than one row in the same timestep bin; check for duplicate or near-duplicate timestamps', call.=FALSE)
  }
  if(any(duplicated(bin_light))) {
    stop('light has more than one row in the same timestep bin; check for duplicate or near-duplicate timestamps', call.=FALSE)
  }

  shift_idx_up <- match(bin_down, bin_up)
  shift_idx_light <- match(bin_down, bin_light)

  # --- assemble onto downstream's row structure, no row ever dropped ------

  out <- data.frame(
    solar.time = downstream$timestamp,
    DO.obs.up = upstream$DO.obs[shift_idx_up],
    DO.sat.up = upstream$DO.sat[shift_idx_up],
    DO.obs.down = downstream$DO.obs,
    DO.sat.down = downstream$DO.sat,
    light = light$light[shift_idx_light],
    depth = downstream$depth,
    temp.water = downstream$temp.water,
    travel.time = downstream$travel.time)

  # --- fill short gaps, while light is still raw --------------------------
  #
  # order matters: the within-day proportion computed below divides each row's
  # travel-time window sum by that day's total light. Filling first means both
  # the numerator and the denominator are computed from a complete series;
  # interpolating the proportion afterwards could not repair the denominator,
  # which would still be short by the missing terms
  out <- mm_fill_gaps_2s(out, max_gap_hours=max_gap_hours)

  out$light <- mm_lag_light_2s(out$solar.time, out$light, out$travel.time)

  template <- mm_data(solar.time, DO.obs.up, DO.sat.up, DO.obs.down, DO.sat.down,
                       light, depth, temp.water, travel.time)
  for(col in names(template)) {
    out[[col]] <- u(out[[col]], get_units(template[[col]]))
  }
  # light here is the within-day proportion, not raw PAR -- see
  # metab_bayes_2s()'s data default, which overrides units the same way
  out$light <- u(v(out$light), NA)

  as.data.frame(out[names(template)])
}

#' Error clearly if a data.frame is missing any of a set of required columns
#'
#' @param df data.frame to check
#' @param required_cols character vector of column names \code{df} must have
#' @param arg_name used in the error message to identify which input failed
#' @keywords internal
#' @noRd
mm_check_required_cols <- function(df, required_cols, arg_name) {
  missing_cols <- setdiff(required_cols, names(df))
  if(length(missing_cols) > 0) {
    stop(paste0(arg_name, ' is missing these columns: ', paste(missing_cols, collapse=', ')), call.=FALSE)
  }
}
