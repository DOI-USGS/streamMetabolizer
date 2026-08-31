#' Default and maximum travel-time ceiling for two-station models
#'
#' Bounds how far back in time a modeled day may reach for upstream
#' observations; days whose travel time exceeds the ceiling are dropped
#' rather than modeled. \code{mm_max_travel_time_default} is the canonical
#' default (currently 10); \code{mm_max_travel_time_cap} is kept well below
#' the 24-hour day window so a single day's travel time can't consume the
#' entire lead-in period.
#'
#' @keywords internal
#' @name mm_max_travel_time
NULL

#' @rdname mm_max_travel_time
mm_max_travel_time_default <- 10

#' @rdname mm_max_travel_time
mm_max_travel_time_cap <- 12

#' Validate a travel-time ceiling
#'
#' Shared by \code{specs()} and \code{\link{mm_align_2s}} so the two can't
#' disagree about what ceiling is legal.
#'
#' @param max_travel_time_hours the value to check.
#' @keywords internal
mm_check_max_travel_time_hours <- function(max_travel_time_hours) {
  if(!is.numeric(max_travel_time_hours) || length(max_travel_time_hours) != 1 ||
     is.na(max_travel_time_hours)) {
    stop('max_travel_time_hours must be a single non-NA number', call.=FALSE)
  }
  if(max_travel_time_hours <= 0) {
    stop('max_travel_time_hours must be > 0', call.=FALSE)
  }
  if(max_travel_time_hours > mm_max_travel_time_cap) {
    stop(paste0(
      'max_travel_time_hours must be <= ', mm_max_travel_time_cap, ' hours; a ceiling ',
      'approaching the 24-hour two-station day window would allow a single day\'s travel ',
      'time to consume the entire lead-in period'), call.=FALSE)
  }
  invisible(NULL)
}

#' Hour at which a two-station day begins
#'
#' Two-station days run 06:00-06:00 (24 hours), following Bishop et al.
#' (2026) -- distinct from one-station's overlapping 4 AM-28-hour diel
#' window (\code{day_start}/\code{day_end}); the two must not be conflated.
#'
#' @keywords internal
mm_day_start_2s <- 6

#' Assign two-station day labels to timestamps
#'
#' Labels each timestamp with the two-station day it belongs to: the day
#' beginning at \code{mm_day_start_2s} (06:00) and running 24 hours. 05:45
#' belongs to the *previous* calendar date; 06:00 to the current one.
#'
#' @param solar.time POSIXct vector of timestamps, in UTC (as enforced by
#'   \code{\link{mm_validate_data}}).
#' @return a Date vector of two-station day labels, the same length as
#'   \code{solar.time}
#' @keywords internal
mm_date_2s <- function(solar.time) {
  as.Date(solar.time - as.difftime(mm_day_start_2s, units='hours'), tz='UTC')
}

#' Fixed epoch anchoring the two-station snap-to-bin grid
#'
#' Shared by \code{\link{mm_bin_index_2s}} (and, through it,
#' \code{\link{mm_snap_to_bin_2s}} and \code{\link{mm_lag_2s}}) so two series
#' snapped/binned independently -- e.g. upstream and downstream -- land on
#' the same absolute bins, rather than each anchored to its own first
#' timestamp.
#'
#' @keywords internal
mm_bin_epoch_2s <- as.POSIXct('1970-01-01 00:00:00', tz='UTC')

#' Per-row timestep-bin index relative to the fixed two-station epoch
#'
#' The single definition of "which bin does this timestamp fall in," shared
#' by \code{\link{mm_snap_to_bin_2s}} and \code{\link{mm_lag_2s}}.
#'
#' @param solar.time POSIXct vector of timestamps, in UTC.
#' @param timestep_days the nominal timestep, in days (see
#'   \code{\link{mm_get_timestep}}).
#' @return integer vector, the same length as \code{solar.time}
#' @keywords internal
mm_bin_index_2s <- function(solar.time, timestep_days) {
  offset_days <- as.numeric(difftime(solar.time, mm_bin_epoch_2s, units='days'))
  # round() to an integer index, rather than comparing timestamps directly,
  # means two rows match only by exact bin equality -- no tolerance to tune,
  # and no sensitivity to the sub-microsecond floating-point noise that
  # timestamp arithmetic introduces
  round(offset_days / timestep_days)
}

#' Snap timestamps onto a common per-series bin grid
#'
#' Concatenated sensor deployments can be internally regular but mutually
#' phase-shifted (e.g. one logging on the \code{:01/:16/:31/:46} minute
#' marks, the next on \code{:04/:19/:34/:49}), which otherwise looks
#' identical to a real gap. Rounding every timestamp to the nearest bin of
#' the series' nominal timestep removes that ambiguity: both phase patterns
#' land on the same grid, and a true gap becomes empty bins rather than a
#' boundary to track separately (see issue #475).
#'
#' Bins are anchored to a fixed epoch (\code{\link{mm_bin_epoch_2s}}), so two
#' series snapped independently -- e.g. upstream and downstream -- land on
#' the same absolute bins. Snapping distorts timestamps by up to half a
#' timestep; this is the same class of approximation already accepted for
#' travel-time rounding in \code{\link{mm_lag_2s}}.
#'
#' @param solar.time POSIXct vector of timestamps, in UTC, sorted ascending.
#' @return POSIXct vector, the same length as \code{solar.time}, of the
#'   snapped timestamps.
#' @keywords internal
mm_snap_to_bin_2s <- function(solar.time) {

  if(length(solar.time) < 2) {
    stop('need at least 2 rows of data to determine a timestep', call.=FALSE)
  }

  # 'modal' format is robust to the occasional large gap between deployments
  # (an outlier, not the mode). mm_get_timestep(format='modal') sometimes
  # wraps its answer in a length-1 list (a pre-existing quirk of that
  # function, not addressed here); unlist before checking length so a single
  # detected timestep is never mistaken for a detection failure
  timestep_days <- unlist(mm_get_timestep(solar.time, format='modal'))
  if(length(timestep_days) != 1 || is.na(timestep_days)) {
    stop('could not detect a single nominal timestep to snap to', call.=FALSE)
  }

  bin_index <- mm_bin_index_2s(solar.time, timestep_days)
  mm_bin_epoch_2s + as.difftime(bin_index * timestep_days, units='days')
}

#' Check that no two rows share a timestep bin
#'
#' @param bin integer vector of bin indices (see \code{\link{mm_bin_index_2s}}).
#' @param what name of the series being checked, used in the error message.
#' @keywords internal
mm_check_unique_bins_2s <- function(bin, what) {
  if(any(duplicated(bin))) {
    stop(
      what, ' has more than one row in the same timestep bin; check for ',
      'duplicate or near-duplicate timestamps', call.=FALSE)
  }
  invisible(NULL)
}

#' Detect and validate a two-station bin grid
#'
#' Establishes the timestep bin grid that the two-station functions match
#' rows on: detects the nominal timestep, confirms \code{solar.time} is
#' already on the snap-to-bin grid, and confirms no two rows land in the
#' same bin.
#'
#' @param solar.time POSIXct vector of timestamps, in UTC, already snapped via
#'   \code{\link{mm_snap_to_bin_2s}}.
#' @param what name of the series, used in error messages.
#' @param caller name of the calling function, named in the off-grid error so
#'   the user knows where to apply \code{\link{mm_snap_to_bin_2s}}.
#' @param check_sorted if TRUE, also require the bins to be ascending. Runs
#'   between the grid and duplicate checks, so unsorted input is reported as
#'   such even when it also contains duplicates.
#' @return a list with \code{timestep_days} (modal timestep, in days) and
#'   \code{bin} (integer bin index per row).
#' @keywords internal
mm_bin_grid_2s <- function(solar.time, what='solar.time', caller, check_sorted=FALSE) {

  if(length(solar.time) < 2) {
    stop('need at least 2 rows of data to determine a timestep', call.=FALSE)
  }

  # detects the timestep independently rather than calling mm_snap_to_bin_2s():
  # that function can't assume its input is already snapped -- it's what
  # produces snapped output -- so it can't rely on the round-trip check below
  timestep_days <- unlist(mm_get_timestep(solar.time, format='modal'))
  if(length(timestep_days) != 1 || is.na(timestep_days)) {
    stop('could not detect a single nominal timestep for ', what, call.=FALSE)
  }

  bin <- mm_bin_index_2s(solar.time, timestep_days)

  # Reconstruct the snapped values from the timestep detected above, rather
  # than calling mm_snap_to_bin_2s() again, to avoid redundantly detecting
  # the same timestep a second time. Already-snapped data is a fixed point
  # of mm_snap_to_bin_2s(), so this round-trip is a no-op iff solar.time is
  # already on the grid. 1-second tolerance is far above the ~microsecond
  # floating-point noise from timestamp arithmetic and far below any real
  # timestep.
  snapped <- mm_bin_epoch_2s + as.difftime(bin * timestep_days, units='days')
  if(max(abs(as.numeric(snapped) - as.numeric(solar.time))) > 1) {
    stop(
      'solar.time is not on a snap-to-bin grid; call mm_snap_to_bin_2s() ',
      'on it before ', caller, '() (see issue #475)', call.=FALSE)
  }

  if(check_sorted && is.unsorted(bin)) {
    stop('data must be sorted ascending by solar.time', call.=FALSE)
  }
  mm_check_unique_bins_2s(bin, what)

  list(timestep_days = timestep_days, bin = bin)
}

#' Compute the per-row upstream lag for two-station data
#'
#' The upstream observation matching downstream row \code{i} was recorded
#' \code{lag[i] <- round(travel.time[i] / timestep_days)} timesteps earlier,
#' where \code{timestep_days} is the modal timestep of \code{solar.time}.
#' This is the single definition of that calculation -- \code{mm_align_2s}
#' and \code{\link{mm_validate_data_2station}} both route through it, so the
#' two cannot drift apart.
#'
#' \code{shift_idx[i]} is \code{NA} (folded into \code{has_leadin = FALSE})
#' when no row exists at row \code{i}'s target timestep bin -- whether
#' because the shift reaches before the start of the data or into a real
#' gap partway through the series. Both cases mean the same thing to every
#' caller: this row cannot be modeled for lack of a real upstream match.
#'
#' @param solar.time POSIXct vector of timestamps, in UTC, sorted ascending,
#'   already snapped to a single nominal timestep via
#'   \code{\link{mm_snap_to_bin_2s}}.
#' @param travel.time numeric vector of reach travel times, in days, the same
#'   length as \code{solar.time}.
#' @return a list with \code{timestep_days} (modal timestep, in days),
#'   \code{lag} (per-row lag, in timesteps), \code{shift_idx} (per-row index
#'   into the input rows supplying that row's upstream values, or \code{NA}
#'   where no row exists at the target timestep bin), and \code{has_leadin}
#'   (logical, \code{!is.na(shift_idx)})
#' @keywords internal
mm_lag_2s <- function(solar.time, travel.time) {

  # requires solar.time already snapped to a grid: otherwise a computed
  # target_bin might not match any real timestamp's bin even when a
  # matching row conceptually exists, silently under-counting lead-in
  grid <- mm_bin_grid_2s(solar.time, caller='mm_lag_2s')
  timestep_days <- grid$timestep_days
  bin_index <- grid$bin

  lag <- round(travel.time / timestep_days)
  target_bin <- bin_index - lag
  # Match by timestep bin, not by subtracting lag[i] from i's row position:
  # a row offset only equals a time offset when every row is exactly one
  # nominal timestep apart, an assumption gaps or phase-shifted deployments
  # break, silently pairing the wrong upstream row with a downstream one.
  # Bin matching instead looks up whichever row actually sits lag[i]
  # timesteps back in time.
  shift_idx <- match(target_bin, bin_index)

  list(
    timestep_days = timestep_days,
    lag = lag,
    shift_idx = shift_idx,
    has_leadin = !is.na(shift_idx))
}

#' Align two-station data onto its 06:00 day window
#'
#' Applies, in order: the per-row lead-in test, two-station day labeling, the
#' travel-time ceiling, and the whole-day completeness requirement. The
#' result describes which input rows are modeled, which input rows supply
#' their upstream values, and how those modeled rows partition into days.
#'
#' Rows are dropped for lack of lead-in individually, not by trimming a
#' fixed prefix sized to the dataset's worst-case lag: a row is modeled
#' whenever \code{\link{mm_lag_2s}} finds it a real upstream match, whether
#' the shift would otherwise reach before the start of the data or into a
#' gap partway through.
#'
#' Whole days are dropped, with a message, in two cases:
#' \itemize{
#'   \item the day's longest travel time exceeds
#'     \code{max_travel_time_hours} (see \code{\link{mm_max_travel_time}});
#'   \item the day does not fill its 24-hour window, i.e. it holds fewer than
#'     \code{round(1 / timestep_days)} modeled rows. Partial days arise at
#'     the edges of any dataset whose bounds don't fall on 06:00, and from
#'     missing sensor readings; they're dropped because the Stan model
#'     expects fixed-shape \code{n_obs x n_days} matrices. Interpolating
#'     small within-day gaps instead of dropping is out of scope here.
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
#'   \code{keep}), \code{n_obs} (modeled rows per day), \code{n_days},
#'   \code{timestep_days}, and \code{removed} (a data.frame of \code{date}
#'   and \code{errors} naming the days dropped here and why)
#' @keywords internal
mm_align_2s <- function(data, max_travel_time_hours=mm_max_travel_time_default) {

  mm_check_max_travel_time_hours(max_travel_time_hours)

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

  # Announcing a dropped day and recording it are one action, not two: a day
  # messaged about but not recorded is exactly the day that disappears from
  # the results with nothing to explain it. headline completes "dropping N
  # day(s) ___"; reason is the stored per-day explanation; detail is the
  # per-day specifics both use
  removed <- mm_no_removed_days_2s()
  drop_days <- function(dates, detail, headline, reason) {
    message(paste0(
      'dropping ', length(dates), ' day(s) ', headline, ': ',
      paste(sprintf('%s (%s)', as.character(dates), detail), collapse=', ')))
    removed <<- rbind(removed, data.frame(
      date=as.Date(dates), errors=paste0(reason, ' (', detail, ')'),
      stringsAsFactors=FALSE))
  }

  # a day none of whose rows have lead-in never enters keep at all. Days before
  # the first modelable row are exempt: that prefix is the lead-in block,
  # supplied on purpose to be drawn from rather than modeled
  all_dates <- mm_date_2s(solar_time)
  no_leadin <- sort(unique(all_dates[all_dates >= all_dates[keep[1]] & !(all_dates %in% date)]))
  if(length(no_leadin) > 0) {
    drop_days(
      no_leadin, sprintf('%d row(s) supplied', tabulate(match(all_dates, no_leadin), length(no_leadin))),
      'with no upstream lead-in at all',
      "no upstream observation at any row's travel-time offset")
  }

  # travel-time ceiling: drop the offending day, not the whole dataset
  max_travel_time_days <- max_travel_time_hours / 24
  worst_by_day <- tapply(travel_time[keep], date, max)
  over_ceiling <- worst_by_day[worst_by_day > max_travel_time_days]
  if(length(over_ceiling) > 0) {
    drop_days(
      names(over_ceiling), sprintf('%.2f hours', unname(over_ceiling) * 24),
      paste0('whose travel.time exceeds the ', max_travel_time_hours, '-hour ceiling'),
      paste0('travel.time exceeds the ', max_travel_time_hours, '-hour ceiling'))
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
    drop_days(
      names(partial), sprintf('%d of %d expected observations', unname(partial), expected_n_obs),
      'that do not fill the 06:00-06:00 window',
      'does not fill the 06:00-06:00 window')
    complete <- !(as.character(date) %in% names(partial))
    keep <- keep[complete]
    date <- date[complete]
  }

  if(length(keep) == 0) {
    stop(paste0(
      'no complete days remain after applying the ', max_travel_time_hours,
      '-hour travel.time ceiling and the 06:00-06:00 day-window requirement'), call.=FALSE)
  }

  removed <- removed[order(removed$date), , drop=FALSE]
  rownames(removed) <- NULL

  list(
    keep = keep,
    shift_idx = lagged$shift_idx[keep],
    date = date,
    n_obs = expected_n_obs,
    n_days = length(unique(date)),
    timestep_days = lagged$timestep_days,
    removed = removed)
}
