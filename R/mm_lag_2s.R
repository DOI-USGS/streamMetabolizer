#' Default hour at which a two-station day begins
#'
#' Two-station days run 24 hours from a start hour chosen during data prep;
#' the default, 6 (06:00-06:00), follows Bishop et al. (2026).
#'
#' @keywords internal
mm_day_start_2s <- 6

#' Validate a two-station day start hour
#'
#' @param day_start_hour the value to check: a single number in [0, 24).
#' @keywords internal
mm_check_day_start_hour_2s <- function(day_start_hour) {
  if(!is.numeric(day_start_hour) || length(day_start_hour) != 1 || is.na(day_start_hour)) {
    stop('day_start_hour must be a single non-NA number', call.=FALSE)
  }
  if(day_start_hour < 0 || day_start_hour >= 24) {
    stop('day_start_hour must be >= 0 and < 24', call.=FALSE)
  }
  invisible(NULL)
}

#' Assign two-station day labels to timestamps
#'
#' Labels each timestamp with the two-station day it belongs to: the 24-hour
#' day beginning at \code{day_start_hour}. With the default start of 6,
#' 05:45 belongs to the *previous* calendar date and 06:00 to the current
#' one.
#'
#' @param solar.time POSIXct vector of timestamps, in UTC (as enforced by
#'   \code{\link{mm_validate_data}}).
#' @param day_start_hour hour of the day, in [0, 24), at which each
#'   two-station day begins. Required, so every caller states the boundary it
#'   uses.
#' @return a Date vector of two-station day labels, the same length as
#'   \code{solar.time}
#' @keywords internal
mm_date_2s <- function(solar.time, day_start_hour) {
  as.Date(solar.time - as.difftime(day_start_hour, units='hours'), tz='UTC')
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
#' @param check_sorted if TRUE, also require the bins to be ascending. Runs
#'   between the grid and duplicate checks, so unsorted input is reported as
#'   such even when it also contains duplicates.
#' @return a list with \code{timestep_days} (modal timestep, in days) and
#'   \code{bin} (integer bin index per row).
#' @keywords internal
mm_bin_grid_2s <- function(solar.time, what='solar.time', check_sorted=FALSE) {

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
    # reachable from exported functions, so point to the exported preparation
    # step rather than the internal snapping helper
    stop(
      'solar.time is not on a regular timestep grid; prepare the data with ',
      'mm_format_data_2s(), which snaps timestamps to a common grid',
      call.=FALSE)
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
  grid <- mm_bin_grid_2s(solar.time)
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

#' Find the days whose light windows reach light in an earlier day
#'
#' Each row's light window is the set of rows its travel-time light sum
#' covers: the rows from its upstream match (see \code{\link{mm_lag_2s}})
#' through the row itself. A day fails when any of its rows has a window
#' reaching a row in an earlier day whose light is above zero (\code{"lit"})
#' or \code{NA} (\code{"unknown"}). A day with both is reported as
#' \code{"lit"} only. Rows without lead-in have no window and are not tested.
#'
#' @param solar.time POSIXct vector of timestamps, in UTC, sorted ascending,
#'   already snapped to a single nominal timestep via
#'   \code{\link{mm_snap_to_bin_2s}}.
#' @param light numeric vector, the same length as \code{solar.time}, of raw
#'   light (not the within-day proportion).
#' @param travel.time numeric vector of reach travel times, in days, the
#'   same length as \code{solar.time}.
#' @param day_start_hour hour of the day, in [0, 24), at which each
#'   two-station day begins.
#' @return a data.frame with one row per failing day, sorted by date:
#'   \code{date} (Date), \code{reason} (\code{"lit"} or \code{"unknown"}),
#'   and \code{n_rows} (the number of that day's rows whose windows fail)
#' @keywords internal
mm_light_check_2s <- function(solar.time, light, travel.time, day_start_hour) {

  mm_check_day_start_hour_2s(day_start_hour)
  if(is.unsorted(solar.time)) {
    stop('data must be sorted ascending by solar.time', call.=FALSE)
  }
  lagged <- mm_lag_2s(solar.time, travel.time)
  # over all rows, not only those with lead-in: window rows are often lead-in
  # rows that are never modeled themselves
  day <- mm_date_2s(solar.time, day_start_hour)

  # sorted rows put each day in one contiguous run, so the earlier-day part
  # of row i's window is rows shift_idx[i] through first_of_day[i] - 1, and
  # cumulative counts give each window's lit and unknown rows in O(n)
  first_of_day <- match(day, day)
  from <- lagged$shift_idx
  to <- first_of_day - 1L
  reaches <- lagged$has_leadin & from <= to
  count_between <- function(flag) {
    cum <- c(0, cumsum(flag))
    out <- numeric(length(flag))
    out[reaches] <- cum[to[reaches] + 1] - cum[from[reaches]]
    out
  }
  unknown <- is.na(light)
  row_lit <- count_between(!unknown & light > 0) > 0
  row_unknown <- count_between(unknown) > 0

  failing <- row_lit | row_unknown
  dates <- sort(unique(day[failing]))
  day_reason <- c('unknown', 'lit')[(dates %in% day[row_lit]) + 1]
  idx <- match(day, dates)
  counted <- failing & ifelse(day_reason[idx] %in% 'lit', row_lit, row_unknown)
  data.frame(
    date=dates,
    reason=day_reason,
    n_rows=tabulate(idx[counted], length(dates)),
    stringsAsFactors=FALSE)
}

#' Align two-station data onto its 24-hour day window
#'
#' Applies, in order: the per-row lead-in test, two-station day labeling, the
#' light-window check, and the whole-day completeness requirement. The
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
#'   \item the day is listed in \code{light_check_failed}: one of its rows'
#'     light windows reaches a row in an earlier day with light above zero or
#'     unknown (see \code{\link{mm_light_check_2s}});
#'   \item the day does not fill its 24-hour window, i.e. it holds fewer than
#'     \code{round(1 / timestep_days)} modeled rows. Partial days arise at
#'     the edges of any dataset whose bounds don't fall on the day start
#'     hour, and from missing rows (a reading logged as \code{NA} in a row
#'     that is present does not make a day partial); they're dropped because
#'     the Stan model expects fixed-shape \code{n_obs x n_days} matrices. Short
#'     gaps are filled beforehand by \code{\link{mm_fill_gaps_2s}}, if at all, not here.
#' }
#'
#' @param data data.frame (units already stripped) containing at least
#'   \code{solar.time} and \code{travel.time}, sorted ascending by
#'   \code{solar.time}.
#' @param day_start_hour hour of the day, in [0, 24), at which each
#'   two-station day begins.
#' @param light_check_failed the days that failed the light-window check, as
#'   returned by \code{\link{mm_light_check_2s}}, or \code{NULL} if the check
#'   was not run.
#' @return a list with \code{keep} (integer indices of the modeled rows, into
#'   the rows of \code{data}), \code{shift_idx} (integer indices, also into
#'   \code{data}, of the rows supplying each modeled row's upstream values),
#'   \code{date} (Date vector of two-station day labels, parallel to
#'   \code{keep}), \code{n_obs} (modeled rows per day), \code{n_days},
#'   \code{timestep_days}, \code{removed} (a data.frame of \code{date}
#'   and \code{errors} naming the days dropped here and why), and
#'   \code{light_check} (TRUE if \code{light_check_failed} was supplied)
#' @keywords internal
mm_align_2s <- function(data, day_start_hour, light_check_failed=NULL) {

  mm_check_day_start_hour_2s(day_start_hour)

  solar_time <- data$solar.time
  travel_time <- data$travel.time

  lagged <- mm_lag_2s(solar_time, travel_time)

  # per-row lead-in: keep every row whose upstream shift lands within the data
  keep <- which(lagged$has_leadin)
  if(length(keep) == 0) {
    # mm_validate_data_2station() raises a more informative version of this
    # error before aligning; this guards direct calls
    stop('no rows have enough upstream lead-in data to be modeled', call.=FALSE)
  }
  date <- mm_date_2s(solar_time[keep], day_start_hour)

  # Announcing a dropped day and recording it are one action, not two: a day
  # messaged about but not recorded is exactly the day that disappears from
  # the results with nothing to explain it. headline completes "dropping N
  # day(s) ___"; reason is the stored per-day explanation; detail is the
  # per-day specifics both use. max_listed caps the days named in the message;
  # every day is still recorded
  removed <- mm_no_removed_days_2s
  drop_days <- function(dates, detail, headline, reason, max_listed=Inf) {
    listed <- utils::head(sprintf('%s (%s)', as.character(dates), detail), max_listed)
    message(paste0(
      'dropping ', length(dates), ' day(s) ', headline, ': ',
      paste(listed, collapse=', '),
      if(length(dates) > length(listed)) paste0(', and ', length(dates) - length(listed), ' more') else ''))
    removed <<- rbind(removed, data.frame(
      date=as.Date(dates), errors=paste0(reason, ' (', detail, ')'),
      stringsAsFactors=FALSE))
  }

  # a day none of whose rows have lead-in never enters keep at all. Days before
  # the first modelable row are exempt: that prefix is the lead-in block,
  # supplied on purpose to be drawn from rather than modeled
  all_dates <- mm_date_2s(solar_time, day_start_hour)
  no_leadin <- sort(unique(all_dates[all_dates >= all_dates[keep[1]] & !(all_dates %in% date)]))
  if(length(no_leadin) > 0) {
    drop_days(
      no_leadin, sprintf('%d row(s) supplied', tabulate(match(all_dates, no_leadin), length(no_leadin))),
      'with no upstream lead-in at all',
      "no upstream observation at any row's travel-time offset")
  }

  # light window: drop each failing day that is still being modeled; days
  # already dropped above, or never modeled, need no second entry
  n_light_dropped <- 0
  if(!is.null(light_check_failed)) {
    failed <- light_check_failed[light_check_failed$date %in% date, , drop=FALSE]
    lit <- failed[failed$reason == 'lit', , drop=FALSE]
    if(nrow(lit) > 0) {
      drop_days(
        lit$date, sprintf('%d row(s)', lit$n_rows),
        'whose travel-time light windows reach light above zero in an earlier day',
        'light window reaches light above zero in an earlier day', max_listed=10)
    }
    unknown <- failed[failed$reason == 'unknown', , drop=FALSE]
    if(nrow(unknown) > 0) {
      drop_days(
        unknown$date, sprintf('%d row(s)', unknown$n_rows),
        'whose travel-time light windows reach light unknown in an earlier day',
        'light unknown in an earlier day of the light window', max_listed=10)
    }
    n_light_dropped <- nrow(failed)
    in_bounds <- !(date %in% failed$date)
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
      'that do not fill the 24-hour day window',
      'does not fill the 24-hour day window')
    complete <- !(as.character(date) %in% names(partial))
    keep <- keep[complete]
    date <- date[complete]
  }

  if(length(keep) == 0) {
    if(n_light_dropped > 0) {
      stop(paste0(
        "Every day was dropped: each day's travel-time light window reaches light from an ",
        "earlier day, which would contaminate the day's light proportions. Check that ",
        "`travel.time` is in days and not too long, and move `day_start_hour` to a time ",
        "with no light in the hours before it."), call.=FALSE)
    }
    stop('no complete days remain: every day was dropped (see the messages above)', call.=FALSE)
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
    removed = removed,
    light_check = !is.null(light_check_failed))
}
