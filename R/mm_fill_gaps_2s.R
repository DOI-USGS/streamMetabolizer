#' @include mm_lag_2s.R
NULL

#' Default and maximum gap-filling tolerance for two-station models
#'
#' The tolerance bounds how much consecutive missing data may be interpolated
#' across before a gap is left in place. Gaps within the tolerance are filled;
#' longer ones are not, and the days holding them are dropped by
#' \code{\link{mm_align_2s}}'s whole-day completeness requirement. There is no
#' middle tier: nothing partially filled ever reaches Stan.
#'
#' \code{mm_max_gap_hours_default} is the canonical default (currently 1).
#' Bishop et al. (2026) interpolate a single missing point and drop the day at
#' two or more consecutive missing points -- 15 minutes at the source data's
#' native timestep. One hour is a deliberate loosening of that rule rather
#' than a restatement of it: on the source record, missing-data runs are
#' sharply bimodal (sub-2-hour blips, or multi-day sensor outages, with
#' little in between), so a one-hour tolerance recovers most of what is
#' recoverable while keeping every interpolated stretch short enough that
#' linear interpolation across a diel DO or light curve stays defensible.
#'
#' \code{mm_max_gap_hours_cap} bounds what a user may ask for. Linear
#' interpolation across more than a couple of hours of a diel curve --
#' spanning dawn, dusk, or the midday peak -- invents structure rather than
#' bridging a blip, and the accompanying day would be better dropped.
#'
#' \code{\link{specs}} cannot reference \code{mm_max_gap_hours_default}
#' directly, for the same reason \code{\link{mm_max_travel_time}} documents:
#' \code{prefer_missing} flags any formal whose unevaluated default is a bare
#' symbol as "no real default", producing a warning that tells users to set
#' the option in \code{\link{revise}()} instead -- backwards, since setting it
#' at \code{specs()}-creation time is the intended use. \code{specs()}'s
#' literal default is kept in sync with this constant by hand; a guard test in
#' \code{test-mm_fill_gaps_2s.R} (search for "max_gap_hours default") fails if
#' the two diverge.
#'
#' @keywords internal
#' @name mm_max_gap
NULL

#' @rdname mm_max_gap
mm_max_gap_hours_default <- 1

#' @rdname mm_max_gap
mm_max_gap_hours_cap <- 2

#' Validate a gap-filling tolerance
#'
#' Shared by \code{\link{mm_fill_gaps_2s}} and \code{\link{specs}} so that the
#' two entry points cannot disagree about what a legal tolerance is:
#' \code{specs()} checks it at specs-creation time (early, as it does for
#' \code{max_travel_time_hours}), and \code{mm_fill_gaps_2s} checks it again
#' because \code{\link{mm_format_data_2s}} takes the tolerance directly rather
#' than through a \code{specs} list.
#'
#' @param max_gap_hours the value to check.
#' @keywords internal
mm_check_max_gap_hours <- function(max_gap_hours) {
  if(!is.numeric(max_gap_hours) || length(max_gap_hours) != 1 || is.na(max_gap_hours)) {
    stop('max_gap_hours must be a single non-NA number', call.=FALSE)
  }
  if(max_gap_hours <= 0) {
    stop('max_gap_hours must be > 0', call.=FALSE)
  }
  if(max_gap_hours > mm_max_gap_hours_cap) {
    stop(paste0(
      'max_gap_hours must be <= ', mm_max_gap_hours_cap, ' hours; linear interpolation ',
      'across a longer stretch of a diel DO or light curve invents structure rather than ',
      'bridging a blip, and the affected day is better dropped'), call.=FALSE)
  }
  invisible(NULL)
}

#' Interpolate short gaps in two-station data
#'
#' Fills gaps in a two-station data.frame that are short enough to bridge by
#' linear interpolation, so that a day marred by a brief sensor dropout can be
#' modeled instead of discarded. Gaps longer than \code{max_gap_hours} are
#' left exactly as they were, and the days holding them are dropped later by
#' \code{\link{mm_align_2s}}'s whole-day completeness requirement -- the
#' policy is two-tier, with no ragged middle ground, because the Stan model's
#' \code{n_obs x n_days} matrices have no masking mechanism.
#'
#' Missing data takes two forms, and this function treats them identically
#' because they describe the same physical event:
#' \itemize{
#'   \item \emph{missing rows}: no row exists at a timestep bin at all, which
#'     is what \code{\link{mm_format_data_2s}} produces when the downstream
#'     sonde logged nothing there. Rows are inserted, and their
#'     \code{solar.time} is computed exactly from the bin grid rather than
#'     interpolated.
#'   \item \emph{missing values}: a row exists but a column is \code{NA} --
#'     from an unmatched upstream observation, or a sensor that logged
#'     \code{NA} rather than skipping the timestep.
#' }
#'
#' Both are measured the same way, in timestep bins between the two
#' observations bracketing the gap, so a gap is judged by how much time it
#' spans rather than by how many rows happen to be present across it. A
#' consequence worth noting: an \code{NA} run that looks short in row terms
#' but straddles a multi-day outage is correctly refused, because its
#' bracketing observations are far apart on the bin grid even though they are
#' adjacent rows.
#'
#' Gaps at the very start or end of the record are never filled: linear
#' interpolation needs an observation on both sides, and extrapolating from
#' one would fabricate data rather than bridge it. Those edge days are already
#' dropped by the completeness requirement.
#'
#' @section Invariant: this function never introduces a row carrying \code{NA}
#'   in any interpolated column. Where a short missing-row gap overlaps a
#'   longer missing-value gap in some column -- a brief downstream dropout
#'   inside a multi-hour upstream outage, say -- the inserted rows cannot be
#'   fully filled, so they are removed again. Rows supplied by the caller are
#'   never dropped, only added to. Without this, a partly-filled row would
#'   satisfy the completeness row count while still carrying \code{NA} into
#'   Stan.
#'
#' @section Interpolating light: the \code{light} column is interpolated like
#'   any other, but what that means depends on what it holds.
#'   \code{\link{mm_format_data_2s}} fills gaps in \emph{raw} light before
#'   converting it to the within-day proportion, which is the meaningful
#'   order: the proportion is then computed from a complete series over a
#'   complete day total. Data formatted by hand arrives with \code{light}
#'   already converted, and interpolating it is an approximation -- an
#'   interpolated proportion is not the proportion that would have been
#'   computed from the missing raw light, because the day total in the
#'   denominator is itself short by the missing terms. It is a good deal
#'   better than the alternative, which is to leave the day unmodelable, but
#'   supplying raw light and letting it be converted is the more accurate
#'   route.
#'
#' @param data data.frame, unitted or not, with \code{solar.time} already
#'   snapped to a single nominal timestep via \code{\link{mm_snap_to_bin_2s}}
#'   and sorted ascending. Every other column is treated as a numeric time
#'   series to be interpolated.
#' @param max_gap_hours the gap-filling tolerance, in hours: runs of missing
#'   data spanning no more than this are interpolated, longer runs are left in
#'   place. See \code{\link{mm_max_gap}}.
#' @return \code{data} with the same columns, column order, and units, and
#'   with rows inserted where a short missing-row gap could be filled.
#' @importFrom stats approx
#' @importFrom unitted u v is.unitted get_units
#' @keywords internal
mm_fill_gaps_2s <- function(data, max_gap_hours=mm_max_gap_hours_default) {

  mm_check_max_gap_hours(max_gap_hours)

  # units may sit on the columns alone (the shape mm_format_data_2s() and
  # two_station_example use) or on the frame as well (a true unitted
  # data.frame); v() strips both, so record which it was and restore in kind.
  # Testing the frame alone would silently drop column units off the common
  # case, since is.unitted() is FALSE for a plain frame of unitted columns
  frame_unitted <- is.unitted(data)
  col_unitted <- vapply(names(data), function(col) is.unitted(data[[col]]), logical(1))
  col_units <- vapply(
    names(data),
    function(col) if(is.unitted(data[[col]])) get_units(data[[col]]) else NA_character_,
    character(1))
  data_v <- v(data)

  if(nrow(data_v) < 2) {
    stop('need at least 2 rows of data to determine a timestep', call.=FALSE)
  }

  timestep_days <- unlist(mm_get_timestep(data_v$solar.time, format='modal'))
  if(length(timestep_days) != 1 || is.na(timestep_days)) {
    stop('could not detect a single nominal timestep for solar.time', call.=FALSE)
  }

  # gaps are measured in timestep bins, so solar.time must already be on the
  # bin grid -- the same precondition, checked the same way, as mm_lag_2s()
  snapped <- mm_snap_to_bin_2s(data_v$solar.time)
  if(max(abs(as.numeric(snapped) - as.numeric(data_v$solar.time))) > 1) {
    stop(
      'solar.time is not on a snap-to-bin grid; call mm_snap_to_bin_2s() ',
      'on it before mm_fill_gaps_2s() (see issue #475)', call.=FALSE)
  }

  bin <- mm_bin_index_2s(data_v$solar.time, timestep_days)
  if(is.unsorted(bin)) {
    stop('data must be sorted ascending by solar.time', call.=FALSE)
  }
  if(any(duplicated(bin))) {
    stop(
      'solar.time has more than one row in the same timestep bin; check for ',
      'duplicate or near-duplicate timestamps', call.=FALSE)
  }

  # the tolerance, as a whole number of timestep bins. the epsilon keeps an
  # exactly-representable ratio (a 1-hour tolerance on a 15-minute timestep is
  # exactly 4) from floor()ing down to 3 on floating-point noise
  max_gap_bins <- floor(max_gap_hours / (timestep_days * 24) + 1e-8)
  if(max_gap_bins < 1) {
    # tolerance finer than the timestep: nothing is fillable, and saying so is
    # more useful than silently returning the input unchanged
    message(paste0(
      'max_gap_hours (', max_gap_hours, ') is shorter than one timestep (',
      signif(timestep_days * 24, 3), ' hours); no gaps can be filled'))
    return(data)
  }

  target_cols <- setdiff(names(data_v), 'solar.time')

  # --- insert candidate rows at short runs of missing bins ----------------
  #
  # only short runs get rows, rather than laying down the full grid and
  # pruning afterwards: a multi-day outage would otherwise become thousands
  # of all-NA rows, and a day made entirely of them would satisfy the
  # completeness row count while carrying nothing but NA into Stan
  step <- diff(bin)
  gap_at <- which(step > 1)
  n_missing <- step[gap_at] - 1
  short_run <- n_missing <= max_gap_bins
  new_bins <- if(any(short_run)) {
    unlist(lapply(which(short_run), function(k) bin[gap_at[k]] + seq_len(n_missing[k])))
  } else {
    numeric(0)
  }
  n_gaps_skipped <- sum(!short_run)

  work <- data_v
  work$.inserted <- FALSE
  work_bin <- bin
  if(length(new_bins) > 0) {
    added <- data_v[rep(NA_integer_, length(new_bins)), , drop=FALSE]
    added$solar.time <- mm_bin_epoch_2s + as.difftime(new_bins * timestep_days, units='days')
    added$.inserted <- TRUE
    work <- rbind(work, added)
    work_bin <- c(bin, new_bins)
    ord <- order(work_bin)
    work <- work[ord, , drop=FALSE]
    work_bin <- work_bin[ord]
    rownames(work) <- NULL
  }

  # --- interpolate every column over gaps within the tolerance ------------
  #
  # for each row, prev_anchor/next_anchor are the bins of the nearest observed
  # values on either side (the constant-interpolation idiom calc_light_merged()
  # uses to measure distance to the nearest observation). Their separation is
  # what the tolerance is tested against, which is what makes a missing row and
  # a missing value the same thing here. rule=1 leaves both anchors NA beyond
  # the ends of the observed range, so leading and trailing gaps are never
  # filled -- deliberately not calc_light_merged()'s rule=2, which would
  # extrapolate a constant off the end of the record
  n_values_filled <- 0L
  for(col in target_cols) {
    x <- work[[col]]
    if(!is.numeric(x)) next
    obs <- !is.na(x)
    if(sum(obs) < 2) next

    interp <- approx(
      x=work_bin[obs], y=x[obs], xout=work_bin, method='linear', rule=1)$y
    prev_anchor <- approx(
      x=work_bin[obs], y=work_bin[obs], xout=work_bin, method='constant', f=0, rule=1)$y
    next_anchor <- approx(
      x=work_bin[obs], y=work_bin[obs], xout=work_bin, method='constant', f=1, rule=1)$y

    fillable <- !obs & !is.na(prev_anchor) & !is.na(next_anchor) &
      (next_anchor - prev_anchor - 1) <= max_gap_bins
    if(any(fillable)) {
      x[fillable] <- interp[fillable]
      work[[col]] <- x
      n_values_filled <- n_values_filled + sum(fillable)
    }
  }

  # --- enforce the no-NA-row invariant ------------------------------------
  incomplete <- rowSums(is.na(work[target_cols])) > 0
  drop_back <- work$.inserted & incomplete
  if(any(drop_back)) {
    n_values_filled <- n_values_filled - sum(!is.na(work[drop_back, target_cols, drop=FALSE]))
    work <- work[!drop_back, , drop=FALSE]
  }

  n_rows_inserted <- sum(work$.inserted)
  work$.inserted <- NULL
  work <- work[names(data_v)]
  rownames(work) <- NULL

  # count only values filled in rows the caller supplied, so the two numbers
  # in the message don't double-count the same interpolation
  n_values_filled <- n_values_filled - n_rows_inserted * length(target_cols)

  if(n_rows_inserted > 0 || n_values_filled > 0) {
    message(paste0(
      'filling gaps of up to ', max_gap_hours, ' hour(s): inserted ', n_rows_inserted,
      ' row(s) at missing timesteps and interpolated ', n_values_filled,
      ' missing value(s)',
      if(n_gaps_skipped > 0) paste0(
        '; ', n_gaps_skipped, ' gap(s) exceed the tolerance and are left in place')
      else ''))
  }

  if(frame_unitted) {
    work <- u(work, col_units[names(work)])
  } else if(any(col_unitted)) {
    for(col in names(work)[col_unitted[names(work)]]) {
      work[[col]] <- u(work[[col]], col_units[[col]])
    }
  }

  work
}
