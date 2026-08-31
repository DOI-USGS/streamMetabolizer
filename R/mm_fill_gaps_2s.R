#' @include mm_lag_2s.R
NULL

#' Default and maximum gap-filling tolerance for two-station models
#'
#' Bounds how much consecutive missing data may be interpolated across
#' before a gap is left in place; longer gaps are left as-is, and the days
#' holding them are dropped by \code{\link{mm_align_2s}}'s whole-day
#' completeness requirement. There is no middle tier: nothing partially
#' filled ever reaches Stan.
#'
#' \code{mm_max_gap_hours_default} is the canonical default (currently 1
#' hour) -- looser than Bishop et al. (2026)'s single-point/15-minute rule.
#' \code{mm_max_gap_hours_cap} bounds what a user may ask for: linear
#' interpolation across more than a couple of hours of a diel curve invents
#' structure rather than bridging a blip.
#'
#' @keywords internal
#' @name mm_max_gap
NULL

#' @rdname mm_max_gap
mm_max_gap_hours_default <- 1
# on the source record, missing-data runs are sharply bimodal (sub-2-hour
# blips, or multi-day sensor outages, with little in between); one hour
# recovers most of what's recoverable while keeping every interpolated
# stretch short enough that linear interpolation across a diel DO or light
# curve stays defensible

#' @rdname mm_max_gap
mm_max_gap_hours_cap <- 2

#' Validate a gap-filling tolerance
#'
#' Shared by \code{specs()} and \code{\link{mm_fill_gaps_2s}} so the two
#' can't disagree about what tolerance is legal.
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
#' Fills gaps short enough to bridge by linear interpolation, so a day
#' marred by a brief sensor dropout can still be modeled. Gaps longer than
#' \code{max_gap_hours} are left as-is and the day is dropped later by
#' \code{\link{mm_align_2s}}'s completeness check -- there's no partial-fill
#' middle ground, since the Stan model's \code{n_obs x n_days} matrices
#' can't mask individual missing values.
#'
#' A gap can be a missing row (no observation at a timestep bin) or a
#' missing value (a row exists but a column is \code{NA}); both are treated
#' as the same event. Gap length is measured in timestep bins between the
#' bracketing observations, not row count, so a run that looks short by row
#' count but spans a multi-day outage is correctly refused.
#'
#' Gaps at the start or end of the record are never filled, since
#' interpolation needs observations on both sides; those edge days are
#' already dropped by the completeness requirement.
#'
#' @section Invariant: This function never returns a row with \code{NA} in
#'   an interpolated column -- a partially-fillable inserted row is dropped
#'   rather than kept.
#'
#' @section Interpolating light: The \code{light} column is interpolated
#'   like any other, but accuracy depends on timing:
#'   \code{\link{mm_format_data_2s}} fills raw light before converting to
#'   the within-day proportion (accurate), while hand-formatted data
#'   arrives pre-converted, so interpolating it approximates a proportion
#'   whose denominator is already short the missing terms. Supplying raw
#'   light is the more accurate route.
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

  # gaps are measured in timestep bins, so solar.time must already be on the
  # bin grid -- the same precondition, checked the same way, as mm_lag_2s()
  grid <- mm_bin_grid_2s(data_v$solar.time, caller='mm_fill_gaps_2s', check_sorted=TRUE)
  timestep_days <- grid$timestep_days
  bin <- grid$bin

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
  #
  # gap length is measured in bins (step = diff(bin)), not row count, so a
  # run that looks short in rows is still correctly rejected if it spans a
  # multi-day outage on the bin grid
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
    # solar.time for an inserted row comes from the bin grid directly, not
    # interpolated -- the bin index alone determines the timestamp
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
  # a short missing-row gap can overlap a longer missing-value gap in some
  # column (a brief downstream dropout inside a multi-hour upstream outage,
  # say), leaving an inserted row that can't be fully filled. Without
  # dropping it back out, a partly-filled row would satisfy the completeness
  # row count while still carrying NA into Stan
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
