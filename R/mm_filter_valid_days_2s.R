#' @include mm_lag_2s.R mm_modeled_rows_2s.R
NULL

#' Day-validity tests applicable to two-station models
#'
#' \code{mm_day_tests_2s_default} is the default; \code{mm_day_tests_2s_allowed}
#' bounds what a user may ask for. Both are a chosen subset of
#' \code{\link{mm_is_valid_day}}'s five tests, not the shared default inherited
#' wholesale.
#'
#' \code{full_day} and \code{even_timesteps} are excluded because they would
#' reject nearly every two-station day, not merely fail to help: the first
#' tests a 4 AM / 28-hour diel window, which is not the 24-hour 06:00-06:00
#' window a two-station day occupies, and the second requires evenly spaced
#' timesteps, which two-station tolerates gaps in on purpose.
#' \code{pos_discharge} is allowed but off by default, a no-op until
#' two-station data carries a \code{discharge} column (issue #475).
#'
#' @keywords internal
#' @name mm_day_tests_2s
NULL

#' @rdname mm_day_tests_2s
mm_day_tests_2s_default <- c('complete_data', 'pos_depth')

#' @rdname mm_day_tests_2s
mm_day_tests_2s_allowed <- c('complete_data', 'pos_discharge', 'pos_depth')

#' Validate a two-station day_tests selection
#'
#' Called both when specs are created and again when they are used, since the
#' consumer can be called directly with tests that never passed through a
#' \code{specs} list.
#'
#' @param day_tests the value to check. \code{NULL} or \code{character(0)}
#'   means no tests and is always legal.
#' @keywords internal
mm_check_day_tests_2s <- function(day_tests) {
  # NULL as well as character(0): c() is the package's idiom for "no tests"
  # and evaluates to NULL, so rejecting it would make two-station the odd one out
  if(length(day_tests) == 0) return(invisible(NULL))
  if(!is.character(day_tests)) {
    stop('day_tests must be a character vector', call.=FALSE)
  }
  bad <- setdiff(day_tests, mm_day_tests_2s_allowed)
  if(length(bad) > 0) {
    stop(paste0(
      'day_tests for two-station models may only include ',
      paste(mm_day_tests_2s_allowed, collapse=', '), '; got ',
      paste(bad, collapse=', '), '. full_day and even_timesteps describe ',
      'one-station conventions that two-station does not share (see ?mm_day_tests_2s)'),
      call.=FALSE)
  }
  invisible(NULL)
}

#' Drop two-station days whose modeled data fails the day-validity tests
#'
#' Sits between the day-window alignment and the fit. The alignment owns the
#' structural questions -- does this day fill its window, is its travel time
#' usable -- and this owns the data-quality ones: are the modeled values all
#' present, is depth positive. Without it \code{NA}s reach Stan and fail the
#' whole joint fit, naming neither column nor date (issue #475). Failing days
#' leave the alignment and are returned separately, so the caller can report
#' them rather than let them disappear.
#'
#' \code{\link{mm_filter_valid_days}} is deliberately not reused, though
#' \code{\link{mm_is_valid_day}} underneath it is: that function partitions
#' rows into one-station's overlapping diel windows, which are not two-station
#' days, and it filters a data.frame rather than an alignment.
#'
#' @section Why the tests run on the modeled frame: a two-station day is not a
#'   contiguous slice of \code{data} -- its upstream columns come from
#'   \code{aln$shift_idx}, rows one travel time earlier and routinely in the
#'   preceding day, while the rest come from \code{aln$keep}. Testing a day's
#'   own rows would both miss upstream \code{NA}s it depends on and blame
#'   others on the wrong date.
#'
#' @param data data.frame as validated by \code{\link{mm_validate_data}} for
#'   \code{\link{metab_bayes_2s}}, units optional. The full dataset, not a
#'   per-day slice.
#' @param aln an alignment as returned by \code{mm_align_2s}.
#' @param day_tests character vector of \code{\link{mm_is_valid_day}} tests to
#'   apply; see \code{\link{mm_day_tests_2s}} for which are applicable and why.
#'   \code{NULL} or \code{character(0)} skips testing entirely.
#' @return a list with \code{aln} (the input alignment, with failing dates
#'   removed and \code{n_days} updated) and \code{removed} (a data.frame of
#'   \code{date} and \code{errors}, one row per dropped date)
#' @keywords internal
mm_filter_valid_days_2s <- function(data, aln, day_tests=mm_day_tests_2s_default) {

  mm_check_day_tests_2s(day_tests)

  if(length(day_tests) == 0) {
    return(list(aln=aln, removed=mm_no_removed_days_2s()))
  }

  modeled <- mm_modeled_rows_2s(data, aln)

  # one test per day, over row indices grouped once rather than rescanning the
  # frame per date. ply_date/timestep_days are passed rather than re-derived
  unique_dates <- unique(aln$date)
  rows_by_date <- split(seq_along(aln$date), factor(aln$date, levels=unique_dates))
  validity <- lapply(unique_dates, function(dt) {
    mm_is_valid_day(
      modeled[rows_by_date[[as.character(dt)]], , drop=FALSE],
      day_tests=day_tests,
      ply_date=dt,
      timestep_days=aln$timestep_days)
  })

  invalid <- !vapply(validity, isTRUE, logical(1))
  if(!any(invalid)) {
    return(list(aln=aln, removed=mm_no_removed_days_2s()))
  }

  removed <- data.frame(
    date=unique_dates[invalid],
    errors=vapply(validity[invalid], paste0, character(1), collapse='; '),
    stringsAsFactors=FALSE)

  message(paste0(
    'dropping ', nrow(removed), ' day(s) that fail day_tests: ',
    paste(sprintf('%s (%s)', as.character(removed$date), removed$errors), collapse=', ')))

  # drop whole dates, so each surviving date keeps its full n_obs rows and
  # still occupies a contiguous block for the matrix pivot
  keep_rows <- !(aln$date %in% unique_dates[invalid])
  aln$keep <- aln$keep[keep_rows]
  aln$shift_idx <- aln$shift_idx[keep_rows]
  aln$date <- aln$date[keep_rows]
  aln$n_days <- length(unique(aln$date))

  list(aln=aln, removed=removed)
}

#' An empty removed-days data.frame
#'
#' The zero-row shape returned when nothing was dropped, so that the several
#' stages that can drop a day always combine without special-casing.
#'
#' @keywords internal
mm_no_removed_days_2s <- function() {
  data.frame(date=as.Date(character(0)), errors=character(0), stringsAsFactors=FALSE)
}
