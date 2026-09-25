#' @include mm_lag_2s.R
NULL

#' Day-validity tests applicable to two-station models
#'
#' \code{mm_day_tests_2s_default} is the default; \code{mm_day_tests_2s_allowed}
#' bounds what a user may ask for -- a chosen subset of
#' \code{\link{mm_is_valid_day}}'s five tests, not the shared default
#' inherited wholesale.
#'
#' \code{full_day} and \code{even_timesteps} are excluded because two-station
#' methods don't need regular timesteps: \code{full_day} requires observations
#' at both ends of the day window (24 hours by default), and
#' \code{even_timesteps} requires evenly spaced timesteps, which two-station tolerates
#' gaps in on purpose (gap tolerance is a known benefit of two-station methods!). \code{pos_discharge} is allowed but off by default, a
#' no-op until two-station data carries a \code{discharge} column (issue
#' #475).
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
#' Sits between alignment and the fit: alignment owns the structural
#' questions (does this day fill its window, is its travel time usable), and
#' this owns the data-quality ones (are the modeled values all present, is
#' depth positive). Failing days are removed from the data and returned
#' separately, so the caller can report them rather than let them disappear.
#'
#' \code{\link{mm_filter_valid_days}} is deliberately not reused, though
#' \code{\link{mm_is_valid_day}} underneath it is: that function partitions
#' rows into one-station's overlapping diel windows, which are not
#' two-station days.
#'
#' @param data an \code{aligned_2s} data.frame (see
#'   \code{\link{mm_align_data_2s}}). Each row already holds the upstream
#'   values it is modeled from, so a day's rows are exactly the values tested.
#' @param day_tests character vector of \code{\link{mm_is_valid_day}} tests to
#'   apply; see \code{\link{mm_day_tests_2s}} for which are applicable and why.
#'   \code{NULL} or \code{character(0)} skips testing entirely.
#' @return a list with \code{data} (the input with failing dates removed,
#'   still of class \code{aligned_2s}) and \code{removed} (a data.frame of
#'   \code{date} and \code{errors}, one row per dropped date)
#' @keywords internal
mm_filter_valid_days_2s <- function(data, day_tests=mm_day_tests_2s_default) {

  mm_check_day_tests_2s(day_tests)

  if(length(day_tests) == 0) {
    return(list(data=data, removed=mm_no_removed_days_2s))
  }

  # one test per day, over row indices grouped once rather than rescanning the
  # frame per date. ply_date is passed rather than re-derived from solar.time
  unique_dates <- unique(data$date)
  rows_by_date <- split(seq_len(nrow(data)), factor(data$date, levels=unique_dates))
  validity <- lapply(unique_dates, function(dt) {
    mm_is_valid_day(
      data[rows_by_date[[as.character(dt)]], , drop=FALSE],
      day_tests=day_tests,
      ply_date=dt)
  })

  invalid <- !vapply(validity, isTRUE, logical(1))
  if(!any(invalid)) {
    return(list(data=data, removed=mm_no_removed_days_2s))
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
  list(data=data[!(data$date %in% removed$date), , drop=FALSE], removed=removed)
}

#' An empty removed-days data.frame
#'
#' The zero-row shape returned when nothing was dropped, so that the several
#' stages that can drop a day always combine without special-casing.
#'
#' @keywords internal
mm_no_removed_days_2s <- data.frame(date=as.Date(character(0)), errors=character(0), stringsAsFactors=FALSE)
