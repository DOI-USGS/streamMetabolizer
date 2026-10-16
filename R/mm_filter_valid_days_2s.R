#' @include mm_lag_2s.R
NULL

#' Day-validity tests applicable to two-station models
#'
#' \code{mm_day_tests_2s_default} is the default; \code{mm_day_tests_2s_allowed}
#' bounds what a user may ask for -- a chosen subset of
#' \code{\link{mm_is_valid_day}}'s five tests, not the shared default
#' inherited wholesale.
#'
#' \code{full_day} and \code{even_timesteps} are excluded. \code{full_day}
#' checks for observations at the one-station day boundaries (4 AM to 4 AM by
#' default), which the two-station day boundary does not match.
#' \code{even_timesteps} is excluded because two-station methods do not
#' need regular timesteps: irregular spacing within or at the ends of a day
#' would not by itself invalidate a two-station day. (The current Stan model
#' does need each day on a complete regular grid, which alignment enforces
#' separately.)
#' \code{pos_discharge} is allowed but off by default; it has no effect unless
#' the data carry a \code{discharge} column, which two-station data currently
#' does not.
#'
#' @keywords internal
#' @name mm_day_tests_2s
NULL

#' @rdname mm_day_tests_2s
mm_day_tests_2s_default <- c('complete_data', 'pos_depth')

#' @rdname mm_day_tests_2s
mm_day_tests_2s_allowed <- c('complete_data', 'pos_discharge', 'pos_depth')

#' @section Minimum rows per day: \code{mm_min_n_obs_2s_default} is the
#'   default \code{min_n_obs}, the fewest observations a day must have: 20,
#'   which admits hourly or finer data. Data with fewer observations
#'   per day stop the run with an error rather than having days dropped; 0
#'   turns the check off.
#' @rdname mm_day_tests_2s
mm_min_n_obs_2s_default <- 20

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

#' Validate a two-station min_n_obs setting
#'
#' Called both when specs are created and again when they are used, since the
#' consumer can be called directly with a value that never passed through a
#' \code{specs} list.
#'
#' @param min_n_obs the value to check.
#' @keywords internal
#' @noRd
mm_check_min_n_obs_2s <- function(min_n_obs) {
  if(!is.numeric(min_n_obs) || length(min_n_obs) != 1 || !is.finite(min_n_obs) ||
     min_n_obs < 0 || min_n_obs != round(min_n_obs)) {
    stop('min_n_obs must be a single whole number of 0 or more', call.=FALSE)
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
#' rows by the one-station \code{day_start}/\code{day_end} boundaries, not
#' the two-station day boundary.
#'
#' @param data an \code{aligned_2s} data.frame (see
#'   \code{\link{mm_align_data_2s}}). Each row already holds the upstream
#'   values it is modeled from, so a day's rows are exactly the values tested.
#' @param day_tests character vector of \code{\link{mm_is_valid_day}} tests to
#'   apply; see \code{\link{mm_day_tests_2s}} for which are applicable and why.
#'   \code{NULL} or \code{character(0)} skips testing entirely.
#' @param min_n_obs the fewest observations each day must have. Every day of
#'   aligned data holds the same number of rows, so a shortfall is an error
#'   rather than a set of dropped days. Defaults to 20, which admits hourly
#'   or finer data; 0 turns the check off.
#' @return a list with \code{data} (the input with failing dates removed,
#'   still of class \code{aligned_2s}) and \code{removed} (a data.frame of
#'   \code{date} and \code{errors}, one row per dropped date)
#' @keywords internal
mm_filter_valid_days_2s <- function(data, day_tests=mm_day_tests_2s_default,
                                    min_n_obs=mm_min_n_obs_2s_default) {

  mm_check_day_tests_2s(day_tests)
  mm_check_min_n_obs_2s(min_n_obs)

  # aligned data holds the same number of rows on every day, so every day
  # passes or every day fails; stopping says why, where dropping them all
  # would not
  if(nrow(data) > 0) {
    n_per_day <- max(table(data$date))
    if(n_per_day < min_n_obs) {
      stop(mm_min_n_obs_msg_2s(data, n_per_day, min_n_obs), call.=FALSE)
    }
  }

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

# Error text for a min_n_obs shortfall. Days of user-aligned data can be
# shorter than 24 hours, so the row count alone does not reveal the timestep:
# the timestep is named, and finer data suggested, only when consecutive rows
# of one day give the timestep and the days span a full 24 hours at it
mm_min_n_obs_msg_2s <- function(data, n_per_day, min_n_obs) {
  timestep_min <- NA_real_
  if(n_per_day >= 2) {
    st <- as.numeric(v(data$solar.time))
    same_day <- v(data$date)[-1] == v(data$date)[-nrow(data)]
    timestep_min <- stats::median(diff(st)[same_day]) / 60
  }
  # a full day reaches within one timestep of 24 hours; 1-second tolerance
  full_days <- !is.na(timestep_min) &&
    n_per_day * timestep_min > 24 * 60 - timestep_min + 1 / 60

  head <- paste0(
    'data have ', n_per_day, ' observation', if(n_per_day == 1) '' else 's', ' per day')
  need <- paste0('; specs min_n_obs requires at least ', min_n_obs, '. ')
  if(full_days) {
    paste0(
      head, ' (a ', format(signif(timestep_min, 3)), '-minute timestep)', need,
      'Use hourly or finer data, or set min_n_obs=0 in specs().')
  } else {
    paste0(
      head, if(!is.na(timestep_min)) ', covering less than 24 hours', need,
      'Use data with at least ', min_n_obs, ' observations per day, or set ',
      'min_n_obs=0 in specs().')
  }
}

#' An empty removed-days data.frame
#'
#' The zero-row shape returned when nothing was dropped, so that the several
#' stages that can drop a day always combine without special-casing.
#'
#' @keywords internal
mm_no_removed_days_2s <- data.frame(date=as.Date(character(0)), errors=character(0), stringsAsFactors=FALSE)
