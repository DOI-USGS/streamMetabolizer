#' @include mm_lag_2s.R mm_modeled_rows_2s.R mm_filter_valid_days_2s.R
NULL

#' Align two-station data onto its modeled rows and 24-hour days
#'
#' Resolves each downstream observation's upstream partner one travel time
#' earlier and returns one row per modeled observation, labeled with the
#' 24-hour day (starting at \code{day_start_hour}) it belongs to. Rows without upstream lead-in, days whose
#' travel time exceeds \code{max_travel_time_days}, and days that do not fill
#' their 24-hour window are dropped with a message and recorded in the
#' result's \code{removed} attribute.
#'
#' Gaps are not filled here. When data have gaps, call
#' \code{\link{mm_fill_gaps_2s}} first.
#'
#' @section Two-station day window: Two-station days run 24 hours from
#'   \code{day_start_hour}, 06:00-06:00 by default. The one-station models'
#'   default day is also 24 hours, 4 AM to 4 AM (\code{day_start=4},
#'   \code{day_end=28}); the two differ only in where the day boundary falls.
#'   Either way, GPP, ER, and gas exchange are held constant over one day's
#'   data. Days that do not fill the window -- at the edges of a dataset whose
#'   bounds don't fall on the day start hour, or
#'   where observations are missing and the gap was left unfilled (see
#'   \code{\link{mm_fill_gaps_2s}}) -- are dropped with a message and
#'   reported as invalid days in the results of
#'   \code{\link{metab_bayes_2s}}.
#'
#' @section Lead-in and travel-time ceiling: \code{data$travel.time} (the
#'   reach travel time between stations, in days) must be strictly positive,
#'   and at least one row must have enough preceding observations of upstream
#'   DO to cover its own travel time. Rows that lack that lead-in -- whose
#'   look-back one travel time earlier falls before the start of \code{data}
#'   or inside a gap -- are not an error: they serve as lead-in only,
#'   supplying upstream DO for later rows without appearing in the result
#'   themselves.
#'
#'   Travel time is also subject to a ceiling; see the
#'   \code{max_travel_time_days} argument.
#'
#' @param data data.frame, units optional, with the columns
#'   \code{\link{metab_bayes_2s}} expects (see \code{\link{mm_format_data_2s}}),
#'   sorted ascending by \code{solar.time}, with \code{solar.time} on a single
#'   regular timestep grid, as in \code{\link{mm_format_data_2s}} output.
#'   \code{travel.time} may not be \code{NA}.
#' @param max_travel_time_days the travel-time ceiling, in days. Defaults to
#'   0.42 days (10 hours); values above 0.5 days (12 hours) are rejected.
#'   Beyond the ceiling, a day's upstream parcel almost certainly originates
#'   before the day's own start, where the light it experienced no
#'   longer has a well-defined day total to be a proportion of. Days whose
#'   longest travel time exceeds the ceiling are dropped with a message rather
#'   than treated as an error, since the remaining days are unaffected. A
#'   travel time far above the ceiling usually means the column was supplied
#'   in the wrong units -- days are expected, not minutes or hours.
#' @param day_start_hour hour of the day, in [0, 24), at which each 24-hour
#'   two-station day begins. Defaults to 6 (06:00-06:00). Must match the
#'   value given to \code{\link{mm_format_data_2s}}, whose light proportions
#'   depend on it; a different value recorded on \code{data} is an error.
#' @return a data.frame of class \code{aligned_2s}, unitless, with a
#'   \code{date} column followed by the \code{\link{metab_bayes_2s}} data
#'   columns in their standard order, and attributes \code{removed} (a
#'   data.frame of \code{date} and \code{errors} for each dropped day),
#'   \code{max_travel_time_days}, and \code{day_start_hour}.
#' @examples
#' # four 06:00-06:00 days (the default) of the example data. The first day's
#' # opening rows serve only as lead-in, so that day does not fill its window
#' # and is dropped
#' dat <- two_station_example[
#'   unitted::v(two_station_example$solar.time) < as.POSIXct('2008-03-16', tz='UTC'), ]
#' aligned <- mm_align_data_2s(dat)
#' table(aligned$date)
#' attr(aligned, 'removed')
#' @importFrom unitted v
#' @export
mm_align_data_2s <- function(data, max_travel_time_days=mm_max_travel_time_default,
                             day_start_hour=mm_day_start_2s) {

  mm_check_day_start_hour_2s(day_start_hour)
  # checked before validation, which drops the attribute
  formatted_hour <- attr(data, 'day_start_hour')
  if(!is.null(formatted_hour)) mm_check_day_start_hour_2s(formatted_hour)
  if(!is.null(formatted_hour) && formatted_hour != day_start_hour) {
    stop(paste0(
      'data were formatted with day_start_hour=', formatted_hour, ' but are being aligned with ',
      'day_start_hour=', day_start_hour, '; light proportions depend on the day start hour, ',
      'so use the same value for both'), call.=FALSE)
  }
  mm_stop_if_na_travel_time_2s(
    data, 'fill it (mm_fill_gaps_2s() bridges short gaps) or remove those rows before aligning',
    day_start_hour)

  data <- mm_validate_data(data, NULL, 'metab_bayes_2s')$data
  mm_validate_data_2station(data)

  alignment <- mm_align_2s(
    v(data), max_travel_time_days=max_travel_time_days, day_start_hour=day_start_hour)
  frame <- data.frame(date=alignment$date, mm_modeled_rows_2s(data, alignment))

  new_aligned_2s(
    frame, removed=alignment$removed, max_travel_time_days=max_travel_time_days,
    day_start_hour=day_start_hour)
}

#' Mark a user-aligned data.frame as aligned two-station data
#'
#' For data aligned some other way than \code{\link{mm_align_data_2s}}: each
#' row must already hold its downstream observation and the upstream values
#' paired with it. The data.frame is validated as aligned data, and its dropped
#' days and travel-time ceiling are recorded as unknown.
#'
#' @param data data.frame, units optional, with the columns
#'   \code{\link{metab_bayes_2s}} expects and optionally a \code{date} column
#'   of day labels for 24-hour days starting at \code{day_start_hour}. If
#'   \code{date} is absent it is added.
#' @param day_start_hour hour of the day, in [0, 24), at which each 24-hour
#'   two-station day begins. Defaults to 6 (06:00-06:00). The \code{light}
#'   proportions must have been computed with the same day start hour.
#' @return a data.frame of class \code{aligned_2s}, unitless, with a
#'   \code{date} column first and a \code{day_start_hour} attribute.
#' @examples
#' dat <- two_station_example[
#'   unitted::v(two_station_example$solar.time) < as.POSIXct('2008-03-16', tz='UTC'), ]
#' aligned <- suppressMessages(mm_align_data_2s(dat))
#'
#' # a plain data.frame of aligned rows, with or without its date column,
#' # marks back to the same rows
#' plain <- as.data.frame(aligned)
#' identical(as.data.frame(mm_as_aligned_2s(plain)), plain)
#' identical(as.data.frame(mm_as_aligned_2s(plain[names(plain) != 'date'])), plain)
#'
#' # the days dropped during alignment are not known for a marked data.frame
#' attr(mm_as_aligned_2s(plain), 'removed')
#' @importFrom unitted v
#' @export
mm_as_aligned_2s <- function(data, day_start_hour=mm_day_start_2s) {

  mm_check_day_start_hour_2s(day_start_hour)
  df <- as.data.frame(v(data))

  mm_stop_if_na_travel_time_2s(
    df, 'fill it or remove those days before marking the data as aligned', day_start_hour)

  # the timestamp test is left out because it accepts only one of
  # solar.time/date; the aligned-data checks below cover both columns instead
  df <- mm_validate_data(
    df, NULL, 'metab_bayes_2s', data_tests=c('missing_cols','extra_cols','units'))$data
  df <- as.data.frame(v(df))

  if(!('date' %in% names(df))) {
    if(!lubridate::is.POSIXct(df$solar.time)) {
      stop("expecting 'solar.time' to be of class 'POSIXct'", call.=FALSE)
    }
    df <- data.frame(date=mm_date_2s(df$solar.time, day_start_hour), df)
  } else if(lubridate::is.POSIXct(df$solar.time) && lubridate::is.Date(df$date) &&
            !anyNA(df$solar.time) && !anyNA(df$date) &&
            !isTRUE(all(df$date == mm_date_2s(df$solar.time, day_start_hour)))) {
    # checked here as well as by the validator below, so the message can point
    # to this function's argument; malformed columns are left to the validator
    stop(paste0(
      "'date' does not match 24-hour days starting at hour ", day_start_hour,
      "; if the dates were made with a different day start hour, pass that ",
      "hour as day_start_hour"), call.=FALSE)
  }
  aligned <- new_aligned_2s(df, day_start_hour=day_start_hour)

  mm_validate_data_2station(aligned)
  aligned
}

# Stop if travel.time has NAs, naming the affected days and ending with `fix`,
# the caller's advice on what to do. Without this, an NA travel.time fails the
# shared validator's positivity test with a base-R error that names neither
# the column nor the day. A frame without travel.time is left to the
# missing-columns check. day_start_hour sets the day labels in the message;
# when it is NULL (unknown), rows are named instead.
mm_stop_if_na_travel_time_2s <- function(data, fix, day_start_hour) {
  if(!('travel.time' %in% names(data))) return(invisible(NULL))
  na_rows <- is.na(v(data[['travel.time']]))
  if(!any(na_rows)) return(invisible(NULL))
  solar_time <- v(data[['solar.time']])
  where <- if(lubridate::is.POSIXct(solar_time) && !is.null(day_start_hour)) {
    paste('on', paste(unique(mm_date_2s(solar_time[na_rows], day_start_hour)), collapse=', '))
  } else {
    paste('in rows', paste(which(na_rows), collapse=', '))
  }
  stop(paste0('travel.time is NA ', where, '; ', fix), call.=FALSE)
}

# Set the aligned_2s class and attributes on a data.frame. A NULL removed or
# max_travel_time_days leaves that attribute absent, meaning unknown. A NULL
# day_start_hour also leaves it absent, which the aligned-data validator
# rejects: the day labels can't be checked without it.
new_aligned_2s <- function(df, removed=NULL, max_travel_time_days=NULL, day_start_hour=NULL) {
  attr(df, 'removed') <- removed
  attr(df, 'max_travel_time_days') <- max_travel_time_days
  attr(df, 'day_start_hour') <- day_start_hour
  class(df) <- c('aligned_2s', 'data.frame')
  df
}

#' @export
#' @noRd
`[.aligned_2s` <- function(x, ...) {
  out <- NextMethod()
  if(!is.data.frame(out)) return(out)
  new_aligned_2s(
    out, removed=attr(x, 'removed'), max_travel_time_days=attr(x, 'max_travel_time_days'),
    day_start_hour=attr(x, 'day_start_hour'))
}

# removed is always NULL on the result: the dropped days of separately aligned
# pieces can't be reconciled into one record. The ceiling carries over only
# when every piece shares it. The day start hour carries over when every
# piece that has one agrees: a piece without one (e.g. a plain data.frame)
# doesn't conflict, since validation checks the hour against every date
#' @export
#' @noRd
rbind.aligned_2s <- function(..., deparse.level=1, make.row.names=TRUE, stringsAsFactors=FALSE) {
  pieces <- Filter(Negate(is.null), list(...))
  ceilings <- lapply(pieces, attr, 'max_travel_time_days')
  ceil <- if(!any(vapply(ceilings, is.null, logical(1))) &&
                all(vapply(ceilings, identical, logical(1), ceilings[[1]]))) ceilings[[1]] else NULL
  hours <- unique(Filter(Negate(is.null), lapply(pieces, attr, 'day_start_hour')))
  start_hour <- if(length(hours) == 1) hours[[1]] else NULL

  plain <- lapply(pieces, as.data.frame)
  out <- do.call(rbind.data.frame, c(plain, list(
    deparse.level=deparse.level, make.row.names=make.row.names, stringsAsFactors=stringsAsFactors)))
  rownames(out) <- NULL
  new_aligned_2s(out, removed=NULL, max_travel_time_days=ceil, day_start_hour=start_hour)
}

# a plain data.frame carries no alignment record, so the attributes go with
# the class
#' @export
#' @noRd
as.data.frame.aligned_2s <- function(x, ...) {
  attr(x, 'removed') <- NULL
  attr(x, 'max_travel_time_days') <- NULL
  attr(x, 'day_start_hour') <- NULL
  class(x) <- 'data.frame'
  x
}

# dplyr verbs rebuild their result from the first input's attributes. That is
# right when every row came from that input (filter, arrange, mutate), but
# bind_rows() would otherwise stamp the first piece's removed days and ceiling
# onto rows aligned separately, so those become unknown instead. The day
# start hour is kept either way: the aligned-data validator checks it against
# every row's date, so a wrong one can't pass silently
#' @importFrom dplyr dplyr_reconstruct
#' @export
#' @noRd
dplyr_reconstruct.aligned_2s <- function(data, template) {
  out <- NextMethod()
  from_template <-
    'solar.time' %in% names(data) && 'solar.time' %in% names(template) &&
    !anyDuplicated(data$solar.time) && all(data$solar.time %in% template$solar.time)
  if(!from_template) {
    attr(out, 'removed') <- NULL
    attr(out, 'max_travel_time_days') <- NULL
  }
  out
}
