#' @include mm_lag_2s.R mm_modeled_rows_2s.R mm_filter_valid_days_2s.R
NULL

#' Align two-station data onto its modeled rows and 06:00 days
#'
#' Resolves each downstream observation's upstream partner one travel time
#' earlier and returns one row per modeled observation, labeled with the
#' 06:00-06:00 day it belongs to. Rows without upstream lead-in, days whose
#' travel time exceeds \code{max_travel_time_days}, and days that do not fill
#' their 24-hour window are dropped with a message and recorded in the
#' result's \code{removed} attribute.
#'
#' Gaps are not filled here. When data have gaps, call
#' \code{\link{mm_fill_gaps_2s}} first.
#'
#' @section Two-station day window: Two-station days run 06:00-06:00 (24
#'   hours). The one-station models' default day is also 24 hours, 4 AM to 4
#'   AM (\code{day_start=4}, \code{day_end=28}); the two differ only in where
#'   the day boundary falls. Either way, GPP, ER, and gas exchange are held
#'   constant over one day's data. Days that do not fill the
#'   window -- at the edges of a dataset whose bounds don't fall on 06:00, or
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
#'   before the day's own 06:00 start, where the light it experienced no
#'   longer has a well-defined day total to be a proportion of. Days whose
#'   longest travel time exceeds the ceiling are dropped with a message rather
#'   than treated as an error, since the remaining days are unaffected. A
#'   travel time far above the ceiling usually means the column was supplied
#'   in the wrong units -- days are expected, not minutes or hours.
#' @return a data.frame of class \code{aligned_2s}, unitless, with a
#'   \code{date} column followed by the \code{\link{metab_bayes_2s}} data
#'   columns in their standard order, and attributes \code{removed} (a
#'   data.frame of \code{date} and \code{errors} for each dropped day) and
#'   \code{max_travel_time_days}.
#' @examples
#' # four 06:00-06:00 days of the example data. The first day's opening rows
#' # serve only as lead-in, so that day does not fill its window and is dropped
#' dat <- two_station_example[
#'   unitted::v(two_station_example$solar.time) < as.POSIXct('2008-03-16', tz='UTC'), ]
#' aligned <- mm_align_data_2s(dat)
#' table(aligned$date)
#' attr(aligned, 'removed')
#' @importFrom unitted v
#' @export
mm_align_data_2s <- function(data, max_travel_time_days=mm_max_travel_time_default) {

  mm_stop_if_na_travel_time_2s(
    data, 'fill it (mm_fill_gaps_2s() bridges short gaps) or remove those rows before aligning')

  data <- mm_validate_data(data, NULL, 'metab_bayes_2s')$data
  mm_validate_data_2station(data)

  alignment <- mm_align_2s(v(data), max_travel_time_days=max_travel_time_days)
  frame <- data.frame(date=alignment$date, mm_modeled_rows_2s(data, alignment))

  new_aligned_2s(frame, removed=alignment$removed, max_travel_time_days=max_travel_time_days)
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
#'   of 06:00-06:00 day labels. If \code{date} is absent it is added.
#' @return a data.frame of class \code{aligned_2s}, unitless, with a
#'   \code{date} column first.
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
mm_as_aligned_2s <- function(data) {

  df <- as.data.frame(v(data))

  mm_stop_if_na_travel_time_2s(df, 'fill it or remove those days before marking the data as aligned')

  # the timestamp test is left out because it accepts only one of
  # solar.time/date; the aligned-data checks below cover both columns instead
  df <- mm_validate_data(
    df, NULL, 'metab_bayes_2s', data_tests=c('missing_cols','extra_cols','units'))$data
  df <- as.data.frame(v(df))

  if(!('date' %in% names(df))) {
    if(!lubridate::is.POSIXct(df$solar.time)) {
      stop("expecting 'solar.time' to be of class 'POSIXct'", call.=FALSE)
    }
    df <- data.frame(date=mm_date_2s(df$solar.time), df)
  }
  aligned <- new_aligned_2s(df)

  mm_validate_data_2station(aligned)
  aligned
}

# Stop if travel.time has NAs, naming the affected days and ending with `fix`,
# the caller's advice on what to do. Without this, an NA travel.time fails the
# shared validator's positivity test with a base-R error that names neither
# the column nor the day. A frame without travel.time is left to the
# missing-columns check.
mm_stop_if_na_travel_time_2s <- function(data, fix) {
  if(!('travel.time' %in% names(data))) return(invisible(NULL))
  na_rows <- is.na(v(data[['travel.time']]))
  if(!any(na_rows)) return(invisible(NULL))
  solar_time <- v(data[['solar.time']])
  where <- if(lubridate::is.POSIXct(solar_time)) {
    paste('on', paste(unique(mm_date_2s(solar_time[na_rows])), collapse=', '))
  } else {
    paste('in rows', paste(which(na_rows), collapse=', '))
  }
  stop(paste0('travel.time is NA ', where, '; ', fix), call.=FALSE)
}

# Set the aligned_2s class and attributes on a data.frame. A NULL removed or
# max_travel_time_days leaves that attribute absent, meaning unknown.
new_aligned_2s <- function(df, removed=NULL, max_travel_time_days=NULL) {
  attr(df, 'removed') <- removed
  attr(df, 'max_travel_time_days') <- max_travel_time_days
  class(df) <- c('aligned_2s', 'data.frame')
  df
}

#' @export
#' @noRd
`[.aligned_2s` <- function(x, ...) {
  out <- NextMethod()
  if(!is.data.frame(out)) return(out)
  new_aligned_2s(out, removed=attr(x, 'removed'), max_travel_time_days=attr(x, 'max_travel_time_days'))
}

# removed is always NULL on the result: the dropped days of separately aligned
# pieces can't be reconciled into one record
#' @export
#' @noRd
rbind.aligned_2s <- function(..., deparse.level=1, make.row.names=TRUE, stringsAsFactors=FALSE) {
  pieces <- Filter(Negate(is.null), list(...))
  ceilings <- lapply(pieces, attr, 'max_travel_time_days')
  ceil <- if(!any(vapply(ceilings, is.null, logical(1))) &&
                all(vapply(ceilings, identical, logical(1), ceilings[[1]]))) ceilings[[1]] else NULL

  plain <- lapply(pieces, as.data.frame)
  out <- do.call(rbind.data.frame, c(plain, list(
    deparse.level=deparse.level, make.row.names=make.row.names, stringsAsFactors=stringsAsFactors)))
  rownames(out) <- NULL
  new_aligned_2s(out, removed=NULL, max_travel_time_days=ceil)
}

# a plain data.frame carries no alignment record, so the attributes go with
# the class
#' @export
#' @noRd
as.data.frame.aligned_2s <- function(x, ...) {
  attr(x, 'removed') <- NULL
  attr(x, 'max_travel_time_days') <- NULL
  class(x) <- 'data.frame'
  x
}

# dplyr verbs rebuild their result from the first input's attributes. That is
# right when every row came from that input (filter, arrange, mutate), but
# bind_rows() would otherwise stamp the first piece's removed days and ceiling
# onto rows aligned separately, so those become unknown instead
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
