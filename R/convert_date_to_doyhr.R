#' Convert a date to a day of year (1-366) with decimal hours
#' 
#' Inspired by / copied from LakeMetabolizer date2doy
#' 
#' @param date A datetime object as POSIXct or POSIXt
#' @return Numeric value expressing the date as the number of days, with decimal
#'   hours, since 00:00 of December 31st of the preceding year (i.e., January
#'   1st at 00:01 is 1.01)
#' @export
#' @importFrom lubridate floor_date
convert_date_to_doyhr <- function(date) {
  1 + as.numeric(date - floor_date(date, unit="year"), units="days") # days_since_dec31
}

#' Convert a a day of year (1-366) with decimal hours to a date
#' 
#' @param doyhr Numeric value expressing the date as the number of days, with decimal
#'   hours, since 00:00 of January 1 of the year
#' @param year Numeric 4-digit year
#' @param tz The time zone to pass to as.POSIXct()
#' @param origin The origin to pass to as.POSIXct()
#' @param ... Other arguments to pass to as.POSIXct()
#' @return A datetime object as POSIXct
#'   
#' @export
convert_doyhr_to_date <- function(doyhr, year, tz="GMT", origin=as.POSIXct("1970-01-01 00:00:00",tz="GMT"), ...) {
  secs_since_jan1 <- (doyhr-1)*24*60*60
  as.POSIXct(sprintf("%d-01-01 00:00:00",year), tz=tz, origin=origin, ...) + secs_since_jan1
}