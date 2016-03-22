#' Determines if specified datetime is during the daytime
#' Returns T/F indicating whether a datetime occurs during the daytime (sunlight hours)
#'
#' @param datetimes Vector of dates as \code{POSIXct} or \code{POSIXlt} (see \code{\link{DateTimeClasses}}) format, but in SOLAR time
#' @param lat Single latitude value of site. South should be negative, north positive
#'
#' @return a boolean vector of same length as \code{datetimes} 
#'
#' @author
#' Luke A. Winslow
#' @seealso 
#' \link{calc_sun_rise_set}
#' @importFrom LakeMetabolizer is.day
#' @examples
#' calc_is_daytime(datetimes=as.POSIXct(paste('2013-03-31', c('1:00','11:00'))), lat=40.75)
#' @export
calc_is_daytime <- function(datetimes, lat) {
  LakeMetabolizer::is.day(datetimes, lat)
}