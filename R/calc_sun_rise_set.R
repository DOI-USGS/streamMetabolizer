#' Calculates the time of sunrise and sunset
#' 
#' Calculates the time of sunrise and sunset based on latitude and date.
#'
#' @param datetimes Vector of dates as \code{POSIXct} or \code{POSIXlt} (see \code{\link{DateTimeClasses}}) format, but in SOLAR time
#' @param lat Single latitude value of site. South should be negative, north positive
#' @return data.frame of 
#' @importFrom LakeMetabolizer sun.rise.set
#' @keywords internal
#' @examples
#' calc_sun_rise_set(lat=40.75,datetimes=as.POSIXlt('2013-03-31'))
#' @author
#' Luke A. Winslow
#' @seealso 
#' \link{is_daytime}
#' @export
calc_sun_rise_set <- function(datetimes, lat) {
  sun.rise.set <- LakeMetabolizer::sun.rise.set(datetimes, lat)
  # reformat
  data.frame(sunrise=sun.rise.set[,1], sunset=sun.rise.set[,2])
}



