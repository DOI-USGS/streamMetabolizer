#' Convert degrees to radians
#' 
#' @param degrees angle in degrees
#' @return angle in radians
to_radians <- function(degrees) {
  degrees * pi / 180
}

#' Convert radians to degrees
#' 
#' @param radians angle in radians
#' @return angle in degrees
to_degrees <- function(radians) {
  radians * 180 / pi
}

#' Calculate declination angle as in 
#' http://pveducation.org/pvcdrom/properties-of-sunlight/declination-angle which
#' cites http://www.sciencedirect.com/science/article/pii/0038092X69900474
#' 
#' @param jday The day of year as a number between 0 (Jan 1) and 364 (365 also
#'   OK for leap year)
#' @param format The format of both the input and the output. May be "degrees" 
#'   or "radians".
#' @return numeric value or vector, in the units specified by \code{format}, 
#'   indicating the declination angle corresponding to each value supplied in 
#'   \code{jday}.
#' @examples
#' decdf <- data.frame(jday=1:366, 
#'   dec=streamMetabolizer:::calc_declination_angle(1:366))
#' \dontrun{
#'   library(ggplot2)
#'   ggplot(decdf, aes(x=jday, y=dec)) + geom_line()
#' }
calc_declination_angle <- function(jday, format=c("degrees", "radians")) {
  format <- match.arg(format)
  declination.angle <- asin(sin(to_radians(23.439))*sin(to_radians((360/365)*(283+jday))))
  if(format == "degrees") {
    declination.angle <- to_degrees(declination.angle)
  }
  declination.angle
}

#' Calculate hour angle as in 
#' http://education.gsfc.nasa.gov/experimental/July61999siteupdate/inv99Project.Site/Pages/solar.insolation.html.
#' 
#' This is an approximation when hour is in clock time; should actually be given
#' in solar time
#' 
#' @param hour numeric value or vector. hour since [solar] midnight as number
#'   between 0 and 23.999
#' @param format The format of both the input and the output. May be "degrees" 
#'   or "radians".
#' @return numeric value or vector, in the units specified by \code{format},
#'   indicating the angle corresponding to each value supplied in \code{hour}.
#' @examples
#' hourdf <- data.frame(hour=c(0:12,12.5:23.5), 
#'   hragl=streamMetabolizer:::calc_hour_angle(c(0:12,12.5:23.5)))
#' \dontrun{
#'   library(ggplot2)
#'   ggplot(hourdf, aes(x=hour, y=hragl)) + 
#'     geom_hline(yintercept=0, color="gold") + geom_line()
#' }
calc_hour_angle <- function(hour, format=c("degrees", "radians")) {
  format <- match.arg(format)
  hour.angle <- (360/24)*(hour-12)
  if(format=="radians") hour.angle <- to_radians(hour.angle)
  hour.angle
}

#' Calculate zenith angle as in 
#' http://education.gsfc.nasa.gov/experimental/July61999siteupdate/inv99Project.Site/Pages/solar.insolation.html
#' 
#' @param latitude numeric value or vector indicating the site latitude in 
#'   decimal degrees (never radians or deg-min-sec, no matter what \code{format}
#'   is) between -90 (South Pole) and 90 (North Pole).
#' @param declination.angle numeric value or vector, in the units specified by 
#'   \code{format}, indicating the declination angle.
#' @param hour.angle numeric value or vector, in the units specified by 
#'   \code{format}, indicating the angle.
#' @param format The format of both the output. May be "degrees" or "radians".
#' @examples
#' zendf <- data.frame(
#'   lat=rep(c(0,20,40,60), each=24*4),
#'   jday=rep(rep(c(1,101,201,301), each=24), times=4), 
#'   hour=rep(c(0:12,13.5:23.5), times=4*4))
#' zendf <- transform(zendf, 
#'   dec=streamMetabolizer:::calc_declination_angle(jday),
#'   hragl=streamMetabolizer:::calc_hour_angle(hour))
#' zendf <- transform(zendf,
#'   zen=streamMetabolizer:::calc_zenith_angle(lat, dec, hragl))
#' \dontrun{
#'   library(ggplot2)
#'   ggplot(zendf, aes(x=jday, y=zen, color=factor(lat))) + 
#'     geom_line()
#' }
calc_zenith_angle <- function(latitude, declination.angle, hour.angle, format=c("degrees", "radians")) {
  format <- match.arg(format)
  latitude <- to_radians(latitude)
  if(format == "degrees") {
    declination.angle <- to_radians(declination.angle)
    hour.angle <- to_radians(hour.angle)
  }
  zenith.angle <- 
    acos(sin(latitude) * sin(declination.angle) + 
           cos(latitude) * cos(declination.angle) * cos(hour.angle))
  if(format == "degrees") {
    zenith.angle <- to_degrees(zenith.angle)
  }
  zenith.angle
}

#' Calculate mean solar time from clock time and longitude
#' 
#' @param date.time


#' Model solar insolation on a horizontal surface (W/m2 == J/s/m2) as in 
#' http://education.gsfc.nasa.gov/experimental/July61999siteupdate/inv99Project.Site/Pages/solar.insolation.html
#' 
#' @importFrom unitted u
#' @inheritParams calc_declination_angle
#' @inheritParams calc_hour_angle
#' @inheritParams calc_zenith_angle
#' @param max.insolation insolation rate at solar noon, W/m2 == J/s/m2. varies greatly with atmospheric conditions
#' @param attach.units logical. Should the returned vector be a unitted object?
#' @examples
#' insdf <- data.frame(
#'   lat=rep(c(0,20,40,60), each=24*4),
#'   jday=rep(rep(c(1,101,201,301), each=24), times=4), 
#'   hour=rep(c(0:12,13.5:23.5), times=4*4))
#' insdf <- transform(insdf, ins=calc_solar_insolation(jday, hour, lat))
#' \dontrun{
#'   library(ggplot2)
#'   ggplot(insdf, aes(color=factor(jday), y=ins, x=hour)) + 
#'     geom_line() + facet_wrap(~lat)
#' }
#' @export
calc_solar_insolation <- function(jday, hour, latitude, max.insolation=2326, format=c("degrees", "radians"), attach.units=FALSE) {
  format <- match.arg(format)
  declination.angle <- calc_declination_angle(jday, format=format)
  hour.angle <- calc_hour_angle(hour, format=format)
  zenith.angle <- calc_zenith_angle(latitude, declination.angle, hour.angle, format=format)
  if(format=="degrees") zenith.angle <- to_radians(zenith.angle)
  insolation <- max.insolation * cos(zenith.angle)
  insolation <- pmax(insolation, 0)

  if(attach.units) {
    u(insolation,"J s^-1 m^-2")
  } else {
    insolation
  }
}