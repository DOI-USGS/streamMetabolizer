#' Convert degrees to radians
#' 
#' @param degrees angle in degrees
#' @return angle in radians
#' @importFrom unitted u is.unitted verify_units
to_radians <- function(degrees) {
  if(is.unitted(degrees)) {
    verify_units(degrees, "deg")
    degrees * u(pi / 180, "rad deg^-1")
  } else {
    degrees * pi / 180
  }
}

#' Convert radians to degrees
#' 
#' @param radians angle in radians
#' @return angle in degrees
#' @importFrom unitted u is.unitted verify_units
to_degrees <- function(radians) {
  if(is.unitted(radians)) {
    verify_units(radians, "rad")
    radians * u(180 / pi, "deg rad^-1")
  } else {
    radians * 180 / pi
  }
}

#' Calculate declination angle as in Yard et al. (2005)
#' 
#' @param jday The day of year as a number between 0 (Jan 1) and 364 (365 also 
#'   OK for leap year)
#' @param format The format of both the input and the output. May be "degrees" 
#'   or "radians".
#' @return numeric value or vector, in the units specified by \code{format}, 
#'   indicating the declination angle corresponding to each value supplied in 
#'   \code{jday}.
#' @importFrom unitted u is.unitted
#' @examples
#' decdf <- data.frame(jday=1:366, 
#'   dec=streamMetabolizer:::calc_declination_angle(1:366))
#' \dontrun{
#' library(ggplot2)
#' ggplot(unitted::v(decdf), aes(x=jday, y=dec)) + geom_line()
#' }
#' @references Yard, Michael D., Glenn E. Bennett, Steve N. Mietz, Lewis G. 
#'   Coggins Jr., Lawrence E. Stevens, Susan Hueftle, and Dean W. Blinn. 
#'   \emph{Influence of Topographic Complexity on Solar Insolation Estimates for
#'   the Colorado River, Grand Canyon, AZ.} Ecological Modelling 183, no. 2-3
#'   (April 25, 2005): 157-72. doi:10.1016/j.ecolmodel.2004.07.027.
calc_declination_angle <- function(jday, format=c("degrees", "radians")) {
  format <- match.arg(format)
  declination.angle <- u(23.439,"deg")*sin(to_radians((360/365)*(283+v(jday))))
  if(format == "radians") {
    declination.angle <- to_radians(declination.angle)
  }
  if(is.unitted(jday)) declination.angle <- u(declination.angle, substr(format, 1, 3))
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
#' @importFrom unitted u is.unitted
#' @examples
#' hourdf <- data.frame(hour=c(0:12,12.5:23.5), 
#'   hragl=streamMetabolizer:::calc_hour_angle(c(0:12,12.5:23.5)))
#' \dontrun{
#' library(ggplot2)
#' ggplot(hourdf, aes(x=hour, y=hragl)) + 
#'   geom_hline(yintercept=0, color="gold") + geom_line()
#' }
calc_hour_angle <- function(hour, format=c("degrees", "radians")) {
  format <- match.arg(format)
  hour.angle <- (360/24)*(hour-12)
  if(is.unitted(hour)) hour.angle <- u(hour.angle, "deg")
  if(format=="radians") hour.angle <- to_radians(hour.angle)
  if(is.unitted(hour)) hour.angle <- u(hour.angle, substr(format, 1, 3))
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
#' @importFrom unitted u is.unitted get_units
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
#' library(ggplot2)
#' ggplot(unitted::v(zendf), aes(x=hour, y=zen, color=jday, group=jday)) + 
#'   geom_line() + facet_wrap(~lat) + 
#'   ggtitle('zenith angles by latitude (panels) and day of year (colors)')
#' }
calc_zenith_angle <- function(latitude, declination.angle, hour.angle, format=c("degrees", "radians")) {
  format <- match.arg(format)
  if(is.unitted(latitude)) {
    # convert latitude units to plain 'deg' as needed
    latitude <- 
      switch(
        get_units(latitude),
        "deg"=latitude,
        "degN"=u(latitude, "deg"),
        "degS"=u(-latitude, "deg"))
  } 
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


#' Model solar insolation on a horizontal surface (W/m2 == J/s/m2) as in 
#' http://education.gsfc.nasa.gov/experimental/July61999siteupdate/inv99Project.Site/Pages/solar.insolation.html
#' 
#' @importFrom unitted u
#' @param app.solar.time POSIXct vector of date-time values in apparent solar time,
#'   e.g., as returned by \code{convert_UTC_to_solartime(...,
#'   time.type="apparent solar")}
#' @inheritParams calc_declination_angle
#' @inheritParams calc_hour_angle
#' @inheritParams calc_zenith_angle
#' @param max.insolation insolation rate at solar noon, W/m2 == J/s/m2. varies
#'   greatly with atmospheric conditions
#' @param attach.units logical. Should the returned vector be a unitted object?
#' @examples
#' insdf <- data.frame(
#'   lat=rep(c(0,20,40,60), each=24*4),
#'   jday=rep(rep(c(1,101,201,301), each=24), times=4), 
#'   hour=rep(c(0:12,13.5:23.5), times=4*4))
#' insdf <- transform(insdf, datetime=convert_doyhr_to_date(jday + hour/24, year=2004))
#' insdf <- transform(insdf, ins=calc_solar_insolation(datetime, lat))
#' \dontrun{
#' library(ggplot2)
#' ggplot(insdf, aes(color=factor(jday), y=ins, x=hour)) + 
#'   geom_line() + facet_wrap(~lat) + 
#'   ggtitle('solar insolation by latitude (panels) and day of year (colors)')
#' }
#' @export
calc_solar_insolation <- function(
  app.solar.time, latitude, max.insolation=convert_PAR_to_SW(2326), 
  format=c("degrees", "radians"), attach.units=is.unitted(app.solar.time)) {
  
  format <- match.arg(format)
  jday <- floor(convert_date_to_doyhr(app.solar.time)) - 1
  hour <- (convert_date_to_doyhr(app.solar.time) %% 1) * 24
  declination.angle <- calc_declination_angle(jday, format=format)
  hour.angle <- calc_hour_angle(hour, format=format)
  zenith.angle <- calc_zenith_angle(latitude, declination.angle, hour.angle, format=format)
  if(format=="degrees") zenith.angle <- to_radians(zenith.angle)
  insolation <- max.insolation * cos(zenith.angle)
  insolation <- pmax(insolation, 0)

  if(attach.units) {
    u(insolation,"W m^-2")
  } else {
    v(insolation)
  }
}
