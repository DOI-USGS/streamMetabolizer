#' Calculates the average air pressure for a site
#' 
#' Will eventually correct for site elevation, but for now just returns standard
#' pressure
#' 
#' @param elevation the site elevation above sea level in m
#' @param attach.units logical. Should the returned vector be a unitted object?
#' @return a numeric vector of barometric pressures in Pa, with units attached 
#'   if any of the input vectors are unitted.
#'   
#' @importFrom unitted u v get_units verify_units is.unitted
#' @examples
#' library(unitted)
#' calc_air_pressure() # no units checking if no units provided
#' @export
calc_air_pressure <- function(elevation, attach.units=FALSE) {
  if(missing(elevation)) {
    # the mean elevation of the USA is 2500 ft: http://www.infoplease.com/ipa/A0001792.html
    elevation <- u(2500, "ft") * u(0.3048, "m ft^-1")
    # which comes out to very roughly 95000 Pa: https://en.wikipedia.org/wiki/Atmospheric_pressure#Mean_sea_level_pressure
    baro <- u(95000, "Pa")
  } else {
    if(get_units(elevation) != "m") stop("wrong elevation units; need m")
    stop("elevation correction isn't yet implemented")
  }
  if(attach.units) baro else v(baro)
}
