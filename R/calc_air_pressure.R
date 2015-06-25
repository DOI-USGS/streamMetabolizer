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
  
  # Bob's code:
  #   ###Estimate barometric pressure as a function of altitude and standard BP
  #   # since we will not correct for dail change in BP use 29.92 (=760 mm)
  #   ##this is based on the barometric formula
  #   ###temp is degC, alt is m, and bpst is in mm of Hg.  Temp is usually relative to a standard, 15 degC.  This is from Colt's book
  #   bpcalc<- function(bpst, alt, temp) {
  #     bpst*exp((-9.80665*0.0289644*alt)/(8.31447*(273.15+temp)))
  #   }
  #   ##call as
  #   bpcalc(bpst=760,alt=2400,temp=15)
  
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
