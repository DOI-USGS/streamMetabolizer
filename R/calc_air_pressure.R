#' Calculates the average air pressure for a site
#' 
#' Corrects air pressure for air temperature and elevation
#' 
#' @param temp.air air temperature in degrees C. Default is 15 degC.
#' @param elevation the site elevation above sea level in m. Default is the
#'   rough mean elevation of the USA at 2500 ft (from
#'   http://www.infoplease.com/ipa/A0001792.html)
#' @param attach.units logical. Should the returned vector be a unitted object?
#' @return a numeric vector of barometric pressures in mb, with units attached 
#'   if requested.
#'   
#' @importFrom unitted u v get_units verify_units is.unitted
#' @examples
#' calc_air_pressure(15, 762)
#' calc_air_pressure(15, 100)
#' @export
calc_air_pressure <- function(temp.air=u(15, "degC"), elevation=u(762, "m"), attach.units=FALSE) {
  
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
  
  # assume units if not provided
  if(!is.unitted(temp.air)) temp.air <- u(temp.air, "degC")
  if(!is.unitted(elevation)) elevation <- u(elevation, "m")

  # check units
  verify_units(temp.air, "degC")
  verify_units(elevation, "m")
  
  # compute pressure. eqn also at https://en.wikipedia.org/wiki/Barometric_formula
  Pb <- u(760, "mmHg") # standard pressure
  g0 <- u(9.80665, "m s^-2") # gravitational acceleration
  M <- u(0.0289644, "kg mol^-1") # molar mass of Earth's air
  Rst <- u(8.31447, "N m mol^-1 K^-1") * u(1, "kg m s^-2 N^-1") # universal gas constant for air: 8.31432 N*m /(mol*K)
  Ta <- u(273.15, "K") + temp.air*u(1, "K degC^-1") # actual temperature in Kelvins
  baro <- Pb * exp((-1 * g0 * M * elevation)/(Rst * Ta)) * u(1.33322368, "mb mmHg^-1")
  
  # return
  if(attach.units) baro else v(baro)
}
