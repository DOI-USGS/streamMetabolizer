#' @title calculate a vector of dissolved oxygen deficit
#' @description Creates a DO.deficit vector for input into various metabolism
#'   models (see...)
#' @param DO.obs a numeric vector of dissolved oxygen concentration
#'   observations, mgO2 L^-1, or a \linkS4class{unitted} object of dissolved
#'   oxygen concentrations.
#' @param temp.water a numeric vector of water temperature in degrees Celsius, 
#'   or a \linkS4class{unitted} object of water temperatures.
#' @param pressure.air barometric pressure in millibars, or a
#'   \linkS4class{unitted} object of barometric pressure.
#' @param salinity.water a numeric vector of salinity in PSU, or a
#'   \linkS4class{unitted} object of salinity. Defaults to zero. Length must be
#'   one or equal to length of \code{temp.water}.
#' @param ... additional parameters passed to
#'   \code{\link[LakeMetabolizer]{o2.at.sat.base}}
#' @return a vector of DO.deficit values
#' @examples
#' calc_DO_deficit(DO.obs=7, temp.water=25, pressure.air=900, salinity.water=2.43)
#' 
#' library(unitted)
#' calc_DO_deficit(
#'   DO.obs = u(c(7,7.5,7),'mgO2 L^-1'),
#'   temp.water = u(c(25,24.5,18.9), 'degC'),
#'   pressure.air = u(c(900,903,910), 'mb'),
#'   salinity.water = u(2.43, 'PSU'))
#' @export
calc_DO_deficit <- function(DO.obs, temp.water, pressure.air, salinity.water = 0, ...){
  
  DO.equil <- calc_DO_at_sat(temp.water, pressure.air, salinity.water, ...)
  
  # to do: verify incoming units (convert if needed?) and set DO.equil units to mgO2 L^-1
  DO.deficit <- DO.equil-DO.obs
  
  return(DO.deficit)
}