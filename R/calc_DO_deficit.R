#'@title calculate a vector of dissolved oxygen deficit
#'@description Creates a DO.deficit vector for input into various metabolism models (see...)
#'@param DO.obs
#'@param temp.water a numeric vector of water temperature in degrees Celsius.
#'@param pressure.air barometric pressure in millibars
#'@param salinity.water a numeric vector of salinity in PSU. Defaults to zero. 
#' Length must be one or equal to length of temperature.
#'@param ... additional parameters passed to \code{\link[LakeMetabolizer]{o2.at.sat.base}}
#'@return a vector of DO.deficit values 
#'@examples
#'DO.obs = 7
#'temp.water =25
#'pressure.air = 900
#'salinity.water = 2.43
#'DO.deficit <- calc_DO_deficit(DO, temp, prss.mb, salinity.water)
#'@export
calc_DO_deficit <- function(DO.obs, temp.water, pressure.air, salinity.water = 0, ...){
  DO.equil <- o2_at_sat(temp.water, pressure.air, salinity.water, ...)
  DO.deficit <- DO.equil-DO.obs
  return(DO.deficit)
}