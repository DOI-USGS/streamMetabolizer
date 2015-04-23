#' Calculates the equilibrium saturation concentration of oxygen in water at the supplied conditions
#' @param temp.water water temperature in degrees C
#' @param pressure.air air pressure in millibars
#' @param salinity.water a numeric vector of salinity in PSU. Defaults to zero. 
#' Length must be one or equal to length of temperature.
#' @param ... additional parameters passed to \code{\link[LakeMetabolizer]{o2.at.sat.base}}
#' @importFrom LakeMetabolizer o2.at.sat.base
#' @keywords internal
o2_at_sat <- function(temp.water, pressure.air, salinity.water = 0, ...){
  LakeMetabolizer::o2.at.sat.base(temp = temp.water, baro = pressure.air, salinity = salinity.water, ...)
}