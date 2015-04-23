#' Calculates the equilibrium saturation concentration of oxygen in water at the supplied conditions
#' @param temp.water water temperature in degrees C
#' @param pressure.air air pressure in millibars
#' @param salinity.water a numeric vector of salinity in PSU. Defaults to zero. 
#' Length must be one or equal to length of temperature.
#' @param ... additional parameters passed to \code{\link[LakeMetabolizer]{o2.at.sat.base}}
#' @importFrom LakeMetabolizer o2.at.sat.base
#' @importFrom unitted u v get_units 
#' @keywords internal
o2_at_sat <- function(temp.water, pressure.air, salinity.water = u(0,'PSU'), ...){
  
  units <- sapply(list(temp.water, pressure.air, salinity.water), get_units)
  
  with.units <- TRUE
  
  if (any(is.na(units)) | !all(units == c("degC", "mb", "PSU"))){
    with.units <- FALSE
  } 
  
  # units are stripped regardless
  temp.water <- v(temp.water)
  pressure.air <- v(pressure.air)
  salinity.water <- v(salinity.water)
  
  o2.at.sat <- LakeMetabolizer::o2.at.sat.base(temp = temp.water, baro = pressure.air, salinity = salinity.water, ...)
  
  if (with.units) {
    return(u(o2.at.sat, 'mg L^-1'))
  } else {
    return(o2.at.sat)
  }
    
}