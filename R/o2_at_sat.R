#' Calculates the equilibrium saturation concentration of oxygen in water at the supplied conditions
#' @param temp.water a numeric vector of water temperature in degrees Celsius, 
#' or a \linkS4class{unitted} object of water temperatures.
#' @param pressure.air barometric pressure in millibars, 
#' or a \linkS4class{unitted} object of barometric pressure.
#' @param salinity.water a numeric vector of salinity in PSU, 
#' or a \linkS4class{unitted} object of salinity. Defaults to zero. 
#' @param ... additional parameters passed to \code{\link[LakeMetabolizer]{o2.at.sat.base}}
#' @importFrom LakeMetabolizer o2.at.sat.base
#' @importFrom unitted u v get_units 
#' @keywords internal
o2_at_sat <- function(temp.water, pressure.air, salinity.water = u(0,'PSU'), ...){

  with.units <- any(sapply(list(temp.water, pressure.air, salinity.water), is.unitted))
  
  if (with.units){
    # if any units are set, they all must be set and must be correct
    verify_units(temp.water, "degC")
    verify_units(pressure.air, "mb")
    verify_units(salinity.water, "PSU")
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
