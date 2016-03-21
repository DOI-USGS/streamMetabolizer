#' Calculates the equilibrium saturation concentration of oxygen in water at the
#' supplied conditions
#' 
#' @param temp.water a numeric vector of water temperature in degrees Celsius, 
#'   or a \linkS4class{unitted} object of water temperatures.
#' @param pressure.air barometric pressure in millibars, or a 
#'   \linkS4class{unitted} object of barometric pressure.
#' @param salinity.water a numeric vector of salinity in PSU, or a 
#'   \linkS4class{unitted} object of salinity. Defaults to zero.
#' @param model character. One of 'garcia-benson', 'garcia', 'weiss', or
#'   'benson', but 'garcia-benson' is recommended.
#' @param ... additional parameters passed to 
#'   \code{\link[LakeMetabolizer]{o2.at.sat.base}}
#' @return a numeric vector of dissolved oxygen equilibrium saturation 
#'   concentrations, in mg/L, with units attached if any of the input vectors 
#'   are unitted.
#'   
#' @importFrom LakeMetabolizer o2.at.sat.base
#' @importFrom unitted u v get_units verify_units is.unitted
#' @examples
#' calc_DO_at_sat(temp=21, press=1000.1, sal=0) # no units checking if no units provided
#' library(unitted)
#' calc_DO_at_sat(temp=u(21,"degC"), press=u(1000.1,"mb"), sal=u(0,"PSU")) # units are checked
#' @export
calc_DO_at_sat <- function(temp.water, pressure.air, salinity.water = u(0,'PSU'), model='garcia-benson', ...){

  with.units <- any(sapply(list(temp.water, pressure.air), is.unitted)) || (if(!missing(salinity.water)) is.unitted(salinity.water) else FALSE)
  
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
  
  o2.at.sat <- LakeMetabolizer::o2.at.sat.base(temp = temp.water, baro = pressure.air, salinity = salinity.water, model = model, ...)
  
  if (with.units) {
    return(u(o2.at.sat, 'mgO2 L^-1'))
  } else {
    return(o2.at.sat)
  }
    
}
