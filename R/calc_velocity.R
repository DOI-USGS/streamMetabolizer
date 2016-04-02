#' Estimate velocity from discharge and hydraulic geometry coefficients
#' 
#' Uses the relationship \eqn{U=k*Q^m} (parameter names and definitions as in 
#' Leopold and Maddock, 1953; default values for k and m as in Raymond et al. 
#' 2012)
#' 
#' @param Q discharge (m^3 s^-1)
#' @param k coefficient representing velocity at unit discharge (usually m/s; e in Raymond et al.)
#' @param m exponent in velocity-discharge relation (unitless; f in Raymond et al.)
#' @return v (= V = U), stream flow velcoity, in the same units as k
#' @examples
#' Qs <- seq(1,9,2)
#' calc_velocity(Q=Qs)
#' calc_velocity(Q=Qs, k=0.4)
#' library(unitted)
#' calc_velocity(Q=u(Qs, "m^3 s^-1"), m=u(40))
#' calc_velocity(Q=u(Qs, "m^3 s^-1"), k=u(0.36, "m s^-1"))
#' @references Raymond, Peter A., Christopher J. Zappa, David Butman, Thomas L. 
#'   Bott, Jody Potter, Patrick Mulholland, Andrew E. Laursen, William H. 
#'   McDowell, and Denis Newbold. \emph{Scaling the gas transfer velocity and 
#'   hydraulic geometry in streams and small rivers}. Limnology & Oceanography: 
#'   Fluids & Environments 2 (2012): 41-53.
#'   
#'   Leopold, L.B., and Thomas Maddock Jr. \emph{The Hydraulic Geometry of
#'   Stream Channels and Some Physiographic Implications}. Report. Professional
#'   Paper, 1953. USGS Publications Warehouse.
#'   http://pubs.er.usgs.gov/publication/pp252.
#'   
#' @importFrom unitted u v verify_units
#' @export
calc_velocity <- function(Q, k=u(0.194,"m s^-1"), m=u(0.285,"")) {
  
  with.units <- is.unitted(Q) || (if(!missing(k)) is.unitted(k) else FALSE) || (if(!missing(m)) is.unitted(m) else FALSE)
  
  if(with.units) {
    # if any units are set, they all must be set and must be correct
    verify_units(Q, "m^3 s^-1")
    if(!(get_units(k) %in% paste0(c("m","cm","mm","ft","in"), " s^-1"))) warning("c has unknown depth units (",get_units(k),")")
    verify_units(m, "")
  } else {
    # if no units are explicitly set, then make sure c and f aren't using their unitted defaults
    k <- v(k)
    m <- v(m)
  }
  
  # The exponential form (below) is equivalent to this log-log form:
  #   U <- exp(log(v(k)) + log(v(Q)) * m)
  #   if(with.units) U <- u(U, get_units(k))
  
  # Do the calculation, overriding the Q units (if any) for this empirical
  # equation. If k is unitted, units will be carried through.
  k * v(Q) ^ m
}