#' Estimate depth from discharge and hydraulic geometry coefficients
#' 
#' Uses the relationship \eqn{d=c*Q^f} (parameter names and definitions as in 
#' Leopold and Maddock, 1953; default values for c and f as in Raymond et al. 
#' 2012)
#' 
#' @param Q discharge (m^3 s^-1)
#' @param c coefficient representing depth at unit discharge (usually m)
#' @param f exponent in depth-discharge relation (unitless)
#' @return d, stream depth, in the same units as c
#' @examples
#' Qs <- seq(1,9,2)
#' calc_depth(Q=Qs)
#' calc_depth(Q=Qs, f=0.4)
#' library(unitted)
#' calc_depth(Q=u(Qs, "m^3 s^-1"), c=u(40,"cm"))
#' calc_depth(Q=u(Qs, "m^3 s^-1"), f=u(0.36))
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
calc_depth <- function(Q, c=u(0.409,"m"), f=u(0.294,"")) {

  with.units <- is.unitted(Q) || (if(!missing(c)) is.unitted(c) else FALSE) || (if(!missing(f)) is.unitted(f) else FALSE)

  if(with.units) {
    # if any units are set, they all must be set and must be correct
    verify_units(Q, "m^3 s^-1")
    if(!(get_units(c) %in% c("m","cm","mm","ft","in"))) warning("c has unknown depth units (",get_units(c),")")
    verify_units(f, "")
  } else {
    # if no units are explicitly set, then make sure c and f aren't using their unitted defaults
    c <- v(c)
    f <- v(f)
  }

  # The exponential form (below) is equivalent to this log-log form:
  #   d <- exp(log(v(c)) + log(v(Q)) * f)
  #   if(with.units) d <- u(d, get_units(c))
    
  # Do the calculation, overriding the Q units (if any) for this empirical
  # equation. If c is unitted, units will be carried through.
  c * v(Q) ^ f
}