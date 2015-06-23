#' The workhorse DO modeling function
#' 
#' Accepts GPP, ER, etc. and returns DO.mod. Used in many functions, including 
#' metab_mle() and predict_DO.metab_mle()
#' 
#' @param GPP.daily One GPP rate per day (mg O2 / L / d)
#' @param ER.daily One ER rate per day (mg O2 / L / d), always non-positive.
#' @param date.time date-time values in POSIXct format
#' @param DO.obs dissolved oxygen concentration observations, \eqn{mg O[2] 
#'   L^{-1}}{mg O2 / L}
#' @param DO.sat dissolved oxygen concentrations if the water were at 
#'   equilibrium saturation \eqn{mg O[2] L^{-1}}{mg O2 / L}. Calculate using 
#'   \link{calc_DO_at_sat}
#' @param depth stream depth, \eqn{m}{m}.
#' @param k.O2 gas exchange coefficients, \eqn{m d^{-1}}{m / d}.
#' @param light photosynthetically active radiation, \eqn{\mu mol\ m^{-2} 
#'   s^{-1}}{micro mols / m^2 / s}
#' @param DO.mod.1 the first DO.obs value, to which the first DO.mod value will 
#'   be set
#' @param n number of DO.mod values to produce
#' @export
#' @examples
#' calc_DO_mod(10, -13, 2.5, 14, 1, rep(12,100), 
#'   rep(1/100,100), rep(1/100,100), rep(1/100,100), 11, 100)
calc_DO_mod <- function(
  GPP.daily, ER.daily, k.600.daily, DO.sat, depth, temp.water, frac.GPP, frac.ER, frac.D, DO.mod.1, n,
  obs.err.sd, proc.err.sd, obs.err.phi, proc.err.phi) {
  
  # partition GPP and ER into their timestep-specific rates (mg/L/timestep at
  # each timestep)
  GPP <- GPP.daily * frac.GPP / depth
  ER <- ER.daily * frac.ER / depth
  K <- convert_k600_to_kGAS(k.600.daily, temperature=temp.water, gas="O2") * frac.D
  
  # make sure anything in the following loop has n observations
  if(length(DO.sat) != n) DO.sat <- rep(DO.sat, length.out=n) # this would be odd
  
  # Model DO with given params
  DO.mod <- numeric(n)
  DO.mod[1] <- DO.mod.1
  for(i in 2:n) {
    DO.mod[i] <- DO.mod[i-1] +
      GPP[i] + 
      ER[i] + 
      K[i] * (DO.sat[i] - DO.mod[i-1])
  } 
  
  # Return
  DO.mod
}