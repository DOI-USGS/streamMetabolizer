#' The workhorse DO modeling function
#' 
#' Accepts GPP, ER, etc. and returns DO.mod. Used in many functions, including 
#' metab_mle() and predict_DO.metab_mle()
#' 
#' @param GPP.daily One GPP rate per day (mg O2 / L / d)
#' @param ER.daily One ER rate per day (mg O2 / L / d), always non-positive.
#' @param K600.daily One K600 per day (1 / d)
#' @param DO.sat dissolved oxygen concentrations if the water were at 
#'   equilibrium saturation \eqn{mg O[2] L^{-1}}{mg O2 / L}. Calculate using 
#'   \link{calc_DO_at_sat}
#' @param depth stream depth, \eqn{m}{m}.
#' @param temp.water stream temperature in degC
#' @param frac.GPP the fraction of daily GPP to apply to each timestep
#' @param frac.ER the fraction of daily ER to apply to each timestep
#' @param frac.D the fraction of daily D to apply to each timestep
#' @param DO.mod.1 the first DO.obs value, to which the first DO.mod value will 
#'   be set
#' @param n number of DO.mod values to produce
#' @export
#' @examples
#' calc_DO_mod(10, -13, 2.5, 14, 1, rep(12,100), 
#'   rep(1/100,100), rep(1/100,100), rep(1/100,100), 11, 100)
calc_DO_mod <- function(
  GPP.daily, ER.daily, K600.daily, DO.sat, depth, temp.water, frac.GPP, frac.ER, frac.D, DO.mod.1, n) {
  
  # partition GPP and ER into their timestep-specific rates (mg/L/timestep at
  # each timestep)
  GPP <- GPP.daily * frac.GPP / depth
  ER <- ER.daily * frac.ER / depth
  K <- convert_k600_to_kGAS(K600.daily, temperature=temp.water, gas="O2") * frac.D
  
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