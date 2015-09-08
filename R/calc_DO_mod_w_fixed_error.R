#' \code{calc_DO_mod_w_fixed_error} simulates DO with vectors of observation
#' and/or process error.
#' 
#' @rdname calc_DO_mod
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
#' @param ... additional arguments passed to other variants on calc_DO_mod
#' @param err.obs A vector of observation errors of length n. Observation errors
#'   are those applied to DO.mod after generating the full time series of 
#'   modeled values.
#' @param err.proc A vector process errors of length n. Process errors are 
#'   applied at each time step, and therefore propagate into the next timestep.
#' @export
#' @examples
#' fr <- rep(1/100,100) # shorthand for readability of next lines
#' preds <- data.frame(
#'  t=1:100,
#'  no_err = calc_DO_mod_w_fixed_error(
#'   10, -13, 2.5, 14, 1, rep(12,100), fr, fr, fr, 11, 100, 0, 0),
#'  obs_err_const = calc_DO_mod_w_fixed_error(
#'   10, -13, 2.5, 14, 1, rep(12,100), fr, fr, fr, 11, 100, 0.1, 0),
#'  obs_err_rnorm = calc_DO_mod_w_fixed_error(
#'   10, -13, 2.5, 14, 1, rep(12,100), fr, fr, fr, 11, 100, stats::rnorm(100,0,0.1), 0),
#'  proc_err_const = calc_DO_mod_w_fixed_error(
#'   10, -13, 2.5, 14, 1, rep(12,100), fr, fr, fr, 11, 100, 0, 0.01),
#'  proc_err_rnorm = calc_DO_mod_w_fixed_error(
#'   10, -13, 2.5, 14, 1, rep(12,100), fr, fr, fr, 11, 100, 0, stats::rnorm(100,0,0.05)),
#'  all_err = calc_DO_mod_w_fixed_error(
#'   10, -13, 2.5, 14, 1, rep(12,100), fr, fr, fr, 11, 100, 
#'   stats::rnorm(100,0,0.1), stats::rnorm(100,0,0.05))
#' )
#' all(preds$no_err == calc_DO_mod(10, -13, 2.5, 14, 1, rep(12,100), 
#'   rep(1/100,100), rep(1/100,100), rep(1/100,100), 11, 100))
#' \dontrun{
#' library(ggplot2); library(tidyr)
#' ggplot(gather(preds, type, value, no_err, obs_err_const, obs_err_rnorm, 
#'               proc_err_const, proc_err_rnorm, all_err), 
#'   aes(x=t, y=value, color=type)) + geom_line() + theme_bw()
#' }
calc_DO_mod_w_fixed_error <- function(
  GPP.daily, ER.daily, K600.daily, DO.sat, depth, temp.water, frac.GPP, frac.ER, frac.D, DO.mod.1, n,
  err.obs=rep(0, n), err.proc=rep(0, n), ...) {
  
  # partition GPP and ER into their timestep-specific rates (mg/L/timestep at
  # each timestep)
  GPP <- GPP.daily * frac.GPP / depth
  ER <- ER.daily * frac.ER / depth
  K <- convert_k600_to_kGAS(K600.daily, temperature=temp.water, gas="O2") * frac.D
  
  # make sure anything in the following loop has n observations
  if(length(DO.sat) != n) DO.sat <- rep(DO.sat, length.out=n) # this would be odd
  if(length(err.proc) != n) err.proc <- rep(err.proc, length.out=n) # this is more probable. err.obs can stay length 1 if it is
  
  # Model DO with given params
  DO.mod <- numeric(n)
  DO.mod[1] <- DO.mod.1
  for(i in 2:n) {
    DO.mod[i] <- DO.mod[i-1] +
      GPP[i] + 
      ER[i] + 
      K[i] * (DO.sat[i] - DO.mod[i-1]) +
      err.proc[i]
  } 
  DO.mod <- DO.mod + err.obs
  
  # Return
  DO.mod
}