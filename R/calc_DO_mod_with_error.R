#' Simulates DO with error
#' 
#' Accepts GPP, ER, etc. and returns DO.mod with error added (observation and/or
#' process error, with or without autocorrelation)
#' 
#' @inheritParams calc_DO_mod
#' @param obs.err.sd The sd of observation error, or 0 for no observation error.
#'   Observation errors are those applied to DO.mod after generating the full 
#'   time series of modeled values.
#' @param obs.err.phi The autocorrelation coefficient of the observation errors,
#'   or 0 for uncorrelated errors.
#' @param proc.err.sd The sd of process error, or 0 for no process error. 
#'   Process errors are applied at each time step, and therefore propagate into 
#'   the next timestep.
#' @param proc.err.phi The autocorrelation coefficient of the process errors, or
#'   0 for uncorrelated errors.
#' @export
#' @examples
#' fr <- rep(1/100,100) # shorthand for readability of next lines
#' preds <- data.frame(
#'  t=1:100,
#'  no_err = calc_DO_mod_with_error(
#'   10, -13, 2.5, 14, 1, rep(12,100), fr, fr, fr, 11, 100, 0, 0, 0, 0),
#'  obs_err = calc_DO_mod_with_error(
#'   10, -13, 2.5, 14, 1, rep(12,100), fr, fr, fr, 11, 100, 0.1, 0, 0, 0),
#'  ac_obs_err = calc_DO_mod_with_error(
#'   10, -13, 2.5, 14, 1, rep(12,100), fr, fr, fr, 11, 100, 0.1, 0.6, 0, 0),
#'  proc_err = calc_DO_mod_with_error(
#'   10, -13, 2.5, 14, 1, rep(12,100), fr, fr, fr, 11, 100, 0, 0, 0.03, 0),
#'  ac_proc_err = calc_DO_mod_with_error(
#'   10, -13, 2.5, 14, 1, rep(12,100), fr, fr, fr, 11, 100, 0, 0, 0.03, 0.6),
#'  all_err = calc_DO_mod_with_error(
#'   10, -13, 2.5, 14, 1, rep(12,100), fr, fr, fr, 11, 100, 0.1, 0.6, 0.03, 0.6)
#' )
#' all(preds$no_err == calc_DO_mod(10, -13, 2.5, 14, 1, rep(12,100), 
#'   rep(1/100,100), rep(1/100,100), rep(1/100,100), 11, 100))
#' library(ggplot2); library(tidyr)
#' ggplot(gather(preds, type, value, no_err, obs_err, ac_obs_err), 
#'   aes(x=t, y=value, color=type)) + geom_line()
#' ggplot(gather(preds, type, value, no_err, proc_err, ac_proc_err), 
#'   aes(x=t, y=value, color=type)) + geom_line()
#' ggplot(gather(preds, type, value, no_err, all_err), 
#'   aes(x=t, y=value, color=type)) + geom_line()
calc_DO_mod_with_error <- function(
  GPP.daily, ER.daily, K600.daily, DO.sat, depth, temp.water, frac.GPP, frac.ER, frac.D, DO.mod.1, n,
  obs.err.sd, obs.err.phi, proc.err.sd, proc.err.phi) {
  
  # partition GPP and ER into their timestep-specific rates (mg/L/timestep at
  # each timestep)
  GPP <- GPP.daily * frac.GPP / depth
  ER <- ER.daily * frac.ER / depth
  K <- convert_k600_to_kGAS(K600.daily, temperature=temp.water, gas="O2") * frac.D
  
  # make sure anything in the following loop has n observations
  if(length(DO.sat) != n) DO.sat <- rep(DO.sat, length.out=n) # this would be odd
  
  # compute errors to add to modeled data
  obs.err <- as.numeric(stats::filter(rnorm(n, 0, obs.err.sd), filter=obs.err.phi, method="recursive"))
  proc.err <- as.numeric(stats::filter(rnorm(n, 0, proc.err.sd), filter=proc.err.phi, method="recursive"))
    
  # Model DO with given params
  DO.mod <- numeric(n)
  DO.mod[1] <- DO.mod.1
  for(i in 2:n) {
    DO.mod[i] <- DO.mod[i-1] +
      GPP[i] + 
      ER[i] + 
      K[i] * (DO.sat[i] - DO.mod[i-1]) +
      proc.err[i]
  } 
  DO.mod <- DO.mod + obs.err
  
  # Return
  DO.mod
}