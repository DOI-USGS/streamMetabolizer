#' #' \code{calc_DO_mod_w_sim_error} simulates DO with new, randomly generated
#' vectors of observation and/or process error each time
#' 
#' @rdname calc_DO_mod
#'   
#' @inheritParams calc_DO_mod_w_fixed_error
#' @param err.obs.sigma The sd of observation error, or 0 for no observation 
#'   error. Observation errors are those applied to DO.mod after generating the 
#'   full time series of modeled values.
#' @param err.obs.phi The autocorrelation coefficient of the observation errors,
#'   or 0 for uncorrelated errors.
#' @param err.proc.sigma The sd of process error, or 0 for no process error. 
#'   Process errors are applied at each time step, and therefore propagate into 
#'   the next timestep.
#' @param err.proc.phi The autocorrelation coefficient of the process errors, or
#'   0 for uncorrelated errors.
#' @export
#' @importFrom stats rnorm
#' @examples
#' fr <- rep(1/100,100) # shorthand for readability of next lines
#' preds <- data.frame(
#'  t=1:100,
#'  no_err = calc_DO_mod_w_sim_error(
#'   10, -13, 2.5, 14, 1, rep(12,100), fr, fr, fr, 11, 100, 0, 0, 0, 0),
#'  obs_err = calc_DO_mod_w_sim_error(
#'   10, -13, 2.5, 14, 1, rep(12,100), fr, fr, fr, 11, 100, 0.1, 0, 0, 0),
#'  ac_obs_err = calc_DO_mod_w_sim_error(
#'   10, -13, 2.5, 14, 1, rep(12,100), fr, fr, fr, 11, 100, 0.1, 0.6, 0, 0),
#'  proc_err = calc_DO_mod_w_sim_error(
#'   10, -13, 2.5, 14, 1, rep(12,100), fr, fr, fr, 11, 100, 0, 0, 0.03, 0),
#'  ac_proc_err = calc_DO_mod_w_sim_error(
#'   10, -13, 2.5, 14, 1, rep(12,100), fr, fr, fr, 11, 100, 0, 0, 0.03, 0.6),
#'  all_err = calc_DO_mod_w_sim_error(
#'   10, -13, 2.5, 14, 1, rep(12,100), fr, fr, fr, 11, 100, 0.1, 0.6, 0.03, 0.6)
#' )
#' all(preds$no_err == calc_DO_mod(10, -13, 2.5, 14, 1, rep(12,100), 
#'   rep(1/100,100), rep(1/100,100), rep(1/100,100), 11, 100))
#' \dontrun{
#' library(ggplot2); library(tidyr)
#' ggplot(gather(preds, type, value, no_err, obs_err, ac_obs_err), 
#'   aes(x=t, y=value, color=type)) + geom_line() + theme_bw()
#' ggplot(gather(preds, type, value, no_err, proc_err, ac_proc_err), 
#'   aes(x=t, y=value, color=type)) + geom_line() + theme_bw()
#' ggplot(gather(preds, type, value, no_err, all_err), 
#'   aes(x=t, y=value, color=type)) + geom_line() + theme_bw()
#' }
calc_DO_mod_w_sim_error <- function(
  GPP.daily, ER.daily, K600.daily, DO.sat, depth, temp.water, frac.GPP, frac.ER, frac.D, DO.mod.1, n,
  err.obs.sigma, err.obs.phi, err.proc.sigma, err.proc.phi, ODE_method="pairmeans", ...) {
  
  # compute errors to add to modeled data
  err.obs <- as.numeric(stats::filter(rnorm(n, 0, err.obs.sigma), filter=err.obs.phi, method="recursive"))
  err.proc <- as.numeric(stats::filter(rnorm(n, 0, err.proc.sigma), filter=err.proc.phi, method="recursive"))
    
  # Model DO with given params
  calc_DO_mod_w_fixed_error(
    GPP.daily=GPP.daily, ER.daily=ER.daily, K600.daily=K600.daily, 
    DO.sat=DO.sat, depth=depth, temp.water=temp.water, 
    frac.GPP=frac.GPP, frac.ER=frac.ER, frac.D=frac.D, DO.mod.1=DO.mod.1, n=n,
    err.obs=err.obs, err.proc=err.proc, ODE_method=ODE_method, ...)
}