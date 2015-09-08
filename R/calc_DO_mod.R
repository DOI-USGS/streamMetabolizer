#' Calculate a time series of DO concentrations from GPP, ER, K600, and other
#' inputs
#' 
#' \code{calc_DO_mod} simulates DO with no observation or process error.
#' 
#' Accepts GPP, ER, etc. and returns DO.mod. Used in many functions, including 
#' metab_mle() and predict_DO.metab_mle()
#' 
#' @inheritParams calc_DO_mod_w_fixed_error
#' @export
#' @examples
#' calc_DO_mod(10, -13, 2.5, 14, 1, rep(12,100), 
#'   rep(1/100,100), rep(1/100,100), rep(1/100,100), 11, 100)
calc_DO_mod <- function(
  GPP.daily, ER.daily, K600.daily, DO.sat, depth, temp.water, frac.GPP, frac.ER, frac.D, DO.mod.1, n, ...) {
  
  # Model DO with given params
  calc_DO_mod_w_fixed_error(
    GPP.daily=GPP.daily, ER.daily=ER.daily, K600.daily=K600.daily, 
    DO.sat=DO.sat, depth=depth, temp.water=temp.water, 
    frac.GPP=frac.GPP, frac.ER=frac.ER, frac.D=frac.D, DO.mod.1=DO.mod.1, n=n,
    err.obs=rep(0, n), err.proc=rep(0, n))
}