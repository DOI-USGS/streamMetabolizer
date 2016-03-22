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
#' n = 24
#' DO <- calc_DO_mod(GPP.daily=3, ER.daily=-9, K600.daily=2.5, 
#'   DO.sat=11, depth=1, temp.water=rep(12,n), 
#'   frac.GPP=with(data.frame(sinpow=sin(seq(0,pi,length.out=n))^4), sinpow/sum(sinpow)), 
#'   frac.ER=rep(1/n,n), frac.D=rep(1/n,n), 
#'   DO.mod.1=8, n=n, ODE_method="Euler")
#' \dontrun{
#' plot(DO)
#' }
calc_DO_mod <- function(
  GPP.daily, ER.daily, K600.daily, DO.sat, depth, temp.water, frac.GPP, frac.ER, frac.D, DO.mod.1, n, ODE_method="pairmeans", ...) {
  
  # Model DO with given params
  calc_DO_mod_w_fixed_error(
    GPP.daily=GPP.daily, ER.daily=ER.daily, K600.daily=K600.daily, 
    DO.sat=DO.sat, depth=depth, temp.water=temp.water, 
    frac.GPP=frac.GPP, frac.ER=frac.ER, frac.D=frac.D, DO.mod.1=DO.mod.1, n=n,
    err.obs=rep(0, n), err.proc=rep(0, n), ODE_method=ODE_method, ...)
}