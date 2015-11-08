#' Name the Bayesian model with the desired features
#' 
#' @seealso The converse of this function is \code{\link{mm_parse_name}}.
#'   
#' @param type the model type, corresponding to the model fitting function 
#'   (\code{\link{metab_bayes}}, \code{\link{metab_mle}}, etc.)
#' @param pooling Should the model pool information among days to get more 
#'   consistent daily estimates?
#' @param err_obs_iid logical. Should IID observation error be included? If not,
#'   the model will be fit to the differences in successive DO measurements, 
#'   rather than to the DO measurements themselves.
#' @param err_proc_acor logical. Should autocorrelated process error (with the 
#'   autocorrelation term phi fitted) be included?
#' @param err_proc_iid logical. Should IID process error be included?
#' @param ode_method The method to use in solving the ordinary differential 
#'   equation for DO. Euler: dDOdt from t=1 to t=2 is solely a function of GPP, 
#'   ER, DO, etc. at t=1. pairmeans: dDOdt from t=1 to t=2 is a function of the 
#'   mean values of GPP, ER, etc. across t=1 and t=2.
#' @param deficit_src From what DO estimate (observed or modeled) should the DO 
#'   deficit be computed?
#' @param bayes_software Which software are we generating code for?
#' @export
mm_name <- function(
  type=c('bayes','mle','night','sim'), 
  pooling='none',
  err_obs_iid=c(TRUE, FALSE),
  err_proc_acor=c(FALSE, TRUE),
  err_proc_iid=c(FALSE, TRUE),
  ode_method=c('pairmeans','Euler','NA'),
  deficit_src=c('DO_mod','DO_obs','NA'),
  bayes_software=c('stan','jags','nlm','lm','rnorm')) {
  
  # choose/check arguments
  type <- match.arg(type)
  pooling <- match.arg(pooling)
  err_obs_iid <- if(!is.logical(err_obs_iid)) stop("need err_obs_iid to be a logical of length 1") else err_obs_iid[1]
  err_proc_acor <- if(!is.logical(err_proc_acor)) stop("need err_proc_acor to be a logical of length 1") else err_proc_acor[1]
  err_proc_iid <- if(!is.logical(err_proc_iid)) stop("need err_proc_iid to be a logical of length 1") else err_proc_iid[1]
  ode_method <- match.arg(ode_method)
  deficit_src <- match.arg(deficit_src)
  bayes_software <- match.arg(bayes_software)
  
  # make the name
  paste0(
    c(bayes='b', mle='m', night='n', sim='s')[[type]], '_',
    c(none='np', partial='pp')[[pooling]], '_',
    if(err_obs_iid) 'oi', if(err_proc_acor) 'pc', if(err_proc_iid) 'pi', '_',
    c(Euler='eu', pairmeans='pm', 'NA'='')[[ode_method]], '_',
    c(DO_mod='km', DO_obs='ko', 'NA'='')[[deficit_src]], '.',
    bayes_software)
}