#' Find the name of a model by its features
#' 
#' A \code{model_name} concisely specifies the structure of a metabolism model. 
#' From a \code{model_name}, an appropriate set of model specifications 
#' (parameters and runtime options) can be generated with \code{\link{specs}}. 
#' From a complete \code{specs} list, a metabolism model can be run with 
#' \code{\link{metab}}.
#' 
#' While the \code{Usage} shows all valid values for each argument, not all 
#' argument combinations are valid; the combination will also be checked if 
#' \code{check_validity==TRUE}. For arguments not explicitly specified, defaults
#' depend on the value of \code{type}: any argument that is not explicitly
#' supplied (besides \code{type} and \code{check_validity}) will default to the
#' values indicated by \code{mm_parse_name(mm_valid_names(type)[1])}.
#' 
#' @details
#' 
#' \subsection{pool_K600}{
#' 
#' Here are the essential model lines (in Stan language) that separate the four 
#' K pooling options.
#' 
#' \tabular{l}{ \strong{\code{pool_K600 = 'none'}}\cr \code{K600_daily ~ 
#' normal(K600_daily_mu, K600_daily_sigma)} }
#' 
#' \tabular{l}{ \strong{\code{pool_K600 = 'normal'}}\cr \code{K600_daily ~ 
#' normal(K600_daily_mu, K600_daily_sigma)}\cr \code{K600_daily_mu ~ 
#' normal(K600_daily_mu_mu, K600_daily_mu_sigma)}\cr \code{K600_daily_sigma ~ 
#' gamma(K600_daily_sigma_shape, K600_daily_sigma_rate)} }
#' 
#' \tabular{l}{ \strong{\code{pool_K600 = 'linear'}}\cr \code{K600_daily_pred <-
#' K600_daily_beta[1] + K600_daily_beta[2] * discharge_daily}\cr 
#' \code{K600_daily ~ normal(K600_daily_pred, K600_daily_sigma)}\cr 
#' \code{K600_daily_beta ~ normal(K600_daily_beta_mu, K600_daily_beta_sigma)}\cr
#' \code{K600_daily_sigma ~ gamma(K600_daily_sigma_shape, 
#' K600_daily_sigma_rate)} }
#' 
#' \tabular{l}{ \strong{\code{pool_K600 = 'binned'}}\cr \code{K600_daily_pred <-
#' K600_daily_beta[Q_bin_daily]}\cr \code{K600_daily ~ normal(K600_daily_pred, 
#' K600_daily_sigma)}\cr \code{K600_daily_beta ~ normal(K600_daily_beta_mu, 
#' K600_daily_beta_sigma)}\cr \code{K600_daily_sigma ~ 
#' gamma(K600_daily_sigma_shape, K600_daily_sigma_rate)} }
#' 
#' }
#' 
#' @seealso The converse of this function is \code{\link{mm_parse_name}}.
#'   
#' @param type the model type, corresponding to the model fitting function 
#'   (\code{\link{metab_bayes}}, \code{\link{metab_mle}}, etc.)
#' @param pool_K600 Should the model pool information among days to get more 
#'   consistent daily estimates for K600? Options: \code{'none'}=no pooling of 
#'   K600; \code{'normal'}=\eqn{K600 ~ N(mu, sigma)}; \code{'linear'}=\eqn{K600 
#'   ~ N(B[0] + B[1]*Q, sigma)}; \code{'binned'}=\eqn{K600 ~ N(B[Q_bin], sigma)}
#'   where \eqn{mu ~ N(mu_mu, mu_sigma)} and \eqn{sigma ~ N(sigma_mu, 
#'   sigma_sigma)}. See Details for more.
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
#' @param engine Which software are we generating code for?
#' @param check_validity if TRUE, checks the resulting name against 
#'   mm_valid_names(type).
#' @import dplyr
#' @export
#' @examples 
#' mm_name('mle')
#' mm_name('night')
#' mm_name('sim')
#' mm_name('bayes')
mm_name <- function(
  type=c('mle','bayes','night','Kmodel','sim'), 
  #pool_GPP='none', pool_ER='none', pool_eoi='alldays', pool_epc='alldays', pool_epi='alldays',
  pool_K600=c('none','normal','linear','binned'),
  err_obs_iid=c(TRUE, FALSE),
  err_proc_acor=c(FALSE, TRUE),
  err_proc_iid=c(FALSE, TRUE),
  ode_method=c('pairmeans','Euler','NA'),
  deficit_src=c('DO_mod','DO_obs','NA'),
  engine=c('stan','jags','nlm','lm','mean','loess','rnorm'),
  check_validity=TRUE) {
  
  # determine type
  type <- match.arg(type)
  
  # set type-specific defaults where values weren't specified
  . <- '.dplyr.var'
  if(type != 'Kmodel') {
    relevant_args <- names(formals(mm_name)) %>% .[!(. %in% c('type','check_validity'))]
  } else {
    # only one argument allowed for Kmodel
    relevant_args <- 'engine' 
    # directly specify all the rest
    pool_K600='none'
    pool_all='none'
    err_obs_iid=FALSE
    err_proc_acor=FALSE
    err_proc_iid=FALSE
    ode_method='NA'
    deficit_src='NA'
  }
  given_args <- names(match.call()[-1])
  missing_args <- relevant_args[!(relevant_args %in% given_args)]
  if(length(missing_args) > 0) {
    default_args <- mm_parse_name(mm_valid_names(type)[1])
    for(ms in missing_args) {
      assign(ms, default_args[[ms]])
    }
  }
  
  # check arguments and throw errors as needed
  if(type != 'Kmodel') {
    pool_K600 <- match.arg(pool_K600)
    pool_all <- if(pool_K600 == 'none') 'none' else 'partial'
    if(!is.logical(err_obs_iid) || length(err_obs_iid) != 1) stop("need err_obs_iid to be a logical of length 1")
    if(!is.logical(err_proc_acor) || length(err_proc_acor) != 1) stop("need err_proc_acor to be a logical of length 1")
    if(!is.logical(err_proc_iid) || length(err_proc_iid) != 1) stop("need err_proc_iid to be a logical of length 1")
    ode_method <- match.arg(ode_method)
    deficit_src <- match.arg(deficit_src)
    if(type=='bayes' && !err_obs_iid && deficit_src == 'DO_mod') stop("for bayesian models, if there's no err_obs, deficit_src must be DO_obs")
  } else {
    if(any(!(given_args %in% c('type','engine','check_validity'))))
       stop("for Kmodel, only type, engine, and check_validity may be specified")
  }
  engine <- match.arg(engine)
  if(!(engine %in% list(bayes=c('stan','jags'), mle='nlm', night='lm', Kmodel=c('mean','lm','loess'), sim='rnorm')[[type]]))
    stop("mismatch between type (",type,") and engine (",engine,")")
  
  # make the name
  mmname <- paste0(
    c(bayes='b', mle='m', night='n', Kmodel='K', sim='s')[[type]], '_',
    c(none='', normal='Kn', linear='Kl', binned='Kb')[[pool_K600]],
    c(none='np', partial='')[[pool_all]], '_',
    if(err_obs_iid) 'oi', if(err_proc_acor) 'pc', if(err_proc_iid) 'pi', '_',
    c(Euler='eu', pairmeans='pm', 'NA'='')[[ode_method]], '_',
    c(DO_mod='km', DO_obs='ko', 'NA'='')[[deficit_src]], '.',
    engine)
  
  # check validity if requested
  check_validity <- if(!is.logical(check_validity)) stop("need check_validity to be a logical of length 1") else check_validity[1]
  if(isTRUE(check_validity)) mm_validate_name(mmname)
  
  # return
  mmname
}