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
#' Here are the essential model lines (in Stan language) that distinguish the K 
#' pooling options.
#' 
#' \tabular{ll}{
#'   \strong{\code{pool_K600}} \tab \strong{Model code}\cr
#'   
#'   \code{none} \tab \code{K600_daily ~ normal(K600_daily_mu, K600_daily_sigma)} \cr
#'   
#'   \code{normal} \tab \code{K600_daily ~ normal(K600_daily_mu, K600_daily_sigma)}\cr
#'   \tab \code{K600_daily_mu ~ normal(K600_daily_mu_mu, K600_daily_mu_sigma)}\cr
#'   \tab \code{K600_daily_sigma ~ gamma(K600_daily_sigma_shape, K600_daily_sigma_rate)}\cr
#'   
#'   \code{linear} \tab \code{K600_daily_pred <- K600_daily_beta[1] + K600_daily_beta[2] * discharge_daily}\cr 
#'   \tab \code{K600_daily ~ normal(K600_daily_pred, K600_daily_sigma)}\cr 
#'   \tab \code{K600_daily_beta ~ normal(K600_daily_beta_mu, K600_daily_beta_sigma)}\cr
#'   \tab \code{K600_daily_sigma ~ gamma(K600_daily_sigma_shape, K600_daily_sigma_rate)}\cr
#'   
#'   \code{binned} \tab \code{K600_daily_pred <- K600_daily_beta[Q_bin_daily]}\cr
#'   \tab \code{K600_daily ~ normal(K600_daily_pred, K600_daily_sigma)}\cr
#'   \tab \code{K600_daily_beta ~ normal(K600_daily_beta_mu, K600_daily_beta_sigma)}\cr
#'   \tab \code{K600_daily_sigma ~ gamma(K600_daily_sigma_shape, K600_daily_sigma_rate)}\cr
#'   
#'   \code{complete} \tab [This option refers to complete pooling via \code{metab_Kmodel} in conjunction with preceding\cr
#'   \tab estimates of K (e.g., by \code{metab_mle} or \code{metab_night}) and
#'   subsequent estimates of GPP and ER\cr
#'   \tab (e.g., by \code{metab_mle} with daily K600 values specified)]\cr
#'  }
#' }
#' 
#' @seealso The converse of this function is \code{\link{mm_parse_name}}.
#'   
#' @param type character. The model type. Options:
#'   \itemize{
#'     \item \code{mle}: maximum likelihood estimation (see also \code{\link{metab_mle}}) 
#'     \item \code{bayes}: bayesian hierarchical models \code{\link{metab_bayes}}
#'     \item \code{night}: nighttime regression (see also \code{\link{metab_night}}) 
#'     \item \code{Kmodel}: regression of \emph{daily} estimates of
#'     \code{K600.daily} versus discharge, time, etc., usually for 3-phase 
#'     estimation of K alone (by MLE or nighttime regression), K vs discharge 
#'     (using this model), and then GPP and ER with fixed K (by MLE) (see also 
#'     \code{\link{metab_Kmodel}})
#'     \item \code{sim}: simulation of \code{DO.obs} 'data' for testing other models (see also \code{\link{metab_sim}})
#'   }
#' @param pool_K600 character. [How] should the model pool information among 
#'   days to get more consistent daily estimates for K600? Options (see Details
#'   for more):
#'   \itemize{
#'     \item \code{none}: no pooling of K600
#'     \item \code{normal}: \eqn{K600 ~ N(mu, sigma)}
#'     \item \code{linear}: \eqn{K600 ~ N(B[0] + B[1]*Q, sigma)}
#'     \item \code{binned}: \eqn{K600 ~ N(B[Q_bin], sigma)} where \eqn{mu ~
#'     N(mu_mu, mu_sigma)} and \eqn{sigma ~ N(sigma_mu, sigma_sigma)}
#'     \item \code{complete}: applicable only for \code{type='Kmodel'}, which is
#'     generally used in conjunction with preceding estimates of K (e.g., by
#'     \code{type='mle'} or \code{type='night'}) and subsequent estimates of GPP
#'     and ER (e.g., by \code{type='mle'} with daily K600 values specified)
#'   }
#' @param err_obs_iid logical. Should IID observation error be included? If not,
#'   the model will be fit to the differences in successive DO measurements, 
#'   rather than to the DO measurements themselves.
#' @param err_proc_acor logical. Should autocorrelated process error (with the 
#'   autocorrelation term phi fitted) be included?
#' @param err_proc_iid logical. Should IID process error be included?
#' @param ode_method character. The method to use in solving the ordinary
#'   differential equation for DO. Options:
#'   \itemize{
#'     \item \code{euler}, formerly \code{Euler}: the final change in DO from
#'     t=1 to t=2 is solely a function of GPP, ER, DO, etc. at t=1
#'     \item \code{trapezoid}, formerly \code{pairmeans}: the final change in DO
#'     from t=1 to t=2 is a function of the mean values of GPP, ER, etc. across
#'     t=1 and t=2.
#'     \item for \code{type='mle'}, options also include \code{rk2} and any 
#'     character method accepted by \code{\link[deSolve]{ode}} in the 
#'     \code{deSolve} package (\code{lsoda}, \code{lsode}, \code{lsodes},
#'     \code{lsodar}, \code{vode}, \code{daspk}, \code{rk4}, \code{ode23},
#'     \code{ode45}, \code{radau}, \code{bdf}, \code{bdf_d}, \code{adams},
#'     \code{impAdams}, and \code{impAdams_d}; note that many of these have not
#'     been well tested in the context of \code{streamMetabolizer} models)
#'   }
#' @param GPP_fun character. Function dictating how gross primary productivity
#'   (GPP) varies within each day. Options:
#'   \itemize{ 
#'     \item \code{linlight}: GPP is a linear function of light with an
#'     intercept at 0 and a slope that varies by day.
#'       \cr \code{GPP(t) = GPP.daily * light(t) / mean.light}
#'       \itemize{
#'         \item \code{GPP.daily}: the daily mean GPP, which is partitioned into
#'         timestep-specific rates according to the fraction of that day's
#'         average light that occurs at each timestep (specifically,
#'         \code{mean.light} is the mean of the first 24 hours of the date's
#'         data window)
#'       } 
#'     \item \code{satlight}: GPP is a saturating function of light.
#'       \cr \code{GPP(t) = Pmax * tanh(alpha * light(t) / Pmax)}
#'       \itemize{
#'         \item \code{Pmax}: the maximum possible GPP
#'         \item \code{alpha}: a descriptor of the rate of increase of GPP as a function of light
#'       }
#'     \item \code{satlightq10temp}: GPP is a saturating function of light and
#'     an exponential function of temperature.
#'       \cr \code{GPP(t) = Pmax * tanh(alpha * light(t) / Pmax) * 1.036 ^ (temp.water(t) - 20)}
#'       \itemize{
#'         \item \code{Pmax}: the maximum possible GPP
#'         \item \code{alpha}: a descriptor of the rate of increase of GPP as a function of light
#'       }
#'     \item \code{NA}: applicable only to \code{type='Kmodel'}, for which GPP is not estimated
#'   }
#' @param ER_fun character. Function dictating how ecosystem respiration (ER)
#'   varies within each day. Options:
#'   \itemize{
#'     \item \code{constant}: ER is constant over every timestep of the day.
#'       \cr \code{ER(t) = ER.daily} 
#'       \itemize{
#'         \item \code{ER.daily}: the daily mean ER, which is equal to instantaneous ER at all times
#'       }
#'     \item \code{q10temp}: ER at each timestep is an exponential function of
#'     the water temperature and a temperature-normalized base rate.
#'       \cr \code{ER(t) = ER20 * 1.045 ^ (temp.water(t) - 20)}
#'       \itemize{
#'         \item \code{ER20}: the value of ER when \code{temp.water} is 20 degrees C
#'       }
#'     \item \code{NA}: applicable only to \code{type='Kmodel'}, for which ER is not estimated
#'   }
#' @param deficit_src character. From what DO estimate (observed or modeled)
#'   should the DO deficit be computed? Options:
#'   \itemize{
#'     \item \code{DO_mod}: the DO deficit at time t will be {(DO.sat(t) -
#'     DO_mod(t))}, the difference between the equilibrium-saturation value and
#'     the current best estimate of the true DO concentration at that time
#'     \item \code{DO_obs}: the DO deficit at time t will be {(DO.sat(t) -
#'     DO.obs(t))}, the difference between the equilibrium-saturation value and
#'     the measured DO concentration at that time
#'     \item \code{DO_obs_filter}: applicable only to \code{type='night'}: a
#'     smoothing filter is applied over the measured DO.obs values before
#'     applying nighttime regression
#'     \item \code{NA}: applicable only to \code{type='Kmodel'}, for which DO deficit is not estimated
#'   }
#' @param engine character. With which function or software should the model
#'   fitting be done?
#'   \itemize{
#'     \item for \code{type='mle'}: \code{nlm} only (the default)
#'     \item for \code{type='bayes'}: \code{stan} (the default) or \code{jags}, 
#'     both of which are external software packages that run MCMC chains for 
#'     Bayesian models. See http://mc-stan.org and
#'     http://mcmc-jags.sourceforge.net, respectively
#'     \item for \code{type='night'}: \code{lm} only (the default)
#'     \item for \code{type='Kmodel'}: \code{mean}, \code{lm}, or \code{loess}
#'     enable different types of relationships between daily K600 and its
#'     predictors (nothing, discharge, time, etc.)
#'     \item for \code{type='sim'}: \code{rnorm} only (the default)
#'   }
#' @param check_validity logical. if TRUE, this function checks the resulting
#'   name against \code{mm_valid_names(type)}.
#' @import dplyr
#' @export
#' @examples 
#' mm_name('mle')
#' mm_name('mle', GPP_fun='satlight', ER_fun='q10temp')
#' mm_name('night')
#' mm_name('sim', err_proc_acor=TRUE)
#' mm_name('bayes', pool_K600='binned')
mm_name <- function(
  type=c('mle','bayes','night','Kmodel','sim'), 
  #pool_GPP='none', pool_ER='none', pool_eoi='alldays', pool_epc='alldays', pool_epi='alldays',
  pool_K600=c('none','normal','linear','binned','complete'),
  err_obs_iid=c(TRUE, FALSE),
  err_proc_acor=c(FALSE, TRUE),
  err_proc_iid=c(FALSE, TRUE),
  ode_method=c('trapezoid','euler','rk2','lsoda','lsode','lsodes','lsodar','vode','daspk',
               'rk4','ode23','ode45','radau','bdf','bdf_d','adams','impAdams','impAdams_d',
               'Euler','pairmeans','NA'),
  GPP_fun=c('linlight','satlight','satlightq10temp','NA'),
  ER_fun=c('constant','q10temp','NA'),
  deficit_src=c('DO_mod','DO_obs','DO_obs_filter','NA'),
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
    pool_K600='complete'
    pool_all='complete'
    err_obs_iid=FALSE
    err_proc_acor=FALSE
    err_proc_iid=FALSE
    ode_method='NA'
    GPP_fun='NA'
    ER_fun='NA'
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
  
  # check arguments and throw errors as needed. these checks define the names 
  # that are possible to create; will be supplemented by call to mm_valid_names 
  # to see if a specific arg combo is actually implemented
  if(type != 'Kmodel') {
    pool_K600 <- match.arg(pool_K600)
    pool_all <- if(pool_K600 == 'none') 'none' else 'partial'
    if(!is.logical(err_obs_iid) || length(err_obs_iid) != 1) stop("need err_obs_iid to be a logical of length 1")
    if(!is.logical(err_proc_acor) || length(err_proc_acor) != 1) stop("need err_proc_acor to be a logical of length 1")
    if(!is.logical(err_proc_iid) || length(err_proc_iid) != 1) stop("need err_proc_iid to be a logical of length 1")
    ode_method <- match.arg(ode_method)
    if(ode_method %in% c('Euler','pairmeans'))
      warning("for ode_method, 'Euler' and 'pairmeans' are deprecated in favor of 'euler' and 'trapezoid'")
    GPP_fun <- match.arg(GPP_fun)
    ER_fun <- match.arg(ER_fun)
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
    c(none='', normal='Kn', linear='Kl', binned='Kb', complete='Kc')[[pool_K600]],
    c(none='np', partial='', complete='')[[pool_all]], '_',
    if(err_obs_iid) 'oi', if(err_proc_acor) 'pc', if(err_proc_iid) 'pi', '_',
    c(Euler='Eu', pairmeans='pm', trapezoid='tr', rk2='r2', 
      lsoda='o1', lsode='o2', lsodes='o3', lsodar='o4', vode='o5', daspk='o6', euler='eu', rk4='o8', 
      ode23='o9', ode45='o10', radau='o11', bdf='o12', bdf_d='o13', adams='o14', impAdams='o15', impAdams_d='o16',
      'NA'='')[[ode_method]], '_',
    c(linlight='pl', satlight='ps', satlightq10temp='pq', 'NA'='')[[GPP_fun]],
    c(constant='rc', q10temp='rq', 'NA'='')[[ER_fun]],
    c(DO_mod='km', DO_obs='ko', DO_obs_filter='kf', 'NA'='')[[deficit_src]], 
    '.', engine)
  
  # check validity if requested
  check_validity <- if(!is.logical(check_validity)) stop("need check_validity to be a logical of length 1") else check_validity[1]
  if(isTRUE(check_validity)) mm_validate_name(mmname)
  
  # return
  mmname
}