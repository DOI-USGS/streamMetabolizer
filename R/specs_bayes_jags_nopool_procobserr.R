#' \code{specs_bayes_jags_nopool_procobserr} - a JAGS model with no pooling and 
#' both process and observation error. Compatible \code{model_file} options are 
#' \code{c('nopool_procobserr_pairmeans.jags', 'nopool_procobserr_Euler.jags')}.
#' 
#' @rdname specs_bayes
#'   
#' @param model_file character. The model definition file to use. The file may 
#'   be specified either as a file path relative to the streamMetabolizer 
#'   models/bayes directory (the first assumption; this directory can be found 
#'   with \code{system.file("models/bayes", package="streamMetabolizer")}) or as
#'   an absolute path or a path relative to the current working directory (the 
#'   second assumption, if the first assumption turns up no files of the given 
#'   name). For example, the default is 
#'   \code{"nopool_procobserr_pairmeans.jags"}. The containing folder is in this
#'   case \code{"models/bayes"}. The suffix, 'jags', will determine which MCMC 
#'   software package is used. The file name (in this case 
#'   \code{"nopool_procobserr_pairmeans.jags"}) will determine not only the 
#'   model file to use but also which variables are packaged and sent to the 
#'   MCMC software.
#' @param bayes_fun character in \code{c('bayes_1ply', 'bayes_all')} indicating 
#'   whether the data should be split into daily chunks first ('bayes_1ply') or 
#'   passed to the model fitting function in one big chunk ('bayes_all')
#' @param bayes_software character in \code{c('jags','stan')} indicating the 
#'   software package to use for the MCMC process
#' @param keep_mcmcs TRUE, FALSE, or (for nopool models) a vector of dates 
#'   (coerced with as.Date if character, etc.) indicating whether to keep all of
#'   the mcmc model objects (TRUE), none of them (FALSE), or specific dates. The
#'   default is FALSE because these objects can be very large.
#'   
#' @param GPP_daily_mu The mean of a dnorm distribution for GPP_daily, the daily
#'   rate of gross primary production
#' @param GPP_daily_sigma The standard deviation of a dnorm distribution for 
#'   GPP_daily, the daily rate of gross primary production
#' @param ER_daily_mu The mean of a dnorm distribution for ER_daily, the daily 
#'   rate of ecosystem respiration
#' @param ER_daily_sigma The standard deviation of a dnorm distribution for 
#'   ER_daily, the daily rate of ecosystem respiration
#' @param K600_daily_mu The mean of a dnorm distribution for K600_daily, the 
#'   daily rate of reaeration
#' @param K600_daily_sigma The standard deviation of a dnorm distribution for 
#'   K600_daily, the daily rate of reaeration
#' @param err_proc_phi_min The lower bound on a dunif distribution for 
#'   err_proc_phi, the process error autocorrelation coefficient
#' @param err_proc_phi_max The upper bound on a dunif distribution for 
#'   err_proc_phi, the process error autocorrelation coefficient
#' @param err_proc_sigma_min The lower bound on a dunif distribution for 
#'   err_proc_sigma, the standard deviation of the process error
#' @param err_proc_sigma_max The upper bound on a dunif distribution for 
#'   err_proc_sigma, the standard deviation of the process error
#' @param err_obs_sigma_min The lower bound on a dunif distribution for 
#'   err_obs_sigma, the standard deviation of the observation error
#' @param err_obs_sigma_max The upper bound on a dunif distribution for 
#'   err_obs_sigma, the standard deviation of the observation error
#'   
#' @inheritParams prepdata_bayes
#' @inheritParams mcmc_bayes
#'   
#' @export
specs_bayes_jags_nopool_procobserr <- function(
  
  # model setup (model_path will be added in metab_bayes)
  model_file = 'nopool_procobserr_pairmeans.jags',
  bayes_fun = 'bayes_1ply',
  bayes_software = 'jags',
  keep_mcmcs = FALSE,
  
  # hyperparameters
  GPP_daily_mu = 10,
  GPP_daily_sigma = 10,
  ER_daily_mu = -10,
  ER_daily_sigma = 10,
  K600_daily_mu = 10,
  K600_daily_sigma = 10,
  
  err_proc_phi_min = 0,
  err_proc_phi_max = 1,
  err_proc_sigma_min = 0,
  err_proc_sigma_max = 0.0005,
  err_obs_sigma_min = 0,
  err_obs_sigma_max = 0.5,
  
  # inheritParams prepdata_bayes
  priors = FALSE,
  
  # inheritParams mcmc_bayes
  params_out = c("GPP_daily", "ER_daily", "K600_daily", "err_obs_sigma", "err_proc_sigma", "err_proc_phi"),
  n_chains = 4, 
  n_cores = 1, 
  adapt_steps = 100, 
  burnin_steps = 40, 
  num_saved_steps = 400, 
  thin_steps = 1,
  verbose = FALSE
  
) {
  
  as.list(environment())
  
}
