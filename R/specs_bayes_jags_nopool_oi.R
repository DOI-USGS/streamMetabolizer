#' MCMC estimation by JAGS with no pooling and only observation error
#' 
#' Compatible \code{model_file} options are
#' \code{c('nopool_oi_pairmeans.jags', 'nopool_oi_Euler.jags')}.
#' 
#' @inheritParams specs_all
#'   
#' @export
specs_bayes_jags_nopool_oi <- function(
  
  # model setup (model_path will be added in metab_bayes)
  model_file = 'nopool_oi_pairmeans.jags',
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
  
  err_obs_sigma_min = 0,
  err_obs_sigma_max = 0.5,
  
  # inheritParams prepdata_bayes
  priors = FALSE,
  
  # inheritParams mcmc_bayes
  params_out = c("GPP_daily", "ER_daily", "K600_daily", "err_obs_sigma"),
  n_chains = 4, 
  n_cores = 1, 
  adapt_steps = 250, 
  burnin_steps = 250, 
  num_saved_steps = 500, 
  thin_steps = 1,
  verbose = FALSE
  
) {
  
  as.list(environment())
  
}