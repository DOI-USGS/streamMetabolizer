#' MCMC estimation by JAGS with no pooling and both process and observation error
#' 
#' Compatible \code{model_file} options are 
#' \code{c('nopool_oipc_pairmeans.jags', 'nopool_oipc_Euler.jags')}.
#' 
#' @inheritParams specs_all
#'   
#' @export
specs_bayes_jags_nopool_oipc <- function(
  
  # model setup (model_path will be added in metab_bayes)
  model_file = 'nopool_oipc_pairmeans.jags',
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
  
  err_proc_acor_phi_min = 0,
  err_proc_acor_phi_max = 1,
  err_proc_acor_sigma_min = 0,
  err_proc_acor_sigma_max = 5,
  err_obs_iid_sigma_min = 0,
  err_obs_iid_sigma_max = 5,
  
  # inheritParams prepdata_bayes
  priors = FALSE,
  
  # inheritParams mcmc_bayes
  params_out = c("GPP_daily", "ER_daily", "K600_daily", "err_obs_iid_sigma", "err_proc_acor_sigma", "err_proc_acor_phi"),
  n_chains = 4, 
  n_cores = 4, 
  adapt_steps = 250, 
  burnin_steps = 250, 
  saved_steps = 500, 
  thin_steps = 1,
  verbose = FALSE
  
) {
  
  as.list(environment())
  
}
