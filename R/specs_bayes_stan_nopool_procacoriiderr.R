#' \code{specs_bayes_stan_nopool_obserr} - a Stan model with no pooling and both
#' process and observation error. Compatible \code{model_file} options are 
#' \code{c('nopool_procobserr_pairmeans.stan', 'nopool_procobserr_Euler.stan')}.
#' 
#' @rdname specs_bayes
#'   
#' @param err_proc_acor_phi_min lower bound on the autocorrelation coefficient
#'   for the autocorrelated component of process & observation error
#' @param err_proc_acor_phi_max upper bound on the autocorrelation coefficient
#'   for the autocorrelated component of process & observation error
#' @param err_proc_acor_sigma_min lower bound on the standard deviation of the 
#'   autocorrelated component of process & observation error
#' @param err_proc_acor_sigma_max upper bound on the standard deviation of the 
#'   autocorrelated component of process & observation error
#' @param err_proc_iid_sigma_min lower bound on the standard deviation of the 
#'   uncorrelated (IID) component of process & observation error
#' @param err_proc_iid_sigma_max upper bound on the standard deviation of the 
#'   uncorrelated (IID) component of process & observation error
#' @inheritParams specs_bayes_jags_nopool_procobserr
#' @inheritParams prepdata_bayes
#' @inheritParams mcmc_bayes
#'   
#' @export
specs_bayes_stan_nopool_procacoriiderr <- function(
  
  # model setup (model_path will be added in metab_bayes)
  model_file = 'nopool_procacoriiderr_pairmeans.stan',
  bayes_fun = 'bayes_1ply',
  bayes_software = 'stan',
  keep_mcmcs = FALSE,
  
  # hyperparameters
  GPP_daily_mu = 10,
  GPP_daily_sigma = 10,
  ER_daily_mu = -10,
  ER_daily_sigma = 10,
  K600_daily_mu = 10,
  K600_daily_sigma = 10,
  
  err_proc_acor_phi_min = 0,
  err_proc_acor_phi_max = 1, # ??? look up BUGS::car.normal. could err_proc_acor_phi be always 1, with non-1 part in err_proc_iid?
  err_proc_acor_sigma_min = 0,
  err_proc_acor_sigma_max = 5, # sdp~dunif(0,0.1)
  err_proc_iid_sigma_min = 0,
  err_proc_iid_sigma_max = 50, # sdo~dunif(0,0.1)
  
  # inheritParams prepdata_bayes
  priors = FALSE,
  
  # inheritParams mcmc_bayes
  params_out = c("GPP_daily", "ER_daily", "K600_daily", "err_proc_acor_phi", "err_proc_acor_sigma", "err_proc_iid_sigma"), #"DO_mod_1", 
  n_chains = 4, 
  n_cores = 1, 
  burnin_steps = 500, 
  num_saved_steps = 500, 
  thin_steps = 1,
  verbose = FALSE
  
) {
  
  as.list(environment())
  
}