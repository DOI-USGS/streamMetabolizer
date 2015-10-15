#' \code{specs_bayes_stan_nopool_obserr} - a Stan model with no pooling and both
#' process and observation error. Compatible \code{model_file} options are 
#' \code{c('nopool_procobserr_pairmeans.stan', 'nopool_procobserr_Euler.stan')}.
#' 
#' @rdname specs_bayes
#'   
#' @inheritParams specs_bayes_jags_nopool_procobserr
#' @inheritParams prepdata_bayes
#' @inheritParams mcmc_bayes
#'   
#' @export
specs_bayes_stan_nopool_procobserr <- function(
  
  # model setup (model_path will be added in metab_bayes)
  model_file = 'nopool_obserr_pairmeans.stan', # or 'nopool_obserr_Euler.stan'
  bayes_fun = 'bayes_1ply',
  bayes_software = 'stan',
  
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
  max_cores = 4, 
  adapt_steps = 100, 
  burnin_steps = 40, 
  num_saved_steps = 400, 
  thin_steps = 1,
  verbose = TRUE
  
) {
  
  list(
    model_file = model_file,
    bayes_fun = bayes_fun,
    bayes_software = bayes_software,
    
    GPP_daily_mu = GPP_daily_mu,
    GPP_daily_sigma = GPP_daily_sigma,
    ER_daily_mu = ER_daily_mu,
    ER_daily_sigma = ER_daily_sigma,
    K600_daily_mu = K600_daily_mu,
    K600_daily_sigma = K600_daily_sigma,
    
    err_obs_sigma_min = err_obs_sigma_min,
    err_obs_sigma_max = err_obs_sigma_max,
    
    priors = priors,
    
    params_out = params_out,
    max_cores = max_cores, 
    adapt_steps = adapt_steps, 
    burnin_steps = burnin_steps, 
    num_saved_steps = num_saved_steps, 
    thin_steps = thin_steps,
    verbose = verbose
  )
  
}