#' \code{specs_bayes_jags_nopool_obserr} - a JAGS model with no pooling and only
#' observation error. Compatible \code{model_file} options are
#' \code{c('jags/nopool_obserr_pairmeans.txt', 'jags/nopool_obserr_Euler.txt')}.
#' 
#' @rdname specs_bayes
#'   
#' @inheritParams specs_bayes_jags_nopool_procobserr
#' @inheritParams prepjags_bayes
#' @inheritParams runjags_bayes
#'   
#' @export
specs_bayes_jags_nopool_obserr <- function(
  
  # model setup (model_path will be added in metab_bayes)
  model_file = 'jags/nopool_obserr_pairmeans.txt', # or 'jags/nopool_obserr_Euler.txt'
  bayes_fun = 'bayes_1ply',
  bayes_software = 'jags',
  
  # hyperparameters
  GPP.daily.mu = 10,
  GPP.daily.sigma = 10,
  ER.daily.mu = -10,
  ER.daily.sigma = 10,
  K600.daily.mu = 10,
  K600.daily.sigma = 10,
  
  err.obs.sigma.min = 0,
  err.obs.sigma.max = 0.5,
  
  # inheritParams prepjags_bayes
  priors = FALSE,
  
  # inheritParams runjags_bayes
  params_out = c("GPP.daily", "ER.daily", "K600.daily", "err.obs.sigma"),
  max_cores = 4, 
  adapt_steps = 100, 
  burnin_steps = 40, 
  num_saved_steps = 400, 
  thin_steps = 1
  
) {
  
  list(
    model_file = model_file,
    bayes_fun = bayes_fun,
    bayes_software = bayes_software,
    
    GPP.daily.mu = GPP.daily.mu,
    GPP.daily.sigma = GPP.daily.sigma,
    ER.daily.mu = ER.daily.mu,
    ER.daily.sigma = ER.daily.sigma,
    K600.daily.mu = K600.daily.mu,
    K600.daily.sigma = K600.daily.sigma,
    
    err.obs.sigma.min = err.obs.sigma.min,
    err.obs.sigma.max = err.obs.sigma.max,
    
    priors = priors,
    
    params_out = params_out,
    max_cores = max_cores, 
    adapt_steps = adapt_steps, 
    burnin_steps = burnin_steps, 
    num_saved_steps = num_saved_steps, 
    thin_steps = thin_steps
  )
  
}