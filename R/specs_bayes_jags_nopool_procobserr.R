#' \code{specs_bayes_jags_nopool_procobserr} - a JAGS model with no pooling and 
#' both process and observation error
#' 
#' @rdname specs_bayes
#'   
#' @param model_file character. The model definition file to use. The file may 
#'   be specified either as a file path relative to the streamMetabolizer 
#'   models/bayes directory (the first assumption; this directory can be found 
#'   with \code{system.file("models/bayes", package="streamMetabolizer")}) or as
#'   an absolute path or a path relative to the current working directory (the 
#'   second assumption, if the first assumption turns up no files of the given 
#'   name). For example, the default is \code{"jags/metab_bayes_simple.txt"}. 
#'   The containing folder (in this case \code{"jags"}) will determine which 
#'   MCMC software package is used. The file name (in this case 
#'   \code{"metab_bayes_simple.txt"}) will determine not only the model file to 
#'   use but also which variables are packaged and sent to the MCMC software.
#' @param bayes_fun character in \code{c('bayes_1ply', 'bayes_all')} indicating 
#'   whether the data should be split into daily chunks first ('bayes_1ply') or 
#'   passed to the model fitting function in one big chunk ('bayes_all')
#' @param bayes_software character in \code{c('jags','stan')} indicating the 
#'   software package to use for the MCMC process
#'   
#' @param GPP.daily.mu The mean of a dnorm distribution for GPP.daily, the daily
#'   rate of gross primary production
#' @param GPP.daily.sigma The standard deviation of a dnorm distribution for 
#'   GPP.daily, the daily rate of gross primary production
#' @param ER.daily.mu The mean of a dnorm distribution for ER.daily, the daily 
#'   rate of ecosystem respiration
#' @param ER.daily.sigma The standard deviation of a dnorm distribution for 
#'   ER.daily, the daily rate of ecosystem respiration
#' @param K600.daily.mu The mean of a dnorm distribution for K600.daily, the 
#'   daily rate of reaeration
#' @param K600.daily.sigma The standard deviation of a dnorm distribution for 
#'   K600.daily, the daily rate of reaeration
#' @param err.proc.phi.min The lower bound on a dunif distribution for 
#'   err.proc.phi, the process error autocorrelation coefficient
#' @param err.proc.phi.max The upper bound on a dunif distribution for 
#'   err.proc.phi, the process error autocorrelation coefficient
#' @param err.proc.sigma.min The lower bound on a dunif distribution for 
#'   err.proc.sigma, the standard deviation of the process error
#' @param err.proc.sigma.max The upper bound on a dunif distribution for 
#'   err.proc.sigma, the standard deviation of the process error
#' @param err.obs.sigma.min The lower bound on a dunif distribution for 
#'   err.obs.sigma, the standard deviation of the observation error
#' @param err.obs.sigma.max The upper bound on a dunif distribution for 
#'   err.obs.sigma, the standard deviation of the observation error
#'   
#' @inheritParams prepjags_bayes
#' @inheritParams runjags_bayes
#'   
#' @export
#' @family model_specs
specs_bayes_jags_nopool_procobserr <- function(
  
  # model setup (model_path will be added in metab_bayes)
  model_file = 'jags/nopool_procobserr.txt',
  bayes_fun = 'bayes_1ply',
  bayes_software = 'jags',
  
  # hyperparameters
  GPP.daily.mu = 10,
  GPP.daily.sigma = 10,
  ER.daily.mu = -10,
  ER.daily.sigma = 10,
  K600.daily.mu = 10,
  K600.daily.sigma = 10,
  
  err.proc.phi.min = 0,
  err.proc.phi.max = 1,
  err.proc.sigma.min = 0,
  err.proc.sigma.max = 0.0005,
  err.obs.sigma.min = 0,
  err.obs.sigma.max = 0.5,
  
  # inheritParams prepjags_bayes
  priors = FALSE,
  
  # inheritParams runjags_bayes
  params_out = c("GPP.daily", "ER.daily", "K600.daily", "err.obs.sigma", "err.proc.sigma", "err.proc.phi"),
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
    
    err.proc.phi.min = err.proc.phi.min,
    err.proc.phi.max = err.proc.phi.max,
    err.proc.sigma.min = err.proc.sigma.min,
    err.proc.sigma.max = err.proc.sigma.max,
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
