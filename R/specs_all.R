#' List all possible arguments and the specs lists in which they occur
#' 
#' @param GPP_init the inital value of daily GPP to use in the NLM fitting 
#'   process
#' @param ER_init the inital value of daily ER to use in the NLM fitting process
#' @param K600_init the inital value of daily K600 to use in the NLM fitting 
#'   process. Ignored if K600 is supplied in data_daily, except for those dates 
#'   where K600 is NA. If there are any such dates, K600_init must have a 
#'   numeric (non-NA) value, as this will be used to estimate K600 for those 
#'   dates.
#' @inheritParams negloglik_1ply
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
#'   
#' @param err_obs_iid_sigma_min The lower bound on a dunif distribution for 
#'   err_obs_iid_sigma, the standard deviation of the observation error
#' @param err_obs_iid_sigma_max The upper bound on a dunif distribution for 
#'   err_obs_iid_sigma, the standard deviation of the observation error
#'   
#' @param err_proc_acor_phi_min lower bound on the autocorrelation coefficient
#'   for the autocorrelated component of process [& sometimes observation] error
#' @param err_proc_acor_phi_max upper bound on the autocorrelation coefficient
#'   for the autocorrelated component of process [& sometimes observation] error
#' @param err_proc_acor_sigma_min lower bound on the standard deviation of the 
#'   autocorrelated component of process [& sometimes observation] error
#' @param err_proc_acor_sigma_max upper bound on the standard deviation of the 
#'   autocorrelated component of process [& sometimes observation] error
#' @param err_proc_iid_sigma_min lower bound on the standard deviation of the 
#'   uncorrelated (IID) component of process [& sometimes observation] error
#' @param err_proc_iid_sigma_max upper bound on the standard deviation of the 
#'   uncorrelated (IID) component of process [& sometimes observation] error
#'   
#' @return a data.frame with arguments as rows, columns as spec_ functions, and 
#'   values as logicals: is this an argument to this function?
#'   
#' @inheritParams prepdata_bayes
#' @inheritParams mcmc_bayes
#' 
#' @inheritParams calc_DO_mod_w_sim_error
#' @param sim.seed NA to specify that each call to predict_DO should generate 
#'   new values, or an integer, as in the \code{seed} argument to 
#'   \code{\link{set.seed}}, specifying the seed to set before every execution
#'   of predict_DO
#'   
#' @import dplyr
#' 
#' @examples 
#' specs_table <- specs_all()
#' specs_table$spec
#' specs_table[specs_table$specs_bayes_jags_nopool_oi, 'spec']
#'   
#' @export
specs_all <- function(
  
  ## MLE
  
  # initial values
  GPP_init, 
  ER_init, 
  K600_init,
  
  # inheritParams negloglik_1ply
  calc_DO_fun,
  ODE_method,
  
  
  ## Bayes
  
  # model setup
  model_file,
  bayes_fun,
  bayes_software,
  keep_mcmcs,
  
  # hyperparameters
  GPP_daily_mu,
  GPP_daily_sigma,
  ER_daily_mu,
  ER_daily_sigma,
  K600_daily_mu,
  K600_daily_sigma,
  
  err_obs_iid_sigma_min,
  err_obs_iid_sigma_max,
  
  err_proc_acor_phi_min,
  err_proc_acor_phi_max,
  err_proc_acor_sigma_min,
  err_proc_acor_sigma_max,
  err_proc_iid_sigma_min,
  err_proc_iid_sigma_max,
  
  # inheritParams prepdata_bayes
  priors,
  
  # inheritParams mcmc_bayes
  params_out,
  n_chains,
  n_cores,
  adapt_steps,
  burnin_steps,
  saved_steps,
  thin_steps,
  verbose,
  
  
  ## Sim
  
  # inheritParams calc_DO_mod_w_sim_error
  err.obs.sigma,
  err.obs.phi, 
  err.proc.sigma, 
  err.proc.phi,
  #ODE_method, # already defined in negloglik_1ply (MLE)
  
  # sim data predictability
  sim.seed
  
) {
  
  # list the arguments
  all_specs <- names(formals(specs_all))
  
  # list the specs_ functions
  . <- '.dplyr.var'
  spec_funs <- ls("package:streamMetabolizer") %>%
    grep("^specs_", ., value=TRUE) %>% 
    grep("all|funs", ., invert=TRUE, value=TRUE)
  
  # create a data.frame with arguments as rows, columns as spec_ functions, and 
  # values as logicals: is this an argument to this function?
  bind_cols(
    data.frame(spec=all_specs, stringsAsFactors=FALSE),
    as.data.frame(
      lapply(setNames(spec_funs,spec_funs), function(fun) {
        fun_args <- names(formals(fun))
        if(any(missing_args <- !(fun_args %in% all_specs))) 
          warning("arguments to ",fun," not defined in specs_all: ", paste0(fun_args[missing_args], collapse=", "))
        sapply(all_specs, function(spec) {
          spec %in% fun_args
        })
      })
    )) %>%
    as.data.frame(stringsAsFactors=FALSE)
  
}