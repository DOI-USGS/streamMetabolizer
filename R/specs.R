#' Generate a coherent list of model specs
#' 
#' Generates an internally consistent list of model specifications that may be 
#' passed to \code{metab_bayes}, \code{metab_mle}, etc. via the 
#' \code{model_specs} argument. This help file gives the definitive list of all 
#' possible model specs, but only a subset of these are relevant to any given 
#' \code{model_name}. See the 'Relevant arguments' section below. Irrelevant 
#' arguments for the given \code{model_name} should not be explicitly passed 
#' into this function (but don't worry - we'll just stop and tell you if you 
#' make a mistake). Relevant arguments for the given \code{model_name} either 
#' have default values or do not (see Usage). Relevant arguments without a 
#' default should rarely be overridden, because their values will be determined 
#' based on other arguments. Relevant arguments that do have a default can, and 
#' often should, be overridden to tailor the model to your needs. Relevant 
#' arguments whose default is a vector, e.g.,
#' 
#' @section Relevant arguments:
#'   
#'   * metab_bayes: Always relevant: \code{model_name, bayes_fun, engine,
#'   keep_mcmcs, GPP_daily_mu, GPP_daily_sigma, ER_daily_mu, ER_daily_sigma,
#'   K600_daily_mu, K600_daily_sigma, priors, params_out, n_chains, n_cores,
#'   burnin_steps, saved_steps, thin_steps, verbose}. If 
#'   \code{mm_parse_name(model_name)$err_obs_iid} then also 
#'   \code{err_obs_iid_sigma_min, err_obs_iid_sigma_max}. If 
#'   \code{mm_parse_name(model_name)$err_proc_acor} then also 
#'   \code{err_proc_acor_phi_min, err_proc_acor_phi_max, 
#'   err_proc_acor_sigma_min, err_proc_acor_sigma_max}. If 
#'   \code{mm_parse_name(model_name)$err_proc_iid} then also 
#'   \code{err_proc_iid_sigma_min, err_proc_iid_sigma_max}. If 
#'   \code{features$engine == 'jags'} then also \code{adapt_steps},
#'   
#'   * metab_mle: \code{model_name, calc_DO_fun, ODE_method, GPP_init, ER_init, 
#'   K600_init}
#'   
#'   * metab_night: \code{model_name}
#'   
#'   * metab_Kmodel: \code{model_name, engine, weights, predictors, filters, 
#'   transforms}
#'   
#'   * metab_sim: \code{model_name, err.obs.sigma, err.obs.phi, err.proc.sigma, 
#'   err.proc.phi, ODE_method, sim.seed}
#'   
#' @param model_name character string identifying the model features. Use 
#'   \code{\link{mm_name}} for valid names. This may be a full model file path 
#'   for custom Bayesian models, as long as basename(model_name) can still be 
#'   parsed correctly with \code{mm_parse_name()}. In that case the file may be 
#'   specified either as a file path relative to the streamMetabolizer models 
#'   directory (the first assumption; this directory can be found with 
#'   \code{system.file("models", package="streamMetabolizer")}) or as an 
#'   absolute path or a path relative to the current working directory (the 
#'   second assumption, if the first assumption turns up no files of the given 
#'   name).
#'   
#' @param engine The software or function to use in fitting the model. Should be
#'   specified via \code{mm_name} rather than here. For type='bayes', one of
#'   \code{c('jags','stan')} indicating the software package to use for the MCMC
#'   process. For type='Kmodel', the name of an interpolation or regression
#'   method relating K to the predictor[s] of choice. One of \code{c("mean",
#'   "lm", "loess")}. For types in \code{c('mle','night','sim')} there's only
#'   one option so it's not included in \code{specs()} (but is nonetheless noted
#'   in the suffix of the model name, e.g., \code{"m_np_oi_pm_km.nlm"} uses
#'   \code{nlm()} for model fitting)
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
#' @param bayes_fun character in \code{c('bayes_1ply', 'bayes_all')} indicating 
#'   whether the data should be split into daily chunks first ('bayes_1ply') or 
#'   passed to the model fitting function in one big chunk ('bayes_all')
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
#' @inheritParams prepdata_Kmodel
#' @inheritParams Kmodel_allply
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
#' @return an internally consistent list of arguments that may be passed to 
#'   \code{metab_bayes}, \code{metab_mle}, etc. as the \code{model_specs} 
#'   argument
#'   
#' @examples 
#' specs(mm_name(type='bayes', err_obs_iid=TRUE, err_proc_acor=TRUE))
#'   
#' @export
specs <- function(
  
  ## All or several models
  
  model_name = mm_name(),
  engine,
  
  
  ## MLE
  
  # initial values
  GPP_init = 10, 
  ER_init = -10, 
  K600_init = 10,
  
  # inheritParams negloglik_1ply
  calc_DO_fun,
  ODE_method,
  
  
  ## Bayes
  
  # model setup
  bayes_fun,
  keep_mcmcs = TRUE,
  
  # hyperparameters
  GPP_daily_mu = 10,
  GPP_daily_sigma = 10,
  ER_daily_mu = -10,
  ER_daily_sigma = 10,
  K600_daily_mu = 10,
  K600_daily_sigma = 10,
  
  err_obs_iid_sigma_min = 0,
  err_obs_iid_sigma_max = 5,
  err_proc_acor_phi_min = 0,
  err_proc_acor_phi_max = 1,
  err_proc_acor_sigma_min = 0,
  err_proc_acor_sigma_max = 5,
  err_proc_iid_sigma_min = 0,
  err_proc_iid_sigma_max = 5,
  
  # inheritParams prepdata_bayes
  priors = FALSE,
  
  # inheritParams mcmc_bayes
  params_out,
  n_chains = 4,
  n_cores = 4,
  adapt_steps = switch(mm_parse_name(model_name)$engine, jags=250, NA),
  burnin_steps = switch(mm_parse_name(model_name)$engine, stan=500, jags=250, NA),
  saved_steps = switch(mm_parse_name(model_name)$engine, stan=500, jags=1000, NA),
  thin_steps = 1,
  verbose = FALSE,
  
  
  ## Kmodel
  
  #inheritParams prepdata_Kmodel
  weights = c("K600/CI"),
  filters = c(CI.max=NA, discharge.daily.max=NA, velocity.daily.max=NA),
  
  #inheritParams Kmodel_allply
  predictors = c("discharge.daily"), 
  transforms = c(K600='log', date=NA, velocity.daily="log", discharge.daily="log"),
  
  
  ## Sim
  
  # inheritParams calc_DO_mod_w_sim_error
  err.obs.sigma = 0.1,
  err.obs.phi = 0,
  err.proc.sigma = 0,
  err.proc.phi = 0,
  
  # declared already (for MLE): ODE_method
  
  # sim data predictability
  sim.seed = NA
  
) {
  
  # collect info about the arguments
  required <- 'model_name'
  all_possible <- names(formals(specs))
  not_missing <- names(as.list(match.call())[-1]) # the arguments that were given explicitly
  yes_missing <- all_possible[!(all_possible %in% not_missing)]
  prefer_missing <- all_possible[sapply(formals(specs), is.symbol)] # the arguments w/o defaults
  
  # argument checks
  if(any(required %in% yes_missing))
    stop("missing and required argument: ", paste(required[required %in% yes_missing], collapse=", "))
  if(length(redundant <- not_missing[not_missing %in% prefer_missing]) > 0) 
    warning("argument[s] that should usually not be specified: ", paste(redundant, collapse=", "))
  
  # check the validity of the model_name against the list of officially accepted model names
  mm_validate_name(model_name)
  
  # parse the model_name
  features <- mm_parse_name(model_name)
  
  # logic checks
  if(!(features$pooling %in% c('none'))) stop("models with pooling are not implemented yet")
  
  # collect the defaults + directly specified arguments
  all_specs <- as.list(environment())
  
  # copy/calculate arguments as appropriate to the model
  specs <- list()
  if(features$type == 'bayes') {
    
    # list all needed arguments
    included <- c(
      # model setup
      'model_name', 'bayes_fun', 'engine', 'keep_mcmcs',
      
      # hyperparameters
      'GPP_daily_mu', 'GPP_daily_sigma', 'ER_daily_mu', 'ER_daily_sigma', 'K600_daily_mu', 'K600_daily_sigma',
      if(features$err_obs_iid) c('err_obs_iid_sigma_min', 'err_obs_iid_sigma_max'),
      if(features$err_proc_acor) c('err_proc_acor_phi_min', 'err_proc_acor_phi_max', 'err_proc_acor_sigma_min', 'err_proc_acor_sigma_max'),
      if(features$err_proc_iid) c('err_proc_iid_sigma_min', 'err_proc_iid_sigma_max'),
      
      # inheritParams prepdata_bayes
      'priors',
      
      # inheritParams mcmc_bayes
      'params_out', 'n_chains', 'n_cores', 
      if(features$engine == 'jags') 'adapt_steps', 
      'burnin_steps', 'saved_steps', 'thin_steps', 'verbose'
    )

    # compute some arguments
    if('bayes_fun' %in% yes_missing) {
      all_specs$bayes_fun <- switch(features$pooling, none='bayes_1ply', NA)
    }
    if('engine' %in% yes_missing) {
      all_specs$engine <- features$engine
    }
    if('params_out' %in% yes_missing) {
      all_specs$params_out <- c(
        c('GPP_daily','ER_daily','K600_daily'), 
        if(features$err_obs_iid) 'err_obs_iid_sigma',
        if(features$err_proc_acor) c('err_proc_acor_phi', 'err_proc_acor_sigma'),
        if(features$err_proc_iid) 'err_proc_iid_sigma')
    }
    
    # check for errors/inconsistencies
    model_path <- system.file(paste0("models/", model_name), package="streamMetabolizer")
    if(!file.exists(model_path)) 
      model_path <- model_name
    if(!file.exists(model_path)) 
      warning(suppressWarnings(paste0("could not locate the model file at ", model_path)))
    if(features$engine == "NA") 
      stop('engine must be specified for Bayesian models')
    
  } else if(features$type == 'mle') {
    # list all needed arguments
    included <- c('model_name', 'calc_DO_fun', 'ODE_method', 'GPP_init', 'ER_init', 'K600_init')
    
    if('calc_DO_fun' %in% yes_missing) {
      all_specs$calc_DO_fun <- if(features$err_obs_iid && !features$err_proc_iid) calc_DO_mod else calc_DO_mod_by_diff
    }
    if('ODE_method' %in% yes_missing) {
      all_specs$ODE_method <- features$ode_method
    }
    
  } else if(features$type == 'night') {
    # list all needed arguments
    included <- c('model_name')
    
  } else if(features$type == 'Kmodel') {
    # list all needed arguments
    included <- c('model_name', 'engine', 'weights', 'filters', 'predictors', 'transforms')
    
    if('engine' %in% yes_missing) {
      all_specs$engine <- features$engine
    }

  } else if(features$type == 'sim') {
    # list all needed arguments
    included <- c('model_name', 'err.obs.sigma', 'err.obs.phi', 'err.proc.sigma', 'err.proc.phi', 'ODE_method', 'sim.seed')
    
    if('ODE_method' %in% yes_missing) {
      all_specs$ODE_method <- features$ode_method
    }

  }
  
  # stop if truly irrelevant arguments were given
  if(length(irrelevant <- not_missing[!(not_missing %in% included)]) > 0) 
    stop("irrelevant argument: ", paste(irrelevant, collapse=", "))
  
  # return just the arguments we actually need
  all_specs[included]
  
}