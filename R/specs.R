#' Generate a coherent list of model specs
#' 
#' Generates an internally consistent list of model specifications that may be 
#' passed to \code{metab_bayes}, \code{metab_mle}, etc. via the \code{specs} 
#' argument. This help file gives the definitive list of all possible model 
#' specs, but only a subset of these are relevant to any given 
#' \code{model_name}. See the 'Relevant arguments' section below. Irrelevant 
#' arguments for the given \code{model_name} should not be explicitly passed 
#' into this function (but don't worry - we'll just stop and tell you if you 
#' make a mistake). Relevant arguments for the given \code{model_name} either 
#' have default values or do not (see Usage). Relevant arguments without a 
#' default should rarely be overridden, because their values will be determined 
#' based on other arguments. Relevant arguments that do have a default can, and 
#' often should, be overridden to tailor the model to your needs.
#' 
#' @section Relevant arguments:
#'   
#'   * metab_bayes: Always relevant: \code{model_name, engine, split_dates, 
#'   keep_mcmcs, day_start, day_end, day_tests, GPP_daily_mu, GPP_daily_sigma, 
#'   ER_daily_mu, ER_daily_sigma, priors, params_out, n_chains, n_cores, 
#'   burnin_steps, saved_steps, thin_steps, verbose}. The need for other 
#'   arguments depends on features of the model structure, as from 
#'   \code{mm_parse_name(model_name)}: If \code{$pool_K600=='none'} then 
#'   \code{K600_daily_mu, K600_daily_sigma}. If \code{$pool_K600=='normal'} then
#'   \code{K600_daily_mu_mu, K600_daily_mu_sigma, K600_daily_sigma_shape, 
#'   K600_daily_sigma_rate}. If \code{pool_K600=='linear'} then 
#'   \code{K600_daily_beta_mu, K600_daily_beta_sigma, K600_daily_sigma_shape, 
#'   K600_daily_sigma_rate}. If \code{pool_K600=='binned'} then 
#'   \code{K600_daily_beta_num, K600_daily_beta_cuts, K600_daily_beta_mu, 
#'   K600_daily_beta_sigma, K600_daily_sigma_shape, K600_daily_sigma_rate}. If 
#'   \code{err_obs_iid} then \code{err_obs_iid_sigma_shape, 
#'   err_obs_iid_sigma_rate}. If \code{err_proc_acor} then 
#'   \code{err_proc_acor_phi_shape, err_proc_acor_phi_rate, 
#'   err_proc_acor_sigma_shape, err_proc_acor_sigma_rate}. If 
#'   \code{err_proc_iid} then \code{err_proc_iid_sigma_shape, 
#'   err_proc_iid_sigma_rate}. If \code{engine == 'jags'} then 
#'   \code{adapt_steps}.
#'   
#'   * metab_mle: \code{model_name, day_start, day_end, day_tests, calc_DO_fun, 
#'   ODE_method, GPP_init, ER_init, K600_init}
#'   
#'   * metab_night: \code{model_name, day_start, day_end, day_tests}
#'   
#'   * metab_Kmodel: \code{model_name, engine, day_start, day_end, day_tests, 
#'   weights, filters, predictors, transforms, other_args}. Note that the 
#'   defaults for \code{weights}, \code{predictors}, \code{filters}, and 
#'   \code{transforms} are adjusted according to the \code{engine} implied by 
#'   \code{model_name}.
#'   
#'   * metab_sim: \code{model_name, day_start, day_end, day_tests, 
#'   err.obs.sigma, err.obs.phi, err.proc.sigma, err.proc.phi, ODE_method, 
#'   sim.seed}
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
#' @param engine The software or function to use in fitting the model. Should be
#'   specified via \code{mm_name} rather than here. For type='bayes', one of 
#'   \code{c('jags','stan')} indicating the software package to use for the MCMC
#'   process. For type='Kmodel', the name of an interpolation or regression 
#'   method relating K to the predictor[s] of choice. One of \code{c("mean", 
#'   "lm", "loess")}. For types in \code{c('mle','night','sim')} there's only 
#'   one option so it's not included in \code{specs()} (but is nonetheless noted
#'   in the suffix of the model name, e.g., \code{"m_np_oi_pm_km.nlm"} uses 
#'   \code{nlm()} for model fitting)
#' @inheritParams mm_model_by_ply
#' @inheritParams mm_is_valid_day
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
#' @param split_dates logical indicating whether the data should be split into 
#'   daily chunks first (TRUE) or processed within one big model (FALSE). If 
#'   valid days differ in their timestep length, split_dates will need to be 
#'   TRUE; otherwise, FALSE is generally more efficient. FALSE is also the only 
#'   appropriate solution for a hierarchical model that pools information on 
#'   error, K600, etc. across days.
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
#' @param K600_daily_mu Applies when pool_K600 is 'none'. The mean of a dnorm 
#'   distribution for K600_daily, the daily rate of reaeration
#' @param K600_daily_sigma Applies when pool_K600 is 'none'. The standard 
#'   deviation of a dnorm distribution for K600_daily, the daily rate of 
#'   reaeration
#'   
#' @param K600_daily_mu_mu hyperparameter for pool_K600='normal'. The mean 
#'   parameter (mu_mu) of a normal distribution of mu in K ~ N(mu, sigma), mu ~ 
#'   N(mu_mu, mu_sigma)
#' @param K600_daily_mu_sigma hyperparameter for pool_K600='normal'. The 
#'   standard deviation parameter (mu_sigma) of a normal distribution of mu in K
#'   ~ N(mu, sigma), mu ~ N(mu_mu, mu_sigma)
#'   
#' @param K600_daily_beta_mu hyperparameter for pool_K600 in 
#'   c('linear','binned'). The means of prior distributions for the 
#'   K600_daily_beta paramters. For pool_K600='linear', there are 2 betas 
#'   corresponding to the intercept and slope, respectively, of the linear model
#'   \code{log(K600) ~ K600_daily_beta[1] + K600_daily_beta[2]*log(Q)}. For 
#'   pool_K600='binned', there are K600_daily_beta_num betas each giving the 
#'   predicted K600 when Q_daily is in the corresponding bin (see 
#'   \code{K600_daily_beta_cuts}).
#' @param K600_daily_beta_sigma hyperparameter for pool_K600='linear'. The 
#'   standard deviation parameter of a normally distributed intercept term 
#'   (beta0) in the linear model K ~ N(beta0 + beta1*log(Q)), beta0 ~ 
#'   N(beta0_mu, beta0_sigma)
#'   
#' @param K600_daily_beta_num hyperparameter for pool_K600='binned'. The number 
#'   of bins into which daily discharge values should be grouped. Each bin 
#'   predicts a single value of K600_daily_pred, such that any day on which 
#'   \code{discharge_bin_daily} equals that bin will have \code{K600_daily ~ 
#'   N(K600_daily_beta[discharge_bin_daily], K600_daily_sigma)}
#' @param K600_daily_beta_cuts hyperparameter for pool_K600='binned'. Either (1)
#'   character of length 1 in c('number','interval') indicating how the bin cuts
#'   should be determined, or (2) numeric (as in \code{breaks} in 
#'   \code{\link[base]{cut}}) of length K600_daily_beta_num+1 giving the 
#'   natural-log-space breakpoints defining the bins. For option 1, the 
#'   implementation uses or is equivalent to the corresponding functions 
#'   \code{\link[ggplot2]{cut_interval}} (to cut into bins having equal numeric 
#'   ranges in natural log space) and \code{\link[ggplot2]{cut_number}} (to cut 
#'   into bins having ~equal numbers of ln_discharge_daily observations). For 
#'   option 2, make sure to include the full range of ln_discharge_daily, with 
#'   the first value smaller than all ln_discharge_daily values and the last 
#'   value greater than or equal to all ln_discharge_daily values.
#'   
#' @param K600_daily_sigma_shape hyperparameter for pool_K600 in 
#'   c('normal','linear','binned'). The shape (= alpha = k) parameter of a gamma
#'   distribution of sigma in K ~ N(mu, sigma), sigma ~ gamma(shape, rate)
#' @param K600_daily_sigma_rate hyperparameter for pool_K600 in 
#'   c('normal','linear','binned'). The rate (= beta = 1/theta = inverse scale) 
#'   parameter of a gamma distribution of sigma in K ~ N(mu, sigma), sigma ~ 
#'   gamma(shape, rate)
#'   
#' @param err_obs_iid_sigma_shape The shape parameter on a gamma distribution 
#'   for err_obs_iid_sigma, the standard deviation of the observation error
#' @param err_obs_iid_sigma_rate The rate parameter on a gamma distribution for 
#'   err_obs_iid_sigma, the standard deviation of the observation error
#' @param err_proc_acor_phi_shape The shape parameter on a gamma distribution 
#'   for err_proc_acor_phi, the autocorrelation coefficient for the 
#'   autocorrelated component of process [& sometimes observation] error
#' @param err_proc_acor_phi_rate The rate parameter on a gamma distribution for 
#'   err_proc_acor_phi, the autocorrelation coefficient for the autocorrelated 
#'   component of process [& sometimes observation] error
#' @param err_proc_acor_sigma_shape The shape parameter on a gamma distribution 
#'   for err_proc_acor_sigma, the standard deviation of the autocorrelated 
#'   component of process [& sometimes observation] error
#' @param err_proc_acor_sigma_rate The rate parameter on a gamma distribution 
#'   for err_proc_acor_sigma, the standard deviation of the autocorrelated 
#'   component of process [& sometimes observation] error
#' @param err_proc_iid_sigma_shape The shape parameter on a gamma distribution 
#'   for err_proc_iid_sigma, the standard deviation of the uncorrelated (IID) 
#'   component of process [& sometimes observation] error
#' @param err_proc_iid_sigma_rate The rate parameter on a gamma distribution for
#'   err_proc_iid_sigma, the standard deviation of the uncorrelated (IID) 
#'   component of process [& sometimes observation] error
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
#'   \code{metab_bayes}, \code{metab_mle}, etc. as the \code{specs} argument
#'   
#' @examples
#' specs(mm_name(type='mle', err_obs_iid=FALSE, err_proc_iid=TRUE))
#' specs(mm_name(type='bayes', pool_K600='normal'))
#' @export
specs <- function(
  
  ## All or several models
  
  model_name = mm_name(),
  engine,
  
  # inheritParams mm_model_by_ply
  day_start = 4,
  day_end = 28,
  
  # inheritParams mm_is_valid_day
  day_tests=c('full_day', 'even_timesteps', 'complete_data'),
  
  
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
  split_dates,
  keep_mcmcs = TRUE,
  
  # hyperparameters for non-hierarchical GPP & ER
  GPP_daily_mu = 10,
  GPP_daily_sigma = 10,
  ER_daily_mu = -10,
  ER_daily_sigma = 10,
  
  # hyperparameters for non-hierarchical K600
  K600_daily_mu = 10,
  K600_daily_sigma = 10,
  
  # hyperparameters for hierarchical K600 - normal
  K600_daily_mu_mu = 10,
  K600_daily_mu_sigma = 10,
  
  # hyperparameters for hierarchical K600 - linear. defaults should be reasonably
  # constrained, not too wide. element names are ignored
  K600_daily_beta_mu = c(intercept=10, slope=3),
  K600_daily_beta_sigma = c(intercept=8, slope=2),
  
  # hyperparameters for hierarchical K600 - binned. K ~ beta[Q_bin_daily]
  # (constant for each binned log(Q), i.e., rectangular interpolation among
  # bins). beta_mu and beta_sigma may be length K600_daily_beta_num or length 1
  # (to be replicated to length K600_daily_beta_num)
  K600_daily_beta_num = 5,
  K600_daily_beta_cuts = 'number',
  # K600_daily_beta_mu = rep(10, K600_daily_beta_num), # already declared for linear hierarchy above
  # K600_daily_beta_sigma = rep(10, K600_daily_beta_num), # already declared for linear hierarchy above
  
  # hyperparameters for any hierarchical K600
  K600_daily_sigma_shape = 1,
  K600_daily_sigma_rate = 2,
  
  # hyperparameters for error terms
  err_obs_iid_sigma_shape = 1,
  err_obs_iid_sigma_rate = 10,
  err_proc_acor_phi_shape = 1,
  err_proc_acor_phi_rate = 1000,
  err_proc_acor_sigma_shape = 1,
  err_proc_acor_sigma_rate = 1000,
  err_proc_iid_sigma_shape = 1,
  err_proc_iid_sigma_rate = 100,
  
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
  weights = c("K600/CI"), # 'K600/CI' is argued for in stream_metab_usa issue #64
  filters = c(CI.max=NA, discharge.daily.max=NA, velocity.daily.max=NA),
  
  #inheritParams Kmodel_allply
  predictors = c("discharge.daily"), 
  transforms = c(K600='log', date=NA, velocity.daily="log", discharge.daily="log"),
  other_args = c(),
  
  
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
  redundant <- not_missing[not_missing %in% prefer_missing]
  if('engine' %in% redundant) {
    warning("'engine' should be specified in mm_name() rather than specs()")
    redundant <- redundant[redundant != 'engine']
  }
  if(length(redundant) > 0) {
    warning("argument[s] that should usually not be specified in specs(): ", paste(redundant, collapse=", "))
  }
  
  
  # check the validity of the model_name against the list of officially accepted model names
  mm_validate_name(model_name)
  
  # parse the model_name
  features <- mm_parse_name(model_name)
  
  # collect the defaults + directly specified arguments
  all_specs <- as.list(environment())
  
  # copy/calculate arguments as appropriate to the model
  specs <- list()
  switch(
    features$type,
    'bayes' = {
      
      # list all needed arguments
      included <- c(
        # model setup
        'model_name', 'engine', 'split_dates', 'keep_mcmcs',
        
        # date ply day_tests
        'day_start', 'day_end', 'day_tests',
        
        # hyperparameters - this section should be identical to the
        # hyperparameters section of prepdata_bayes
        c('GPP_daily_mu','GPP_daily_sigma','ER_daily_mu','ER_daily_sigma'),
        switch(
          features$pool_K600,
          none=c('K600_daily_mu', 'K600_daily_sigma'),
          normal=c('K600_daily_mu_mu', 'K600_daily_mu_sigma', 'K600_daily_sigma_shape', 'K600_daily_sigma_rate'),
          linear=c('K600_daily_beta_mu', 'K600_daily_beta_sigma', 'K600_daily_sigma_shape', 'K600_daily_sigma_rate'),
          binned=c('K600_daily_beta_num', 'K600_daily_beta_cuts', 'K600_daily_beta_mu', 'K600_daily_beta_sigma', 'K600_daily_sigma_shape', 'K600_daily_sigma_rate')),
        if(features$err_obs_iid) c('err_obs_iid_sigma_shape', 'err_obs_iid_sigma_rate'),
        if(features$err_proc_acor) c('err_proc_acor_phi_shape', 'err_proc_acor_phi_rate', 'err_proc_acor_sigma_shape', 'err_proc_acor_sigma_rate'),
        if(features$err_proc_iid) c('err_proc_iid_sigma_shape', 'err_proc_iid_sigma_rate'),
        
        # inheritParams prepdata_bayes
        'priors',
        
        # inheritParams mcmc_bayes
        'params_out', 'n_chains', 'n_cores', 
        if(features$engine == 'jags') 'adapt_steps', 
        'burnin_steps', 'saved_steps', 'thin_steps', 'verbose'
      )
      
      # compute some arguments
      if('engine' %in% yes_missing) {
        all_specs$engine <- features$engine
      }
      if('split_dates' %in% yes_missing) {
        all_specs$split_dates <- ifelse(
          features$pool_K600 %in% 'none', FALSE, # pretty sure FALSE is faster. also allows hierarchical error terms
          ifelse(features$pool_K600 %in% c('normal','linear','binned'), FALSE, 
                 stop("unknown pool_K600; unsure how to set split_dates")))
      }
      if(features$pool_K600 == 'binned') {
        # defaults are for linear pool_K600 & need adjustment for binned method
        all_specs$K600_daily_beta_mu <- rep(10, all_specs$K600_daily_beta_num)
        all_specs$K600_daily_beta_sigma <- rep(10, all_specs$K600_daily_beta_num)
      }
      if('params_out' %in% yes_missing) {
        all_specs$params_out <- c(
          c('GPP_daily','ER_daily','K600_daily'),
          switch(
            features$pool_K600,
            none=c(),
            normal=c('K600_daily_mu', 'K600_daily_sigma'),
            linear=c('K600_daily_beta', 'K600_daily_sigma'),
            binned=c('K600_daily_beta', 'K600_daily_sigma')), 
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
      
    },
    'mle' = {
      # list all needed arguments
      included <- c('model_name', 'day_start', 'day_end', 'day_tests', 'calc_DO_fun', 'ODE_method', 'GPP_init', 'ER_init', 'K600_init')
      
      if('calc_DO_fun' %in% yes_missing) {
        all_specs$calc_DO_fun <- if(features$err_obs_iid && !features$err_proc_iid) calc_DO_mod else calc_DO_mod_by_diff
      }
      if('ODE_method' %in% yes_missing) {
        all_specs$ODE_method <- features$ode_method
      }
      
    }, 
    'night' = {
      # list all needed arguments
      included <- c('model_name', 'day_start', 'day_end', 'day_tests')
      
      # some different defaults for night relative to other models
      if('day_start' %in% yes_missing) {
        all_specs$day_start <- 12
      }
      if('day_end' %in% yes_missing) {
        all_specs$day_end <- 36
      }
      if('day_tests' %in% yes_missing) {
        all_specs$day_tests <- c(day_tests, 'include_sunset')
      }
      
    }, 
    'Kmodel' = {
      # list all needed arguments
      included <- c(
        'model_name', 'engine', 'day_start', 'day_end', 'day_tests',
        'weights', 'filters', 'predictors', 'transforms', 'other_args')
      
      if('engine' %in% yes_missing) {
        all_specs$engine <- features$engine
      }
      
      # some different defaults for each engine, because no one set of defaults
      # makes sense for all engines
      #if('weights' %in% yes_missing) all_specs$weights <- c("K600/CI") # same for all, so use default as in Usage
      switch(
        all_specs$engine,
        mean={
          if('filters' %in% yes_missing) all_specs['filters'] <- list(c()) # need special syntax to assign c(). see http://stackoverflow.com/a/7945259/3203184
          if('predictors' %in% yes_missing) all_specs['predictors'] <- list(c())
          if('transforms' %in% yes_missing) all_specs$transforms <- c(K600='log')
          if('other_args' %in% yes_missing) all_specs$other_args <- list(possible_args=NULL)
        },
        lm={
          if('filters' %in% yes_missing) all_specs$filters <- c(CI.max=NA, discharge.daily.max=NA)
          if('predictors' %in% yes_missing) all_specs$predictors <- c("discharge.daily")
          if('transforms' %in% yes_missing) all_specs$transforms <- c(K600='log', discharge.daily="log")
          if('other_args' %in% yes_missing) all_specs$other_args <- list(possible_args=names(formals(lm))[-which(names(formals(lm)) %in% c('formula','data','weights'))])
        },
        loess={
          if('filters' %in% yes_missing) all_specs$filters <- c(CI.max=NA, discharge.daily.max=NA, velocity.daily.max=NA)
          if('predictors' %in% yes_missing) all_specs$predictors <- c("date", "discharge.daily")
          if('transforms' %in% yes_missing) all_specs$transforms <- c(K600='log', date=NA, velocity.daily="log", discharge.daily="log")
          if('other_args' %in% yes_missing) all_specs$other_args <- list(possible_args=names(formals(loess))[-which(names(formals(loess)) %in% c('formula','data','weights'))])
        }
      )
      
    },
    'sim' = {
      # list all needed arguments
      included <- c(
        'model_name', 'day_start', 'day_end', 'day_tests',
        'err.obs.sigma', 'err.obs.phi', 'err.proc.sigma', 'err.proc.phi',
        'ODE_method', 'sim.seed')
      
      if('ODE_method' %in% yes_missing) {
        all_specs$ODE_method <- features$ode_method
      }
    }
  )
  
  # stop if truly irrelevant arguments were given
  if(length(irrelevant <- not_missing[!(not_missing %in% included)]) > 0) 
    stop("irrelevant argument: ", paste(irrelevant, collapse=", "))
  
  # return just the arguments we actually need
  add_specs_class(all_specs[included])
  
}