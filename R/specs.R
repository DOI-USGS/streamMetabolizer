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
#'   keep_mcmcs, keep_mcmc_data, day_start, day_end, day_tests, ER_daily_mu,
#'   ER_daily_sigma, params_in, params_out, n_chains, n_cores, burnin_steps,
#'   saved_steps, thin_steps, verbose}. The need for other arguments depends on
#'   features of the model structure, as from \code{mm_parse_name(model_name)}:
#'   \itemize{ \item If \code{GPP_fun=='linlight'} then \code{GPP_daily_mu,
#'   GPP_daily_sigma}, while if \code{GPP_fun=='satlight'} then
#'   \code{alpha_meanlog, alpha_sdlog, Pmax_mu, Pmax_sigma}. \item If
#'   \code{pool_K600=='none'} then \code{K600_daily_meanlog, K600_daily_sdlog}.
#'   \item If \code{pool_K600=='normal'} then \code{K600_daily_meanlog_meanlog,
#'   K600_daily_meanlog_sdlog, K600_daily_sdlog_sigma}. \item If
#'   \code{pool_K600=='linear'} then \code{lnK600_lnQ_intercept_mu,
#'   lnK600_lnQ_intercept_sigma, lnK600_lnQ_slope_mu, lnK600_lnQ_slope_sigma,
#'   K600_daily_sigma_sigma}. \item If \code{pool_K600=='binned'} then
#'   \code{K600_lnQ_nodes_centers, K600_lnQ_nodediffs_sdlog,
#'   K600_lnQ_nodes_meanlog, K600_lnQ_nodes_sdlog, K600_daily_sigma_sigma}.
#'   \item If \code{err_obs_iid} then \code{err_obs_iid_sigma_scale}. \item If
#'   \code{err_proc_acor} then \code{err_proc_acor_phi_alpha,
#'   err_proc_acor_phi_beta, err_proc_acor_sigma_scale}. \item If
#'   \code{err_proc_iid} then \code{err_proc_iid_sigma_scale}. \item If
#'   \code{err_proc_GPP} then \code{err_mult_GPP_sdlog_sigma}.}
#'
#'   * metab_mle: \code{model_name, day_start, day_end, day_tests,
#'   init.GPP.daily, init.Pmax, init.alpha, init.ER.daily, init.ER20,
#'   init.K600.daily}
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
#'   err_obs_sigma, err_obs_phi, err_proc_sigma, err_proc_phi, sim_seed}. Those
#'   arguments whose period-separated name occurs in the default data_daily
#'   argument to metab(sim) can be specified here as NULL, numeric, or a
#'   function to be called each time \code{predict_DO} or \code{predict_metab}
#'   is called on the model. If given as a function, an argument will be called
#'   with any already-evaluated parameters (including the contents of data_daily
#'   and n, the number of dates) passed in as arguments; for example, K600_daily
#'   can see n, discharge.daily, and GPP_daily can see n, discharge.daily, and
#'   K600.daily.
#'
#' @section MLE Initial Values:
#'
#'   For metab_mle models (maximum likelihood estimation), specification
#'   arguments whose names begin with \code{init} are applicable. Which
#'   arguments are required depends on the value of model_name and can be
#'   determined by calling \code{grep('^init.', names(specs(mname)),
#'   value=TRUE)} once for your model name \code{mname} before supplying any
#'   arguments.
#'
#' @param model_name character string identifying the model features. Use
#'   \code{\link{mm_name}} to create a valid name based on desired attributes,
#'   or \code{\link{mm_valid_names}} to see all valid names. Two alternatives to
#'   the names given by \code{mm_valid_names()} are also accepted: (1) a model
#'   type as accepted by the \code{type} argument to \code{mm_name}, which will
#'   be used to create the default model name for that model type, or (2) a full
#'   model file path for custom Bayesian models, as long as basename(model_name)
#'   can still be parsed correctly with \code{mm_parse_name()} and the file
#'   exists. In that case the file may be specified either as a file path
#'   relative to the streamMetabolizer models directory (the first assumption;
#'   this directory can be found with \code{system.file("models",
#'   package="streamMetabolizer")}) or as an absolute path or a path relative to
#'   the current working directory (the second assumption, if the first
#'   assumption turns up no files of the given name).
#' @param engine The software or function to use in fitting the model. Should be
#'   specified via \code{mm_name} rather than here. For \code{type='bayes'},
#'   always \code{'stan'} indicating the software package to use for the MCMC
#'   process (see http://mc-stan.org/). For types in
#'   \code{c('mle','night','sim')} there's again only one option per model (R
#'   functions; these need not be named here but will be noted in the suffix of
#'   the model name, e.g., \code{"m_np_oi_tr_plrckm.nlm"} uses \code{nlm()} for
#'   model fitting). For type='Kmodel', the name of an interpolation or
#'   regression method relating K to the predictor[s] of choice. One of
#'   \code{c("mean", "lm", "loess")}.
#' @inheritParams mm_model_by_ply
#' @inheritParams mm_is_valid_day
#'
#' @param init.GPP.daily the inital value of daily mean GPP (gO2 d^-1 m^-2) to
#'   use in the NLM fitting process. See the MLE Initial Values section under
#'   Details.
#' @param init.Pmax the initial value of Pmax (gO2 d^-1 m^-2) to use in the GPP
#'   versus light relationship in the NLM fitting process. Pmax is the maximum
#'   GPP value of the GPP-light curve. See the MLE Initial Values section under
#'   Details.
#' @param init.alpha the inital value of alpha (gO2 s d^-1 umol^-1, i.e., units
#'   of GPP/light) to use in the GPP versus light relationship in the NLM
#'   fitting process. alpha is the initial slope of the GPP-light curve. See the
#'   MLE Initial Values section under Details.
#' @param init.ER.daily the inital value of daily mean ER (gO2 d^-1 m^-2) to use
#'   in the NLM fitting process. See the MLE Initial Values section under
#'   Details.
#' @param init.ER20 the initial value of ER20 (gO2 d^-1 m^-2) to use in the ER
#'   versus temperature relationship in the NLM fitting process. ER20 is the
#'   respiration rate at 20 degrees C. See the MLE Initial Values section under
#'   Details.
#' @param init.K600.daily the inital value of daily mean K600 (d^-1) to use in
#'   the NLM fitting process. Ignored if K600 is supplied in data_daily, except
#'   for those dates where K600 is NA. If there are any such dates, K600_init
#'   must have a numeric (non-NA) value, as this will be used to estimate K600
#'   for those dates. See the MLE Initial Values section under Details.
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
#'   default is TRUE because these objects often need inspecting.
#' @param keep_mcmc_data FALSE, TRUE, or (for nopool models) a vector of dates
#'   (coerced with as.Date if character, etc.) indicating whether to keep all of
#'   the mcmc model objects (TRUE), none of them (FALSE), or specific dates. The
#'   default is FALSE because these objects can be very large.
#'
#' @param GPP_daily_mu The mean of a dnorm distribution for GPP_daily, the daily
#'   rate of gross primary production
#' @param GPP_daily_lower The lower bound on every fitted value of GPP_daily,
#'   the daily rate of gross primary production. Use values other than -Inf with
#'   caution, recognizing that sometimes the input data are unmodelable and that
#'   a negative estimate of GPP_daily (when unconstrained) could be your only
#'   indication.
#' @param GPP_daily_sigma The standard deviation of a dnorm distribution for
#'   GPP_daily, the daily rate of gross primary production
#' @param alpha_meanlog The mean of a dlnorm (lognormal) distribution for alpha,
#'   the daily initial slope of the Jassby-Platt saturating curve relating GPP
#'   to light
#' @param alpha_sdlog The standard deviation parameter of a dlnorm (lognormal)
#'   distribution for alpha, the daily initial slope of the Jassby-Platt
#'   saturating curve relating GPP to light.
#' @param Pmax_mu The mean of a dnorm (normal) distribution for Pmax, the daily
#'   maximum GPP value of a Jassby-Platt saturating curve relating GPP to light.
#' @param Pmax_sigma The standard deviation of a dnorm (normal) distribution for
#'   Pmax, the daily maximum GPP value of a Jassby-Platt saturating curve
#'   relating GPP to light.
#'
#' @param ER_daily_mu The mean of a dnorm distribution for ER_daily, the daily
#'   rate of ecosystem respiration
#' @param ER_daily_upper The upper (less negative) bound on every fitted value
#'   of ER_daily, the daily rate of ecosystem respiration. Use values other than
#'   Inf with caution, recognizing that sometimes the input data are unmodelable
#'   and that a positive estimate of ER_daily (when unconstrained) could be your
#'   only indication.
#' @param ER_daily_sigma The standard deviation of a dnorm distribution for
#'   ER_daily, the daily rate of ecosystem respiration
#'
#' @param K600_daily_meanlog Applies when pool_K600 is 'none'. The mean of a
#'   dlnorm distribution for K600_daily, the daily rate of reaeration
#' @param K600_daily_sdlog The lognormal scale parameter (standard deviation) of
#'   a dlnorm distribution having meanlog equal to \code{K600_daily_meanlog}
#'   (when pool_K600 is 'none') or \code{K600_daily_predlog} (when pool_K600 is
#'   'normal_sdfixed') for K600_daily, the daily rate of reaeration as corrected
#'   for temperature and the diffusivity of oxygen
#' @param K600_daily_sigma The standard deviation of a dnorm distribution having
#'   mean equal to \code{exp(K600_daily_predlog)} (applicable when pool_K600 is
#'   'linear_sdfixed' or 'binned_sdfixed') for K600_daily, the daily rate of
#'   reaeration as corrected for temperature and the diffusivity of oxygen
#' @param K600_daily_sdlog_sigma hyperparameter for pool_K600 in c('normal').
#'   The scale (= sigma) parameter of a half-normal distribution of sdlog in K ~
#'   lN(meanlog, sdlog), sdlog ~ halfnormal(0, sigma=sdlog_sigma). Visualize the
#'   PDF of K600_daily_sdlog with \code{\link{plot_distribs}}.
#' @param K600_daily_sigma_sigma hyperparameter for pool_K600 in
#'   c('linear','binned'). The scale (= sigma) parameter of a half-normal
#'   distribution of sigma in K ~ lN(meanlog, sigma), sigma ~ halfnormal(0,
#'   sigma=sigma_sigma). Visualize the PDF of K600_daily_sdlog with
#'   \code{\link{plot_distribs}}.
#'
#' @param K600_daily_meanlog_meanlog hyperparameter for pool_K600='normal'. The
#'   mean parameter (meanlog_meanlog) of a lognormal distribution of meanlog in
#'   K ~ lN(meanlog, sdlog), meanlog ~ lN(meanlog_meanlog, meanlog_sdlog)
#' @param K600_daily_meanlog_sdlog hyperparameter for pool_K600='normal'. The
#'   standard deviation parameter (meanlog_sdlog) of a lognormal distribution of
#'   meanlog in K ~ lN(meanlog, sdlog), meanlog ~ lN(meanlog_meanlog,
#'   meanlog_sdlog)
#'
#' @param lnK600_lnQ_intercept_mu hyperparameter for pool_K600 == 'linear'. The
#'   mean of the prior distribution for the intercept parameter in
#'   \code{log(K600) ~ lnK600_lnQ_intercept + lnK600_lnQ_slope*log(Q)}
#' @param lnK600_lnQ_intercept_sigma hyperparameter for pool_K600 == 'linear'.
#'   The standard deviation of the prior distribution for the intercept
#'   parameter in \code{log(K600) ~ lnK600_lnQ_intercept +
#'   lnK600_lnQ_slope*log(Q)}
#' @param lnK600_lnQ_slope_mu hyperparameter for pool_K600='linear'. The mean of
#'   the prior distribution for the slope parameter in \code{log(K600) ~
#'   lnK600_lnQ_intercept + lnK600_lnQ_slope*log(Q)}
#' @param lnK600_lnQ_slope_sigma hyperparameter for pool_K600='linear'. The
#'   standard deviation of the prior distribution for the slope parameter in
#'   \code{log(K600) ~ lnK600_lnQ_intercept + lnK600_lnQ_slope*log(Q)}
#'
#' @param K600_lnQ_nodes_centers data configuration argument for
#'   pool_K600='binned'. numeric vector giving the natural-log-space centers of
#'   the discharge bins. See also \code{\link{calc_bins}}
#' @param K600_lnQ_nodediffs_sdlog hyperparameter for pool_K600='binned'. The
#'   standard deviations of the differences in estimated K600 between successive
#'   lnQ_nodes (bins), where the means of those differences are always zero
#' @param K600_lnQ_nodes_meanlog hyperparameter for pool_K600='binned'. The
#'   means of lognormal prior distributions for the K600_lnQ_nodes parameters.
#' @param K600_lnQ_nodes_sdlog hyperparameter for pool_K600='binned'. The
#'   standard deviations of lognormal prior distributions for the K600_lnQ_nodes
#'   parameters.
#'
#' @param err_obs_iid_sigma_scale The scale (= sigma) parameter of a half-Cauchy
#'   distribution for err_obs_iid_sigma, the standard deviation of the
#'   observation error. Visualize the PDF of err_obs_iid_sigma with
#'   \code{\link{plot_distribs}}.
#' @param err_proc_acor_phi_alpha The alpha (= shape1) parameter on a beta
#'   distribution for err_proc_acor_phi, the autocorrelation coefficient for the
#'   autocorrelated component of process [& sometimes observation] error.
#'   Visualize the PDF of err_proc_acor_phi with \code{\link{plot_distribs}}.
#' @param err_proc_acor_phi_beta The beta (= shape2) parameter on a beta
#'   distribution for err_proc_acor_phi, the autocorrelation coefficient for the
#'   autocorrelated component of process [& sometimes observation] error.
#'   Visualize the PDF of err_proc_acor_phi with \code{\link{plot_distribs}}.
#' @param err_proc_acor_sigma_scale The scale (= sigma) parameter of a
#'   half-Cauchy distribution for err_proc_acor_sigma, the standard deviation of
#'   the autocorrelated component of process [& sometimes observation] error.
#'   Visualize the PDF of err_proc_acor_sigma with \code{\link{plot_distribs}}.
#' @param err_proc_iid_sigma_scale The scale (= sigma) parameter of a
#'   half-Cauchy distribution for err_proc_iid_sigma, the standard deviation of
#'   the uncorrelated (IID) component of process [& sometimes observation]
#'   error. Visualize the PDF of err_proc_iid_sigma with
#'   \code{\link{plot_distribs}}.
#' @param err_mult_GPP_sdlog_sigma The scale parameter of a half-normal
#'   distribution for err_mult_GPP_sdlog, the scale parameter of the lognormal
#'   distribution of err_mult_GPP. err_mult_GPP is multiplied by light and then
#'   normalized to a daily mean of 1 before being multiplied by GPP_daily to
#'   estimate GPP_inst. The effect is a special kind of process error that is
#'   proportional to light (with noise) and is applied to GPP rather than to
#'   dDO/dt.
#'
#' @param params_in Character vector of hyperparameters to pass from the specs
#'   list into the data list for the MCMC run. Will be automatically generated
#'   during the specs() call; need only be revised if you're using a custom
#'   model that requires different hyperparameters.
#'
#' @inheritParams prepdata_bayes
#' @inheritParams runstan_bayes
#'
#' @inheritParams prepdata_Kmodel
#' @inheritParams Kmodel_allply
#'
#' @param K600_lnQ_cnode_meanlog For a sim model with pool_K600='binned'. The
#'   mean of a lognormal distribution describing the y=K600 value of the middle
#'   (or just past middle) node in the piecewise lnK ~ lnQ relationship
#' @param K600_lnQ_cnode_sdlog For a sim model with pool_K600='binned'. The sd
#'   of a lognormal distribution describing the y=K600 value of the middle (or
#'   just past middle) node in the piecewise lnK ~ lnQ relationship
#' @param K600_lnQ_nodediffs_meanlog For a sim model with pool_K600='binned'.
#'   The average (in log space) difference between ln(K) values of successive
#'   nodes. A non-zero value introduces a trend in K ~ Q.
#' @param lnK600_lnQ_nodes For a sim model with pool_K600='binned'. The values
#'   of lnK600 at each node. The default value of this spec is a function that
#'   computes lnK600s based on simulated K~Q relationships.
#'
#' @param discharge_daily Daily values, or a function to generate daily values,
#'   of mean daily discharge in m^3 s^-1. Fixed values may alternatively be
#'   specified as discharge.daily in the data_daily passed to
#'   \code{\link{metab}}.
#' @param DO_mod_1 Daily values, or a function to generate daily values, of the
#'   first DO.mod value on each date. Fixed values may alternatively be
#'   specified as \code{DO.mod.1} in the \code{data_daily} passed to
#'   \code{\link{metab}}. Or may be implied by a \code{DO.obs} column in
#'   \code{data}, from which the first values on each date will be extracted by
#'   \code{metab()}.
#' @param K600_daily Daily values, or a function to generate daily values, of
#'   the reaeration rate constant K600. Fixed values may alternatively be
#'   specified as \code{K600.daily} in the data_daily passed to
#'   \code{\link{metab}}.
#' @param GPP_daily Daily values, or a function to generate daily values, of the
#'   photosynthesis parameter GPP_daily. Fixed values may alternatively be
#'   specified as \code{GPP.daily} in the data_daily passed to
#'   \code{\link{metab}}.
#' @param Pmax Daily values, or a function to generate daily values, of the
#'   photosynthesis parameter Pmax. Fixed values may alternatively be specified
#'   as \code{Pmax} in the data_daily passed to \code{\link{metab}}.
#' @param alpha Daily values, or a function to generate daily values, of the
#'   photosynthesis parameter alpha. Fixed values may alternatively be specified
#'   as \code{alpha} in the data_daily passed to \code{\link{metab}}.
#' @param ER_daily Daily values, or a function to generate daily values, of the
#'   respiration parameter ER_daily. Fixed values may alternatively be specified
#'   as \code{ER.daily} in the data_daily passed to \code{\link{metab}}.
#' @param ER20 Daily values, or a function to generate daily values, of the
#'   respiration parameter ER20. Fixed values may alternatively be specified as
#'   \code{ER20} in the data_daily passed to \code{\link{metab}}.
#'
#' @param err_obs_sigma Daily values, or a function to generate daily values, of
#'   the sd of observation error, or 0 for no observation error. Observation
#'   errors are those applied to DO.mod after generating the full time series of
#'   modeled values.
#' @param err_obs_phi Daily values, or a function to generate daily values, of
#'   the autocorrelation coefficient of the observation errors, or 0 for
#'   uncorrelated errors.
#' @param err_proc_sigma Daily values, or a function to generate daily values,
#'   of the sd of process error, or 0 for no process error. Process errors are
#'   applied at each time step, and therefore propagate into the next timestep.
#' @param err_proc_phi Daily values, or a function to generate daily values, of
#'   the autocorrelation coefficient of the process errors, or 0 for
#'   uncorrelated errors.
#' @param err_round A single value indicating whether simulated DO.obs should be
#'   rounded to simulate the common practice of only reporting a few significant
#'   figures for DO. Use NA for no effect, or an integer as in the \code{digits}
#'   argument to \code{\link{round}} if simulated DO.obs should be rounded to
#'   the given number of digits beyond \code{.}.
#' @param sim_seed NA to specify that each call to predict_DO should generate
#'   new values, or an integer, as in the \code{seed} argument to
#'   \code{\link{set.seed}}, specifying the seed to set before every execution
#'   of predict_DO and/or predict_metab.
#'
#' @return an internally consistent list of arguments that may be passed to
#'   \code{metab} as the \code{specs} argument
#'
#' @importFrom stats rnorm rlnorm
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
  day_tests=c('full_day', 'even_timesteps', 'complete_data', 'pos_discharge', 'pos_depth'),
  required_timestep=NA,
  
  
  ## MLE
  
  # initial values
  init.GPP.daily = 8, 
  init.Pmax = 10,
  init.alpha = 0.0001,
  init.ER.daily = -10, 
  init.ER20 = -10,
  init.K600.daily = 10,
  
  
  ## Bayes
  
  # model setup
  split_dates,
  keep_mcmcs = TRUE,
  keep_mcmc_data = TRUE,
  
  # hyperparameters for non-hierarchical GPP & ER
  GPP_daily_mu = 3.1,
  GPP_daily_lower = -Inf,
  GPP_daily_sigma = 6.0,
  alpha_meanlog = -4.6,
  alpha_sdlog = 0.5,
  Pmax_mu = 10,
  Pmax_sigma = 7,
  ER_daily_mu = -7.1,
  ER_daily_upper = Inf,
  ER_daily_sigma = 7.1,
  
  # hyperparameters for non-hierarchical K600
  K600_daily_meanlog = log(12),
  
  # hyperparameters for hierarchical K600 - normal
  K600_daily_meanlog_meanlog = log(12),
  K600_daily_meanlog_sdlog = 1.32,
  
  # hyperparameters for hierarchical K600 - linear. defaults should be
  # reasonably constrained, not too wide
  lnK600_lnQ_intercept_mu = 2,
  lnK600_lnQ_intercept_sigma = 2.4,
  lnK600_lnQ_slope_mu = 0,
  lnK600_lnQ_slope_sigma = 0.5,
  
  # hyperparameters for hierarchical K600 - binned. K600_daily ~ 
  # lognormal(K600_daily_nodes_meanlog[lnQ_bin], 
  # K600_daily_nodes_sdlog[lnQ_bin]) with linear interpolation among bins before
  # exponentiating. nodes_meanlog and nodes_sdlog may be length b = 
  # length(K600_daily_lnQ_nodes) or length 1 (to be replicated to length b). 
  # -8:6 covers almost all points in Raymond et al. 2012 and will therefore 
  # always be too broad a range for a single stream. -3:3 will catch some
  # streams to rivers as a first cut, though users should still modify
  K600_lnQ_nodes_centers = -3:3, # the x=lnQ values for the nodes
  K600_lnQ_nodediffs_sdlog = 0.5, # for centers 1 apart; for centers 0.2 apart, use 1/5 of this
  K600_lnQ_nodes_meanlog = rep(log(12), length(K600_lnQ_nodes_centers)), # distribs for the y=K600 values of the nodes
  K600_lnQ_nodes_sdlog = rep(1.32, length(K600_lnQ_nodes_centers)),
  
  # hyperparameters for any K pooling or non-pooling strategy
  K600_daily_sdlog = switch(mm_parse_name(model_name)$pool_K600, none=1, normal_sdfixed=0.05, NA),
  K600_daily_sigma = switch(mm_parse_name(model_name)$pool_K600, linear_sdfixed=10, binned_sdfixed=5, NA),
  K600_daily_sdlog_sigma = switch(mm_parse_name(model_name)$pool_K600, normal=0.05, NA),
  K600_daily_sigma_sigma = switch(mm_parse_name(model_name)$pool_K600, linear=1.2, binned=0.24, NA),
  # normal_sdzero, linear_sdzero, and binned_sdzero all have no parameters for this
  
  # hyperparameters for error terms
  err_obs_iid_sigma_scale = 0.03,
  err_proc_iid_sigma_scale = 5,
  err_proc_acor_phi_alpha = 1,
  err_proc_acor_phi_beta = 1,
  err_proc_acor_sigma_scale = 1,
  err_mult_GPP_sdlog_sigma = 1,
  
  # vector of hyperparameters to include as MCMC data
  params_in,
  
  # inheritParams runstan_bayes
  params_out,
  n_chains = 4,
  n_cores = 4,
  burnin_steps = 500,
  saved_steps = 500,
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
  
  # multi-day simulation parameters. already above for bayes:
  # K600_lnQ_nodes_centers, K600_lnQ_nodediffs_sdlog
  K600_lnQ_cnode_meanlog = log(6), # distrib for the y=K600 values of the middle (or just past middle) node
  K600_lnQ_cnode_sdlog = 1, # distrib for the y=K600 values of the middle (or just past middle) node
  K600_lnQ_nodediffs_meanlog = 0.2, # non-zero introduces a trend in K ~ Q
  lnK600_lnQ_nodes = function(K600_lnQ_nodes_centers, K600_lnQ_cnode_meanlog, K600_lnQ_cnode_sdlog,
                              K600_lnQ_nodediffs_meanlog, K600_lnQ_nodediffs_sdlog, ...) {
    sim_Kb(K600_lnQ_nodes_centers, K600_lnQ_cnode_meanlog, K600_lnQ_cnode_sdlog,
           K600_lnQ_nodediffs_meanlog, K600_lnQ_nodediffs_sdlog)
  },
  
  # daily simulation parameters
  discharge_daily = function(n, ...) rnorm(n, 20, 3),
  DO_mod_1 = NULL,
  K600_daily = function(n, K600_daily_predlog=log(10), ...) pmax(0, rnorm(n, K600_daily_predlog, 4)),
  GPP_daily = function(n, ...) pmax(0, rnorm(n, 8, 4)),
  Pmax = function(n, ...) pmax(0, rnorm(n, 10, 2)),
  alpha = function(n, ...) pmax(0, rnorm(n, 0.0001, 0.00002)),
  ER_daily = function(n, ...) pmin(0, rnorm(n, -10, 5)),
  ER20 = function(n, ...) pmin(0, rnorm(n, -10, 4)),
  
  # sub-daily simulation parameters
  err_obs_sigma = 0.01,
  err_obs_phi = 0,
  err_proc_sigma = 0.2,
  err_proc_phi = 0,
  err_round = NA,
  
  # simulation replicability
  sim_seed = NA
  
) {
  
  # make it easier to enter custom specs by creating the type-specific default if model_name %in% 'mle', etc.
  if(model_name %in% eval(formals(mm_name)$type))
    model_name <- mm_name(type=model_name)
  
  # check the validity of the model_name against the list of officially accepted model names
  mm_validate_name(model_name)
  
  # parse the model_name
  features <- mm_parse_name(model_name, expand=TRUE)
  
  # collect info about the arguments
  required <- 'model_name'
  all_possible <- names(formals(specs))
  not_missing <- names(as.list(match.call())[-1]) # the arguments that were given explicitly
  yes_missing <- all_possible[!(all_possible %in% not_missing)]
  prefer_missing <- setdiff(all_possible[sapply(formals(specs), is.symbol)], 'params_out') # the arguments w/o defaults, mostly
  prefer_not_missing <- if(features$type == 'bayes' && features$GPP_fun == 'satlight') {
    c('alpha_meanlog', 'alpha_sdlog', 'Pmax_mu', 'Pmax_sigma') 
  } else {
    c() # could be made more extensive
  }
  
  # argument checks
  if(any(required %in% yes_missing))
    stop("missing and required argument: ", paste(required[required %in% yes_missing], collapse=", "))
  if(any(prefer_not_missing %in% yes_missing)) {
    warn_about <- prefer_not_missing[prefer_not_missing %in% yes_missing]
    warning("you should specify site-appropriate values for all parameters and especially ", paste(warn_about, collapse=", "))
  }
  redundant <- not_missing[not_missing %in% prefer_missing]
  if('engine' %in% redundant) {
    warning("'engine' should be specified in mm_name() rather than specs()")
    redundant <- redundant[redundant != 'engine']
  }
  if(length(redundant) > 0) {
    warning("argument[s] that should usually be specified in revise() rather than specs(): ", paste(redundant, collapse=", "))
  }
  
  # collect the defaults + directly specified arguments
  all_specs <- as.list(environment())
  
  # copy/calculate arguments as appropriate to the model
  specs <- list()
  switch(
    features$type,
    'bayes' = {
      
      # list the specs that will make it all the way to the Stan model as data
      all_specs$params_in <- c(
        switch(
          features$GPP_fun,
          linlight=c('GPP_daily_mu','GPP_daily_lower','GPP_daily_sigma'),
          satlight=c('alpha_meanlog', 'alpha_sdlog', 'Pmax_mu', 'Pmax_sigma')),
        c('ER_daily_mu','ER_daily_upper','ER_daily_sigma'),
        switch(
          features$pool_K600_type,
          none=c('K600_daily_meanlog'),
          normal=c('K600_daily_meanlog_meanlog', 'K600_daily_meanlog_sdlog'),
          linear=c('lnK600_lnQ_intercept_mu', 'lnK600_lnQ_intercept_sigma', 'lnK600_lnQ_slope_mu', 'lnK600_lnQ_slope_sigma'),
          binned=c('K600_lnQ_nodediffs_sdlog', 'K600_lnQ_nodes_meanlog', 'K600_lnQ_nodes_sdlog')),
        switch(
          features$pool_K600_sd,
          zero=c(),
          fixed=switch(features$pool_K600_type, none=, normal='K600_daily_sdlog', linear=, binned='K600_daily_sigma'),
          fitted=switch(features$pool_K600_type, normal='K600_daily_sdlog_sigma', linear=, binned='K600_daily_sigma_sigma')
        ),
        if(features$err_obs_iid) 'err_obs_iid_sigma_scale',
        if(features$err_proc_acor) c('err_proc_acor_phi_alpha', 'err_proc_acor_phi_beta', 'err_proc_acor_sigma_scale'),
        if(features$err_proc_iid) 'err_proc_iid_sigma_scale',
        if(features$err_proc_GPP) 'err_mult_GPP_sdlog_sigma'
      )
      
      # list all needed arguments
      included <- c(
        # model setup
        'model_name', 'engine', 'split_dates', 'keep_mcmcs', 'keep_mcmc_data',
        
        # date ply day_tests
        'day_start', 'day_end', 'day_tests', 'required_timestep',
        
        # discharge binning parameters are not params_in, though they're 
        # conceptually related and therefore colocated in formals(specs)
        if(features$pool_K600_type == 'binned') c('K600_lnQ_nodes_centers'),
        
        # params_in is both a vector of specs to include and a vector to include in specs
        all_specs$params_in, 'params_in',
        
        # inheritParams runstan_bayes
        'params_out', 'n_chains', 'n_cores', 
        'burnin_steps', 'saved_steps', 'thin_steps', 'verbose'
      )
      
      # compute some arguments
      if('engine' %in% yes_missing) {
        all_specs$engine <- features$engine
      }
      if('split_dates' %in% yes_missing) {
        all_specs$split_dates <- switch(
          features$pool_K600_type,
          'none' = FALSE, # pretty sure FALSE is faster. also allows hierarchical error terms
          'normal'=, 'linear'=, 'binned' = FALSE, 
          stop("unknown pool_K600; unsure how to set split_dates"))
      }
      if(features$pool_K600_type == 'binned') {
        # defaults are for linear pool_K600 & need adjustment for binned method
        all_specs$K600_daily_beta_mu <- rep(10, length(all_specs$K600_daily_lnQ_nodes))
        all_specs$K600_daily_beta_sigma <- rep(10, length(all_specs$K600_daily_lnQ_nodes))
      }
      if('params_out' %in% yes_missing) {
        all_specs$params_out <- c(
          c('GPP', 'ER', 'DO_R2'),
          switch(
            features$GPP_fun,
            linlight=c('GPP_daily'),
            satlight=c('alpha', 'Pmax')),
          c('ER_daily', 'K600_daily'),
          switch(
            features$pool_K600_type,
            none=c(),
            normal=c('K600_daily_predlog'),
            linear=c('K600_daily_predlog', 'lnK600_lnQ_intercept', 'lnK600_lnQ_slope'),
            binned=c('K600_daily_predlog', 'lnK600_lnQ_nodes')), 
          if(features$pool_K600_sd == 'fitted')
            switch(
              features$pool_K600_type,
              normal='K600_daily_sdlog',
              linear=, binned='K600_daily_sigma'),
          if(features$err_obs_iid) c('err_obs_iid_sigma', 'err_obs_iid'),
          if(features$err_proc_acor) c('err_proc_acor', 'err_proc_acor_phi', 'err_proc_acor_sigma'),
          if(features$err_proc_iid) c('err_proc_iid_sigma', 'err_proc_iid'),
          if(features$err_proc_GPP) c('err_proc_GPP', 'GPP_pseudo_R2'))
      }
      
      # check for errors/inconsistencies
      model_path <- tryCatch(
        mm_locate_filename(model_name), 
        error=function(e) {
          warning(e)
          return(model_name)
        })
      if(features$engine == "NA") 
        stop('engine must be specified for Bayesian models')
      
    },
    'mle' = {
      # determine which init values will be needed
      . <- '.dplyr.var'
      init.needs <- paste0('init.', get_param_names(model_name)$required)
      
      # list all needed arguments
      included <- c('model_name', 'day_start', 'day_end', 'day_tests', 'required_timestep', init.needs)

    }, 
    'night' = {
      # list all needed arguments
      included <- c('model_name', 'day_start', 'day_end', 'day_tests', 'required_timestep')
      
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
        'model_name', 'engine', 'day_start', 'day_end', 'day_tests', 'required_timestep',
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
          if('filters' %in% yes_missing) all_specs['filters'] <- list(c()) # need special syntax to assign c(). see https://stackoverflow.com/a/7945259/3203184
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
          if('other_args' %in% yes_missing) all_specs$other_args <- list(possible_args=names(formals('loess'))[-which(names(formals('loess')) %in% c('formula','data','weights'))])
        }
      )
      
    },
    'sim' = {
      # determine which daily parameters will be needed
      par_needs <- gsub('\\.', '_', unlist(get_param_names(model_name)[c('optional','required')]))
      
      # list all needed arguments
      included <- c(
        'model_name', 'day_start', 'day_end', 'day_tests', 'required_timestep',
        switch(
          features$pool_K600,
          none=c(),
          normal=stop("pool_K600='normal' unavailable for now; try 'binned' instead"),
          linear=stop("pool_K600='linear' unavailable for now; try 'binned' instead"), # 'discharge_daily', etc.
          binned=c('K600_lnQ_nodes_centers', 
                   'K600_lnQ_cnode_meanlog', 'K600_lnQ_cnode_sdlog', 'K600_lnQ_nodediffs_meanlog', 'K600_lnQ_nodediffs_sdlog',
                   'lnK600_lnQ_nodes')),
        par_needs, 'err_round', 'sim_seed')
      
      if(features$pool_K600 == 'binned') {
        if('K600_lnQ_nodes_centers' %in% yes_missing) # override the default, which is for 'bayes' rather than 'sim'
          all_specs$K600_lnQ_nodes_centers <- function(discharge.daily, ...) calc_bins(log(discharge.daily), 'width', width=0.2)$bounds
      }
    }
  )
  
  # stop if truly irrelevant arguments were given
  if(length(irrelevant <- not_missing[!(not_missing %in% included)]) > 0) 
    stop("irrelevant argument: ", paste(irrelevant, collapse=", "))
  
  # return just the arguments we actually need
  add_specs_class(all_specs[included])
  
}
