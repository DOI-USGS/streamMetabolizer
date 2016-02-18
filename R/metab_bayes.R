#' @include metab_model-class.R
NULL

#' Basic Bayesian metabolism model fitting function
#' 
#' Fits a Bayesian model to estimate GPP and ER from input data on DO, 
#' temperature, light, etc. See \code{\link{mm_name}} to choose a Bayesian model
#' and \code{\link{specs}} for relevant options for the \code{model_specs}
#' argument.
#' 
#' @author Alison Appling, Bob Hall
#' @inheritParams metab_model_prototype
#' @inheritParams mm_is_valid_day
#' @return A metab_bayes object containing the fitted model.
#' @examples
#' \dontrun{
#' # set the date in several formats
#' start.chron <- chron::chron(dates="08/23/12", times="22:00:00")
#' end.chron <- chron::chron(dates="08/25/12", times="06:00:00")
#' start.posix <- as.POSIXct(format(start.chron, "%Y-%m-%d %H:%M:%S"), tz="GMT")
#' end.posix <- as.POSIXct(format(end.chron, "%Y-%m-%d %H:%M:%S"), tz="GMT")
#' mid.date <- as.Date(start.posix + (end.posix - start.posix)/2)
#' start.numeric <- as.numeric(start.posix - as.POSIXct(format(mid.date, "%Y-%m-%d 00:00:00"),
#'    tz="GMT"), units='hours')
#' end.numeric <- as.numeric(end.posix - as.POSIXct(format(mid.date, "%Y-%m-%d 00:00:00"),
#'   tz="GMT"), units='hours')
#' 
#' # get, format, & subset data
#' vfrench <- streamMetabolizer:::load_french_creek(attach.units=FALSE)
#' vfrenchshort <- vfrench[vfrench$solar.time >= start.posix & vfrench$solar.time <= end.posix, ]
#' 
#' # fit
#' mm <- metab_bayes(data=vfrenchshort, day_start=start.numeric, 
#'   day_end=end.numeric)
#' get_fit(mm)[2,c('GPP_daily_median','ER_daily_median','K600_daily_median')]
#' streamMetabolizer:::load_french_creek_std_mle(vfrenchshort, estimate='PRK')
#' get_fitting_time(mm)
#' plot_DO_preds(predict_DO(mm))
#' }
#' @export
#' @family metab_model
metab_bayes <- function(
  data=mm_data(solar.time, DO.obs, DO.sat, depth, temp.water, light), data_daily=mm_data(NULL), # inheritParams metab_model_prototype
  model_specs=specs('b_np_oi_pm_km.jags'), # inheritParams metab_model_prototype
  info=NULL, day_start=4, day_end=27.99, # inheritParams metab_model_prototype
  tests=c('full_day', 'even_timesteps', 'complete_data') # inheritParams mm_is_valid_day
) {
  
  fitting_time <- system.time({
    # Check data for correct column names & units
    dat_list <- mm_validate_data(data, if(missing(data_daily)) NULL else data_daily, "metab_bayes")
    data <- dat_list[['data']]
    data_daily <- dat_list[['data_daily']]
    
    # Check and parse model file path. First try the streamMetabolizer models
    # dir, then try a regular path, then complain / continue depending on
    # whether we found a file. Add the complete path to model_specs. This is
    # best done here, in the metab_bayes call, so that models can be defined,
    # passed to another computer, and still run successfully.
    model_specs$model_path <- system.file(paste0("models/", model_specs$model_name), package="streamMetabolizer")
    if(!file.exists(model_specs$model_path)) 
      model_specs$model_path <- model_specs$model_name
    if(!file.exists(model_specs$model_path)) 
      stop("could not locate the model file at ", model_specs$model_path)
    
    # check the format of keep_mcmcs (more checks, below, are bayes_fun-specific)
    if(is.logical(model_specs$keep_mcmcs) && length(model_specs$keep_mcmcs) != 1) {
      stop("if keep_mcmcs is logical, it must have length 1")
    }
    
    # model the data
    switch(
      model_specs$bayes_fun,
      'bayes_1ply' = {
        if(!is.logical(model_specs$keep_mcmcs)) {
          model_specs$keep_mcmcs <- as.Date(model_specs$keep_mcmcs)
        }
        # one day at a time, splitting into overlapping 31.5-hr 'plys' for each date
        bayes_all <- mm_model_by_ply(
          bayes_1ply, data=data, data_daily=data_daily, # for mm_model_by_ply
          day_start=day_start, day_end=day_end, # for mm_model_by_ply and mm_is_valid_day
          tests=tests, # for mm_is_valid_day
          model_specs=model_specs) # for bayes_1ply
        if('mcmcfit' %in% names(bayes_all)) {
          bayes_mcmc <- bayes_all$mcmcfit
          names(bayes_mcmc) <- bayes_all$solar.date
          bayes_all$mcmcfit <- NULL
        } else {
          bayes_mcmc <- NULL
        }
      },
      'bayes_all' = {
        # all days at a time, after first filtering out bad days
        if((day_end - day_start) > 24) warning("multi-day models should probably have day_end - day_start <= 24 hours")
        filtered <- mm_filter_valid_days(data, data_daily, day_start=day_start, day_end=day_end, tests=tests)
        bayes_all <- bayes_allply(
          data_all=filtered$data, data_daily_all=filtered$data_daily,
          model_specs=model_specs)
      }, {
        stop("unrecognized bayes_fun")
      }
    )
  })
  
  # Package and return results
  mm <- metab_model(
    model_class="metab_bayes", 
    info=info,
    fit=bayes_all,
    mcmc=bayes_mcmc,
    fitting_time=fitting_time,
    args=list(
      model_specs=model_specs,
      day_start=day_start, day_end=day_end, tests=tests
    ),
    data=data,
    data_daily=data_daily)
  
  # Update data with DO predictions
  mm@data <- predict_DO(mm)
  
  # Return
  mm
}


#### helpers ####

#' Make daily metabolism estimates from input parameters
#' 
#' Called from metab_bayes().
#' 
#' @inheritParams mm_model_by_ply_prototype
#' @inheritParams mm_is_valid_day
#' @inheritParams metab_model_prototype
#' @return data.frame of estimates and MCMC model diagnostics
#' @importFrom stats setNames
#' @keywords internal
bayes_1ply <- function(
  data_ply, data_daily_ply, day_start, day_end, solar_date, # inheritParams mm_model_by_ply_prototype
  tests=c('full_day', 'even_timesteps', 'complete_data'), # inheritParams mm_is_valid_day
  model_specs # inheritParams metab_model_prototype
) {
  
  # Provide ability to skip a poorly-formatted day for calculating 
  # metabolism, without breaking the whole loop. Just collect 
  # problems/errors as a list of strings and proceed. Also collect warnings.
  validity <- mm_is_valid_day(data_ply, day_start=day_start, day_end=day_end, tests=tests)
  stop_strs <- if(isTRUE(validity)) character(0) else validity
  warn_strs <- character(0)
  
  # Calculate metabolism by Bayesian MCMC
  if(length(stop_strs) == 0) {
    bayes_1day <- withCallingHandlers(
      tryCatch({
        # first: try to run the bayes fitting function
        data_list <- prepdata_bayes(
          data=data_ply, data_daily=data_daily_ply, solar_date=solar_date,
          model_specs=model_specs, priors=model_specs$priors)
        model_specs$keep_mcmc <- if(is.logical(model_specs$keep_mcmcs)) {
          isTRUE(model_specs$keep_mcmcs)
        } else {
          isTRUE(solar_date %in% model_specs$keep_mcmcs)
        }
        all_mcmc_args <- c('bayes_software','model_path','params_out','keep_mcmc','n_chains','n_cores','adapt_steps','burnin_steps','saved_steps','thin_steps','verbose')
        do.call(mcmc_bayes, c(
          list(data_list=data_list),
          model_specs[all_mcmc_args[all_mcmc_args %in% names(model_specs)]]))
      }, error=function(err) {
        # on error: give up, remembering error. dummy values provided below
        stop_strs <<- c(stop_strs, err$message)
        NA
      }), warning=function(war) {
        # on warning: record the warning and run again
        warn_strs <<- c(warn_strs, war$message)
        invokeRestart("muffleWarning")
      })
  } 

  # stop_strs may have accumulated during nlm() call. If failed, use dummy data 
  # to fill in the model output with NAs.
  if(length(stop_strs) > 0) {
    bayes_1day <- data.frame(GPP_daily_mean=NA)
  }
  
  # Return, reporting any results, warnings, and errors
  data.frame(bayes_1day,
             warnings=paste0(warn_strs, collapse="; "), 
             errors=paste0(stop_strs, collapse="; "),
             stringsAsFactors=FALSE)
}


#' Make daily metabolism estimates from input parameters using a hierarchical
#' approach.
#' 
#' Called from metab_bayes().
#' 
#' @param data_all data.frame of the form \code{mm_data(solar.time, DO.obs, 
#'   DO.sat, depth, temp.water, light)} and containing data for just one 
#'   estimation-day (this may be >24 hours but only yields estimates for one 
#'   24-hour period)
#' @param data_daily_all data.frame of daily priors, if appropriate to the given
#'   model_path
#' @inheritParams metab_model_prototype
#' @return data.frame of estimates and MCMC model 
#'   diagnostics
#' @keywords internal
bayes_allply <- function(
  data_all, data_daily_all, 
  model_specs # inheritParams metab_model_prototype
) {
  stop("this function is not yet implemented")
}


#### helpers to the helper ####

#' Prepare data for passing to JAGS or Stan
#' 
#' This function accepts exactly one day's worth of data, (one ply, which might 
#' be 24 hrs or 31.5 or so), which should already be validated. It prepares the 
#' data needed to run a Bayesian MCMC method to estimate GPP, ER, and K600.
#' 
#' @inheritParams metab_model_prototype
#' @inheritParams mm_model_by_ply_prototype
#' @inheritParams metab_bayes
#' @param priors logical. Should the data list be modified such that JAGS will 
#'   return priors rather than posteriors?
#' @return list of data for input to runjags_bayes or runstan_bayes
#' @keywords internal
prepdata_bayes <- function(
  data, data_daily, solar_date, # inheritParams metab_model_prototype
  model_specs, # inheritParams metab_bayes
  priors=FALSE
) {
  
  # Useful info for setting the MCMC data
  timestep_days <- suppressWarnings(mean(as.numeric(diff(v(data$solar.time)), units="days"), na.rm=TRUE))
  has_oierr <- grepl('_oi(pc)*(pi)*[[:punct:]]', basename(model_specs$model_path))
  has_pierr <- grepl('_(oi)*(pc)*pi[[:punct:]]', basename(model_specs$model_path))
  has_pcerr <- grepl('_(oi)*pc(pi)*[[:punct:]]', basename(model_specs$model_path))
  
  # Format the data for JAGS/Stan. Stan disallows period-separated names, so
  # change all the input data to underscore-separated. parameters given in
  # model_specs are already underscore-separated for this reason
  data_list = c(
    list(
      
      # Daily
      n = nrow(data),
      DO_obs_1 = data$DO.obs[1],
      
      # Every timestep
      frac_GPP = data$light/sum(data$light[as.character(data$solar.time,"%Y-%m-%d")==as.character(solar_date)]),
      frac_ER = rep(timestep_days, nrow(data)),
      frac_D = rep(timestep_days, nrow(data)), # the yackulic shortcut models rely on this being even over time
      KO2_conv = convert_k600_to_kGAS(k600=1, temperature=data$temp.water, gas="O2"),
      depth = data$depth,
      DO_sat = data$DO.sat,
      DO_obs = data$DO.obs
    ),
    
    model_specs[c(
      # Hyperparameters
      c('GPP_daily_mu','GPP_daily_sigma','ER_daily_mu','ER_daily_sigma','K600_daily_mu','K600_daily_sigma'), # metabolism
      if(has_oierr) c('err_obs_iid_sigma_min','err_obs_iid_sigma_max') else c(), # uncorrelated observation error
      if(has_pcerr) c('err_proc_acor_phi_min','err_proc_acor_phi_max','err_proc_acor_sigma_min','err_proc_acor_sigma_max') else c(), # autocorrelated process error
      if(has_pierr) c('err_proc_iid_sigma_min','err_proc_iid_sigma_max') else c() # uncorrelated component of process error
    )]
  )
  if(priors) {
    data_list <- data_list[-which(names(data_list)=="DO_obs")]
  }
  
  data_list
}

#' Run an MCMC simulation on a formatted data ply
#' 
#' @param data_list a formatted list of inputs to the JAGS model
#' @param bayes_software character string indicating which software to use
#' @param model_path the JAGS model file to use, as a full file path
#' @param params_out a character vector of parameters whose values in the MCMC 
#'   runs should be recorded and summarized
#' @param keep_mcmc logical. If TRUE, the Jags or Stan output object will be 
#'   saved. Be careful; these can be big, and a run with many models might 
#'   overwhelm R's memory.
#' @param n_chains the number of chains to run
#' @param n_cores the number of cores to apply to this run
#' @param adapt_steps the number of steps per chain to use in adapting the model
#' @param burnin_steps the number of steps per chain to run and ignore before
#'   starting to collect MCMC 'data'
#' @param saved_steps the number of MCMC steps per chain to save
#' @param thin_steps the number of steps to move before saving another step. 1
#'   means save all steps.
#' @param verbose logical. give status messages?
#' @return a data.frame of outputs
#' @import parallel
#' @keywords internal
mcmc_bayes <- function(data_list, bayes_software=c('stan','jags'), model_path, params_out, keep_mcmc=FALSE, n_chains=4, n_cores=4, adapt_steps=1000, burnin_steps=4000, saved_steps=40000, thin_steps=1, verbose=FALSE) {
  bayes_software <- match.arg(bayes_software)
  bayes_function <- switch(bayes_software, jags = runjags_bayes, stan = runstan_bayes)
  
  tot_cores <- detectCores()
  if (!is.finite(tot_cores)) { tot_cores <- 1 } 
  n_cores <- min(tot_cores, n_cores)
  message(paste0("MCMC (",bayes_software,"): requesting ",n_chains," chains on ",n_cores," of ",tot_cores," available cores\n"))
  
  bayes_function(
    data_list=data_list, model_path=model_path, params_out=params_out, keep_mcmc=keep_mcmc, n_chains=n_chains, n_cores=n_cores, 
    adapt_steps=adapt_steps, burnin_steps=burnin_steps, saved_steps=saved_steps, thin_steps=thin_steps, verbose=verbose)
}

#' Run JAGS on a formatted data ply
#' 
#' Seems to need to import rjags but does not, for now, because I can't get 
#' rjags to install on the Condor cluster. Including an import rjags line here 
#' allowed runjags to do its job last time I tried.
#' 
#' @inheritParams mcmc_bayes
#' @param ... args passed to other runxx_bayes functions but ignored here
#' @import dplyr
#' @keywords internal
runjags_bayes <- function(data_list, model_path, params_out, keep_mcmc=FALSE, n_chains=4, adapt_steps=1000, burnin_steps=4000, saved_steps=40000, thin_steps=1, verbose=FALSE, ...) {
  
  if(!requireNamespace("runjags", quietly = TRUE)) {
    stop("the runjags package is required for JAGS MCMC models")
  }
  
  inits_fun <- function(chain) {
    list(.RNG.name=
           c("base::Wichmann-Hill",
             "base::Marsaglia-Multicarry",
             "base::Super-Duper",
             "base::Mersenne-Twister")[chain])
    # Let JAGS initialize other parameters automatically
  }
  
  runjags_out <- runjags::run.jags(
    method=c("rjags","parallel","snow")[2],
    model=model_path,
    monitor=params_out,
    data=data_list,
    inits=inits_fun,
    n.chains=n_chains,
    adapt=adapt_steps,
    burnin=burnin_steps,
    sample=saved_steps,
    thin=thin_steps,
    summarise=TRUE,
    plots=FALSE,
    silent.jags=!verbose)
  
  # format output into a 1-row data.frame
  jags_mat <- cbind(runjags_out$summary$statistics[,c('Naive SE','Time-series SE')], 
                    runjags_out$summaries,
                    runjags_out$summary$quantiles) %>% as.matrix() # combine 2 matrices of statistics
  names_params <- rep(rownames(jags_mat), each=ncol(jags_mat)) # the GPP, ER, etc. part of the name
  names_stats <- rep(tolower(gsub(" |-", "_", gsub("%", "pct", colnames(jags_mat)))), times=nrow(jags_mat)) # add the mean, sd, etc. part of the name
  jags_out <- jags_mat %>% t %>% c %>% # get a 1D vector of GPP_daily_mean, GPP_sd, ..., ER_daily_mean, ER_daily_sd, ... etc
    t %>% as.data.frame() %>% # convert from 1D vector to 1-row data.frame
    setNames(paste0(names_params, "_", names_stats))
  
  # add the model object if requested
  if(keep_mcmc == TRUE) {
    jags_out <- mutate(jags_out, mcmcfit=list(runjags_out))
  }

  return(jags_out)
}

#' Run Stan on a formatted data ply
#' 
#' @inheritParams mcmc_bayes
#' @param ... args passed to other runxx_bayes functions but ignored here
#' @import parallel
#' @import dplyr
#' @keywords internal
runstan_bayes <- function(data_list, model_path, params_out, keep_mcmc=FALSE, n_chains=4, n_cores=4, burnin_steps=1000, saved_steps=1000, thin_steps=1, verbose=FALSE, ...) {
  
  # stan() can't find its own function cpp_object_initializer() unless the 
  # namespace is loaded. requireNamespace is somehow not doing this. Thoughts
  # (not solution):
  # https://stat.ethz.ch/pipermail/r-devel/2014-September/069803.html
  if(!suppressPackageStartupMessages(require(rstan))) {
    stop("the rstan package is required for Stan MCMC models")
  }
  
  runstan_out <- rstan::stan(
    file=model_path,
    data=data_list,
    pars=params_out,
    include=TRUE,
    chains=n_chains,
    warmup=burnin_steps,
    iter=saved_steps+burnin_steps,
    thin=thin_steps,
    init="random",
    save_dso=TRUE, # must be true if you're using more than one core
    verbose=verbose,
    open_progress=FALSE,
    cores=n_cores)

  # this is a good place for a breakpoint when running small numbers of models
  # manually (or keep_mcmc also helps with inspection)
  #   show(runstan_out)
  #   rstan::plot(runstan_out)
  #   pairs(runstan_out)
  #   traceplot(runstan_out)
  
  # format output into a 1-row data.frame. see ls('package:rstan')
  stan_mat <- rstan::summary(runstan_out)$summary
  names_params <- rep(rownames(stan_mat), each=ncol(stan_mat)) # the GPP, ER, etc. part of the name
  names_stats <- rep(gsub("%", "pct", colnames(stan_mat)), times=nrow(stan_mat)) # add the mean, sd, etc. part of the name
  stan_out <- stan_mat %>% t %>% c %>% # get a 1D vector of GPP_daily_mean, GPP_sd, ..., ER_daily_mean, ER_daily_sd, ... etc
    t %>% as.data.frame() %>% # convert from 1D vector to 1-row data.frame
    setNames(paste0(names_params, "_", names_stats)) 
  
  # add the model object if requested
  if(keep_mcmc == TRUE) {
    stan_out <- mutate(stan_out, mcmcfit=list(runstan_out))
  }
  
  return(stan_out)
}


#### metab_bayes class ####

#' Metabolism model fitted by Bayesian MCMC
#' 
#' \code{metab_bayes} models use Bayesian MCMC methods to fit values of GPP, ER,
#' and K for a given DO curve.
#' 
#' @exportClass metab_bayes
#' @family metab.model.classes
setClass(
  "metab_bayes", 
  contains="metab_model",
  slots=c(mcmc="ANY")
)

#' Extract any MCMC model objects that were stored with the model
#' 
#' A function specific to metab_bayes models. Returns an MCMC object or, for 
#' nopool models, a list of MCMC objects. These objects are not saved by 
#' default; see \code{keep_mcmcs} argument to \code{\link{specs}} for options.
#' 
#' @param metab_model A Bayesian metabolism model (metab_bayes) from which to 
#'   return the MCMC model object[s]
#' @return The MCMC model object[s]
#' @export
get_mcmc <- function(metab_model) {
  UseMethod("get_mcmc")
}

#' Retrieve any MCMC model object[s] that were saved with a metab_bayes model
#' 
#' @inheritParams get_mcmc
#' @export 
#' @family get_mcmc
get_mcmc.metab_bayes <- function(metab_model) {
  metab_model@mcmc
}