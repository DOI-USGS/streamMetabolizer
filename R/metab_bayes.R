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
#' start.posix <- as.POSIXct(format(start.chron, "%Y-%m-%d %H:%M:%S"), tz="UTC")
#' end.posix <- as.POSIXct(format(end.chron, "%Y-%m-%d %H:%M:%S"), tz="UTC")
#' mid.date <- as.Date(start.posix + (end.posix - start.posix)/2)
#' start.numeric <- as.numeric(start.posix - as.POSIXct(format(mid.date, "%Y-%m-%d 00:00:00"),
#'    tz="UTC"), units='hours')
#' end.numeric <- as.numeric(end.posix - as.POSIXct(format(mid.date, "%Y-%m-%d 00:00:00"),
#'   tz="UTC"), units='hours')
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
  model_specs=specs(mm_name('bayes')), # inheritParams metab_model_prototype
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
    
    # check the format of keep_mcmcs (more checks, below, are split_dates-specific)
    if(is.logical(model_specs$keep_mcmcs)) {
      if(length(model_specs$keep_mcmcs) != 1) {
        stop("if keep_mcmcs is logical, it must have length 1")
      }
    } else if(model_specs$split_dates == FALSE) {
      stop("if split_dates==FALSE, keep_mcmcs must be a single logical value")
    }
    
        # all days at a time, after first filtering out bad days
        if((day_end - day_start) > 24) warning("multi-day models should probably have day_end - day_start <= 24 hours")
        filtered <- mm_filter_valid_days(data, data_daily, day_start=day_start, day_end=day_end, tests=tests)
        bayes_all <- bayes_allply(
          data_all=filtered$data, data_daily_all=filtered$data_daily,
          model_specs=model_specs)
    # model the data. create outputs bayes_all (a data.frame) and bayes_mcmc (an MCMC object from JAGS or Stan)
    if(model_specs$split_dates == TRUE) {
      if(!is.logical(model_specs$keep_mcmcs)) {
        model_specs$keep_mcmcs <- as.Date(model_specs$keep_mcmcs)
      }
      # one day at a time, splitting into overlapping 31.5-hr 'plys' for each date
      bayes_daily <- mm_model_by_ply(
        bayes_1ply, data=data, data_daily=data_daily, # for mm_model_by_ply
        day_start=day_start, day_end=day_end, # for mm_model_by_ply and mm_is_valid_day
        tests=tests, # for mm_is_valid_day
        model_specs=model_specs) # for bayes_1ply
      # if we saved the modeling object[s] in the df, pull them out now
      if('mcmcfit' %in% names(bayes_daily)) {
        bayes_mcmc <- bayes_daily$mcmcfit
        names(bayes_mcmc) <- bayes_daily$date
        bayes_daily$mcmcfit <- NULL
      } else {
        bayes_mcmc <- NULL
      }
      bayes_all <- list(daily=bayes_daily)
      
    } else if(model_specs$split_dates == FALSE) {
    }
  
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
  data_ply, data_daily_ply, day_start, day_end, ply_date, # inheritParams mm_model_by_ply_prototype
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
          data=data_ply, data_daily=data_daily_ply, ply_date=ply_date,
          model_specs=model_specs, priors=model_specs$priors)
        model_specs$keep_mcmc <- if(is.logical(model_specs$keep_mcmcs)) {
          isTRUE(model_specs$keep_mcmcs)
        } else {
          isTRUE(ply_date %in% model_specs$keep_mcmcs)
        }
        all_mcmc_args <- c('engine','model_path','params_out','split_dates','keep_mcmc','n_chains','n_cores','adapt_steps','burnin_steps','saved_steps','thin_steps','verbose')
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

  # stop_strs may have accumulated during prepdata_bayes() or mcmc_bayes()
  # calls. If failed, use dummy data to fill in the model output with NAs.
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
#' @importFrom unitted v
#' @keywords internal
prepdata_bayes <- function(
  data, data_daily, ply_date=NA, # inheritParams metab_model_prototype
  model_specs, # inheritParams metab_bayes
  priors=FALSE
) {
  
  # remove units if present
  data <- v(data)
  data_daily <- v(data_daily)
  if(length(ply_date) != 1) stop("ply_date must have length 1")
  if(!is.na(ply_date)) data$date <- ply_date
  
  # define a function to package 1+ days of obs of a variable into a time x date matrix
  date_table <- table(data$date)
  num_dates <- length(date_table)
  num_daily_obs <- unique(unname(date_table))
  if(length(num_daily_obs) > 1) stop("dates have differing numbers of rows; observations cannot be combined in matrix")
  time_by_date_matrix <- function(vec) {
    matrix(data=vec, nrow=num_daily_obs, ncol=num_dates, byrow=FALSE)
  }
  
  # double-check that our dates are going to line up with the input dates. this 
  # should be redundant w/ above date_table checks, so just being extra careful
  obs_dates <- time_by_date_matrix(as.character(data$date, "%Y-%m-%d"))
  unique_dates <- apply(obs_dates, MARGIN=2, FUN=function(col) unique(col))
  if(!all.equal(unique_dates, names(date_table))) stop("couldn't fit given dates into matrix")
  
  # confirm that every day has the same modal timestep and put a value on that timestep
  obs_times <- time_by_date_matrix(as.numeric(data$solar.time - data$solar.time[1], units='days'))
  unique_timesteps <- unique(apply(obs_times, MARGIN=2, FUN=function(col) unique(round(diff(col), digits=12)))) # 10 digits is 8/1000000 of a second. 14 digits exceeds machine precision for datetimes
  if(length(unique_timesteps) != 1) stop("could not determine a single timestep for all observations")
  timestep_days <- mean(apply(obs_times, MARGIN=2, FUN=function(col) mean(diff(col))))
  
  # parse model name into features for deciding what data to include
  features <- mm_parse_name(basename(model_specs$model_path))
  
  # Format the data for JAGS/Stan. Stan disallows period-separated names, so
  # change all the input data to underscore-separated. parameters given in
  # model_specs are already underscore-separated for this reason
  data_list = c(
    list(
      
      # Overall
      d = num_dates,
      
      # Daily
      n = num_daily_obs, # one value applicable to every day
      DO_obs_1 = array(time_by_date_matrix(data$DO.obs)[1,], dim=num_dates), # duplication of effort below should be small compared to MCMC time
      
      # Every timestep
      frac_GPP = {
        # normalize light by the sum of light in the first 24 hours of the time window
        in_solar_day <- apply(obs_times, MARGIN=2, FUN=function(col) {col - col[1] <= 1} )
        mat_light <- time_by_date_matrix(data$light)
        sweep(mat_light, MARGIN=2, STATS=colSums(mat_light*in_solar_day), FUN=`/`)
      },
      frac_ER  = time_by_date_matrix(timestep_days),
      frac_D   = time_by_date_matrix(timestep_days), # the yackulic shortcut models rely on this being constant over time
      KO2_conv = time_by_date_matrix(convert_k600_to_kGAS(k600=1, temperature=data$temp.water, gas="O2")),
      depth    = time_by_date_matrix(data$depth),
      DO_sat   = time_by_date_matrix(data$DO.sat),
      DO_obs   = time_by_date_matrix(data$DO.obs)
    ),
    
    model_specs[c(
      # Hyperparameters - this section should be identical to the
      # hyperparameters section of specs|bayes
      c('GPP_daily_mu','GPP_daily_sigma','ER_daily_mu','ER_daily_sigma'),
      switch(
        features$pool_K600,
        none=c('K600_daily_mu', 'K600_daily_sigma'),
        normal=c('K600_daily_mu_mu', 'K600_daily_mu_sigma', 'K600_daily_sigma_shape', 'K600_daily_sigma_rate'),
        linear=c('K600_daily_beta0_mu', 'K600_daily_beta0_sigma', 'K600_daily_beta1_mu', 'K600_daily_beta1_sigma', 'K600_daily_sigma_shape', 'K600_daily_sigma_rate'),
        binned=stop('need to think about this one')),
      if(features$err_obs_iid) c('err_obs_iid_sigma_min', 'err_obs_iid_sigma_max'),
      if(features$err_proc_acor) c('err_proc_acor_phi_min', 'err_proc_acor_phi_max', 'err_proc_acor_sigma_min', 'err_proc_acor_sigma_max'),
      if(features$err_proc_iid) c('err_proc_iid_sigma_min', 'err_proc_iid_sigma_max')
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
#' @param engine character string indicating which software to use
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
mcmc_bayes <- function(data_list, engine=c('stan','jags'), model_path, params_out, split_dates, keep_mcmc=FALSE, n_chains=4, n_cores=4, adapt_steps=1000, burnin_steps=4000, saved_steps=40000, thin_steps=1, verbose=FALSE) {
  engine <- match.arg(engine)
  bayes_function <- switch(engine, jags = runjags_bayes, stan = runstan_bayes)
  
  tot_cores <- detectCores()
  if (!is.finite(tot_cores)) { tot_cores <- 1 } 
  n_cores <- min(tot_cores, n_cores)
  message(paste0("MCMC (",engine,"): requesting ",n_chains," chains on ",n_cores," of ",tot_cores," available cores\n"))
  
  bayes_function(
    data_list=data_list, model_path=model_path, params_out=params_out, split_dates=split_dates, keep_mcmc=keep_mcmc, n_chains=n_chains, n_cores=n_cores, 
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
runjags_bayes <- function(data_list, model_path, params_out, split_dates, keep_mcmc=FALSE, n_chains=4, adapt_steps=1000, burnin_steps=4000, saved_steps=40000, thin_steps=1, verbose=FALSE, ...) {
  
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
runstan_bayes <- function(data_list, model_path, params_out, split_dates, keep_mcmc=FALSE, n_chains=4, n_cores=4, burnin_steps=1000, saved_steps=1000, thin_steps=1, verbose=FALSE, ...) {
  
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