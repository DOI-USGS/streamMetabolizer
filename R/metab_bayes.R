#' @include metab_model-class.R
NULL

#' Basic Bayesian metabolism model fitting function
#' 
#' Fits a model to estimate GPP and ER from input data on DO, temperature, 
#' light, etc.
#' 
#' Current and future models: models/bayes/jags/nopool_obserr.txt 
#' models/bayes/jags/nopool_procobserr.txt models/bayes/jags/KfQ_procobserr.txt 
#' models/bayes/jags/PRpoolKfQ_procobserr.txt 
#' models/bayes/jags/PRKpool_procobserr.txt
#' 
#' @author Alison Appling, Bob Hall
#' @inheritParams metab_model_prototype
#' @inheritParams mm_is_valid_day
#' @param model_specs a list of model specifications and parameters for a model.
#'   Although this may be specified manually, it is easier to use a predefined
#'   function from the \code{\link{specs_bayes}} family with a name beginning
#'   with "specs_bayes". The help file for that function lists the necessary
#'   parameters, describes them in detail, and gives default values.
#' @inheritParams runjags_bayes
#' @return A metab_bayes object containing the fitted model.
#' @examples
#' \dontrun{
#'  metab_bayes(data=data.frame(empty="shouldbreak"))
#' }
#' @export
#' @family metab_model
metab_bayes <- function(
  data=mm_data(local.time, DO.obs, DO.sat, depth, temp.water, light), data_daily=mm_data(NULL), info=NULL, day_start=-1.5, day_end=30, # inheritParams metab_model_prototype
  tests=c('full_day', 'even_timesteps', 'complete_data'), # inheritParams mm_is_valid_day
  model_specs=specs_bayes_jags_nopool_obserr() # model_specs in inheritParams runjags_bayes
) {
  
  # Check data for correct column names & units
  dat_list <- mm_validate_data(data, if(missing(data_daily)) NULL else data_daily, "metab_bayes")
  data <- dat_list[['data']]
  data_daily <- dat_list[['data_daily']]
  
  # Check and parse model file path. First try the streamMetabolizer models dir,
  # then try a regular path, then give up / continue depending on whether we
  # found a file. Add the complete path to model_specs
  model_specs$model_path <- system.file(paste0("models/bayes/", model_specs$model_file), package="streamMetabolizer")
  if(!file.exists(model_specs$model_path)) 
    model_specs$model_path <- model_specs$model_file 
  if(!file.exists(model_specs$model_path)) 
    stop(suppressWarnings(paste0(
      "could not locate the model file at either\n",
      normalizePath(file.path(system.file("models", package="streamMetabolizer"), model_specs$model_file)), " or\n",
      normalizePath(model_specs$model_file))))
  
  # model the data
  switch(
    model_specs$bayes_fun,
    'bayes_1ply' = {
      # one day at a time, splitting into overlapping 31.5-hr 'plys' for each date
      bayes_all <- mm_model_by_ply(
        bayes_1ply, data=data, data_daily=data_daily, # for mm_model_by_ply
        day_start=day_start, day_end=day_end, # for mm_model_by_ply and mm_is_valid_day
        tests=tests, # for mm_is_valid_day
        model_specs=model_specs) # for bayes_1ply
    },
    'bayes_all' = {
      # all days at a time, after first filtering out bad days
      filtered <- mm_filter_valid_days(data, data_daily, day_start=6, day_end=30, tests=tests)
      bayes_all <- bayes_allply(
        data_all=filtered$data, data_daily_all=filtered$data_daily,
        model_specs=model_specs)
    }, {
      stop("unrecognized bayes_fun")
    }
  )
  
  # Package and return results
  metab_model(
    model_class="metab_bayes", 
    info=info,
    fit=bayes_all,
    args=list(
      day_start=day_start, day_end=day_end, tests=tests, 
      calc_DO_fun=calc_DO_mod, # metab_bayes always uses the equivalent of calc_DO_mod. specify here for predict_DO
      model_specs=model_specs
    ),
    data=data,
    data_daily=data_daily)
}


#### helpers ####

#' Make daily metabolism estimates from input parameters
#' 
#' Called from metab_bayes().
#' 
#' @inheritParams mm_model_by_ply_prototype
#' @inheritParams mm_is_valid_day
#' @inheritParams metab_bayes
#' @return data.frame of estimates and MCMC model 
#'   diagnostics
#' @keywords internal
bayes_1ply <- function(
  data_ply, data_daily_ply, day_start=-1.5, day_end=30, local_date, # inheritParams mm_model_by_ply_prototype
  tests=c('full_day', 'even_timesteps', 'complete_data'), # inheritParams mm_is_valid_day
  model_specs # inheritParams metab_bayes
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
        data_list <- prepjags_bayes(
          data=data_ply, data_daily=data_daily_ply, local_date=local_date,
          model_specs=model_specs, priors=model_specs$priors)
        do.call(runjags_bayes, c(
          list(data_list=data_list), 
          model_specs[c('model_path','max_cores','adapt_steps','burnin_steps','num_saved_steps','thin_steps')]))
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
    bayes_1day <- setNames(as.list(rep(NA, 6)), 
                           c('mean.GPP.daily', 'mean.ER.daily', 'mean.K600.daily', 
                             'sd.GPP.daily', 'sd.ER.daily', 'sd.K600.daily'))
  }
  
  # Return, reporting any results, warnings, and errors
  data.frame(GPP=bayes_1day$mean.GPP.daily, GPP.sd=bayes_1day$sd.GPP.daily, 
             ER=bayes_1day$mean.ER.daily, ER.sd=bayes_1day$sd.ER.daily, 
             K600=bayes_1day$mean.K600.daily, K600.sd=bayes_1day$sd.K600.daily,
             warnings=paste0(warn_strs, collapse="; "), 
             errors=paste0(stop_strs, collapse="; "),
             stringsAsFactors=FALSE)
}


#' Make daily metabolism estimates from input parameters using a hierarchical
#' approach.
#' 
#' Called from metab_bayes().
#' 
#' @param data_all data.frame of the form \code{mm_data(local.time, DO.obs, 
#'   DO.sat, depth, temp.water, light)} and containing data for just one 
#'   estimation-day (this may be >24 hours but only yields estimates for one 
#'   24-hour period)
#' @param data_daily_all data.frame of daily priors, if appropriate to the given
#'   model_path
#' @inheritParams metab_bayes
#' @return data.frame of estimates and MCMC model 
#'   diagnostics
#' @keywords internal
bayes_allply <- function(
  data_all, data_daily_all, 
  model_specs # inheritParams metab_bayes
) {
  stop("this function is not yet implemented")
}


#### helpers to the helper ####

#' Prepare data for passing to JAGS
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
#' @return list of data for input to runjags
#' @export
prepjags_bayes <- function(
  data, data_daily, local_date, # inheritParams metab_model_prototype
  model_specs, # inheritParams metab_bayes
  priors=FALSE
) {
  
  # Useful info for setting JAGS data
  timestep_days <- suppressWarnings(mean(as.numeric(diff(v(data$local.time)), units="days"), na.rm=TRUE))
  has_procerr <- grepl('_proc[[:alpha:]]*err[[:punct:]]', basename(model_specs$model_path))
  
  # Format the data for JAGS
  data_list = c(
    list(
      
      # Daily
      n = nrow(data),
      
      # Every timestep
      frac.GPP = data$light/sum(data$light[strftime(data$local.time,"%Y-%m-%d")==local_date]),
      frac.ER = rep(timestep_days, nrow(data)),
      frac.D = rep(timestep_days, nrow(data)),
      KO2.conv = convert_k600_to_kGAS(k600=1, temperature=data$temp.water, gas="O2"),
      depth = data$depth,
      DO.sat = data$DO.sat,
      DO.obs = data$DO.obs
    ),
    
    model_specs[c(
      # Hyperparameters
      c('GPP.daily.mu','GPP.daily.sigma','ER.daily.mu','ER.daily.sigma','K600.daily.mu','K600.daily.sigma'), # metabolism
      if(has_procerr) c('err.proc.phi.min','err.proc.phi.max','err.proc.sigma.min','err.proc.sigma.max') else c(), # process error
      c('err.obs.sigma.min','err.obs.sigma.max') # observation error
    )]
  )
  if(priors) {
    data_list <- data_list[-which(names(data_list)=="DO.obs")]
  }
  
  data_list
}

#' Actually run JAGS on a formatted data ply
#' 
#' Seems to need to import rjags but does not, for now, because I can't get
#' rjags to install on the Condor cluster. Including an import rjags line here
#' allowed runjags to do its job last time I tried.
#' 
#' @param data_list a formatted list of inputs to the JAGS model
#' @param model_path the JAGS model file to use, as a full file path
#' @param max_cores the maximum number of cores to apply to this run
#' @param adapt_steps the number of steps to use in adapting the model
#' @param burnin_steps the number of steps to run and ignore before starting to 
#'   collect MCMC 'data'
#' @param num_saved_steps the number of MCMC steps to save
#' @param thin_steps the number of steps to move before saving another step
#' @return a data.frame of outputs
#' @import runjags
#' @import parallel
#' @export
runjags_bayes <- function(data_list, model_path, max_cores=4, adapt_steps=1000, burnin_steps=4000, num_saved_steps=40000, thin_steps=1) {
  
  inits_fun <- function(chain) {
    list(.RNG.name=
           c("base::Wichmann-Hill",
             "base::Marsaglia-Multicarry",
             "base::Super-Duper",
             "base::Mersenne-Twister")[chain])
    # Let JAGS initialize other parameters automatically
  }
  
  parameters = switch(
    basename(model_path),
    'nopool_obserr.txt'=c("GPP.daily", "ER.daily", "K600.daily", "err.obs.sigma"),
    'nopool_procobserr.txt'=c("GPP.daily", "ER.daily", "K600.daily", "err.obs.sigma", "err.proc.sigma")
  )
  
  n_cores = detectCores()
  if (!is.finite(n_cores)) { n_cores = 1 } 
  n_chains = max(3, min(max_cores , max(1, n_cores)))
  message(paste0("Found ",n_cores," cores; requesting ",n_chains," chains.\n"))
  
  runjags_out <- run.jags(
    method=c("rjags","parallel","snow")[2],
    model=model_path,
    monitor=parameters,
    data=data_list,
    inits=inits_fun,
    n.chains=n_chains,
    adapt=adapt_steps,
    burnin=burnin_steps,
    sample=ceiling(num_saved_steps/n_chains),
    thin=thin_steps,
    summarise=TRUE,
    plots=FALSE)
  
  # format output into a data.frame
  jags_out <- as.data.frame(c(
    mean=as.list(runjags_out$summary$statistics[,'Mean']), 
    sd=as.list(runjags_out$summary$statistics[,'SD'])))
  
  return(jags_out)
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
  contains="metab_model"
)

