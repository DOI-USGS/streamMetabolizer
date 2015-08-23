#' @include metab_model-class.R
NULL

#' Basic Bayesian metabolism model fitting function
#' 
#' Fits a model to estimate GPP and ER from input data on DO, 
#'   temperature, light, etc.
#'   
#' @author Alison Appling, Bob Hall
#' @inheritParams metab_model_prototype
#' @inheritParams mm_is_valid_day
#' @inheritParams runjags_bayes
#' @return A metab_bayes object containing the fitted model.
#' @examples
#' \dontrun{
#'  metab_bayes(data=data.frame(empty="shouldbreak"))
#' }
#' @export
#' @family metab_model
metab_bayes <- function(
  data=mm_data(local.time, DO.obs, DO.sat, depth, temp.water, light), data_daily=NULL, info=NULL, day_start=-1.5, day_end=30, # inheritParams metab_model_prototype
  tests=c('full_day', 'even_timesteps', 'complete_data'), # inheritParams mm_is_valid_day
  model_file='metab_bayes_simple.txt', max_cores=4, adapt_steps=10, burnin_steps=40, num_saved_steps=400, thin_steps=1 # inheritParams runjags_bayes
) {
  
  # Check data for correct column names & units
  data <- mm_validate_data(data, "metab_bayes")
  
  # model the data
  if(model_file %in% c('metab_bayes_simple.txt','metab_bayes_procerr.txt')) {
    # one day at a time, splitting into overlapping 31.5-hr 'plys' for each date
    bayes_all <- mm_model_by_ply(
      bayes_1ply, data=data, data_daily=data_daily, # for mm_model_by_ply
      day_start=day_start, day_end=day_end, # for mm_model_by_ply and mm_is_valid_day
      tests=tests, # for mm_is_valid_day
      model_file=model_file, max_cores=max_cores, adapt_steps=adapt_steps, burnin_steps=burnin_steps, num_saved_steps=num_saved_steps, thin_steps=thin_steps) # for bayes_1ply
      
  }
  
  # Package and return results
  metab_model(
    model_class="metab_bayes", 
    info=info,
    fit=bayes_all,
    args=list(day_start=day_start, day_end=day_end, tests=tests,
              model_file=model_file, max_cores=max_cores, adapt_steps=adapt_steps, burnin_steps=burnin_steps, num_saved_steps=num_saved_steps, thin_steps=thin_steps,
              calc_DO_fun=calc_DO_mod # metab_bayes always uses the equivalent of calc_DO_mod. specify here for predict_DO
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
#' @inheritParams runjags_bayes
#' @return data.frame of estimates and \code{\link[stats]{nlm}} model 
#'   diagnostics
#' @keywords internal
bayes_1ply <- function(
  data_ply, data_daily_ply, day_start=-1.5, day_end=30, local_date, # inheritParams mm_model_by_ply_prototype
  tests=c('full_day', 'even_timesteps', 'complete_data'), # inheritParams mm_is_valid_day
  model_file, max_cores=4, adapt_steps=10, burnin_steps=40, num_saved_steps=400, thin_steps=1 # inheritParams runjags_bayes
) {
  
  # Provide ability to skip a poorly-formatted day for calculating 
  # metabolism, without breaking the whole loop. Just collect 
  # problems/errors as a list of strings and proceed. Also collect warnings.
  validity <- mm_is_valid_day(data_ply, day_start=day_start, day_end=day_end, 
                              tests=c("full_day","even_timesteps","complete_data"), 
                              need_complete=c("DO.obs","DO.sat","depth","temp.water","light"))
  stop_strs <- if(isTRUE(validity)) character(0) else validity
  warn_strs <- character(0)
  
  # Calculate metabolism by Bayesian MCMC
  if(length(stop_strs) == 0) {
    bayes.1d <- withCallingHandlers(
      tryCatch({
        # first: try to run the bayes fitting function
        data_list <- prepjags_bayes(data_ply, model_file)
        runjags_bayes(data_list=data_list, model_file, max_cores=max_cores, adapt_steps=adapt_steps, 
                      burnin_steps=burnin_steps, num_saved_steps=num_saved_steps, thin_steps=thin_steps)
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
    bayes.1d <- setNames(as.list(rep(NA, 6)), 
                         c('mean.GPP.daily', 'mean.ER.daily', 'mean.K600.daily', 
                           'sd.GPP.daily', 'sd.ER.daily', 'sd.K600.daily'))
  }
  
  # Return, reporting any results, warnings, and errors
  data.frame(GPP=bayes.1d$mean.GPP.daily, GPP.sd=bayes.1d$sd.GPP.daily, 
             ER=bayes.1d$mean.ER.daily, ER.sd=bayes.1d$sd.ER.daily, 
             K600=bayes.1d$mean.K600.daily, K600.sd=bayes.1d$sd.K600.daily,
             warnings=paste0(warn_strs, collapse="; "), 
             errors=paste0(stop_strs, collapse="; "),
             stringsAsFactors=FALSE)
}

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
    bayes.1d <- setNames(as.list(rep(NA, 6)), 
                         c('mean.GPP.daily', 'mean.ER.daily', 'mean.K600.daily', 
                           'sd.GPP.daily', 'sd.ER.daily', 'sd.K600.daily'))
  }
  
  # Return, reporting any results, warnings, and errors
  data.frame(GPP=bayes.1d$mean.GPP.daily, GPP.sd=bayes.1d$sd.GPP.daily, 
             ER=bayes.1d$mean.ER.daily, ER.sd=bayes.1d$sd.ER.daily, 
             K600=bayes.1d$mean.K600.daily, K600.sd=bayes.1d$sd.K600.daily,
             warnings=paste0(warn_strs, collapse="; "), 
             errors=paste0(stop_strs, collapse="; "),
             stringsAsFactors=FALSE)
}


#### helpers to the helper ####

#' Prepare data for passing to JAGS
#' 
#' This function accepts exactly one day's worth of data, (one ply, which might 
#' be 24 hrs or 31.5 or so), which should already be validated. It prepares the 
#' data needed to run a Bayesian MCMC method to estimate GPP, ER, and K600.
#' 
#' @param data_ply one day's worth of data
#' @param priors logical. Should the data list be modified such that JAGS will
#'   return priors rather than posteriors?
#' @return list of data for input to runjags
#' @export
prepjags_bayes <- function(data_ply, priors=FALSE) {
  
  #Useful info for setting JAGS data
  local.date <- names(which.max(table(as.Date(data_ply$local.time))))
  timestep.days <- suppressWarnings(mean(as.numeric(diff(v(data_ply$local.time)), units="days"), na.rm=TRUE))
  
  # Format the data for Jags
  data_list = list(
    
    # Daily
    n = nrow(data_ply),
    
    # Every timestep
    frac.GPP = data_ply$light/sum(data_ply$light[strftime(data_ply$local.time,"%Y-%m-%d")==local.date]),
    frac.ER = rep(timestep.days, nrow(data_ply)),
    frac.D = rep(timestep.days, nrow(data_ply)),
    KO2.conv = convert_k600_to_kGAS(k600=1, temperature=data_ply$temp.water, gas="O2"),
    depth = data_ply$depth,
    DO.sat = data_ply$DO.sat,
    DO.obs = data_ply$DO.obs,
    
    # Constants
    GPP.daily.mu = 10,
    GPP.daily.tau = 1/(10^2),
    ER.daily.mu = -10,
    ER.daily.tau = 1/(10^2),
    K600.daily.mu = 10,
    K600.daily.tau = 1/(10^2)
  )
  if(priors) {
    data_list <- data_list[-which(names(data_list)=="DO.obs")]
  }
  data_list <- data_list[sapply(data_list, function(data_elem) !is.null(data_elem))]
  
  data_list
}

#' Actually run JAGS on a formatted data ply
#' 
#' Seems to need to import rjags but does not, for now, because I can't get
#' rjags to install on the Condor cluster. Including an import rjags line here
#' allowed runjags to do its job last time I tried.
#' 
#' @param data_list a formatted list of inputs to the JAGS model
#' @param model_file the JAGS model file to use
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
runjags_bayes <- function(data_list, model_file='metab_bayes_simple.txt', max_cores=4, adapt_steps=1000, burnin_steps=4000, num_saved_steps=40000, thin_steps=1) {
  
  inits_fun <- function(chain) {
    list(.RNG.name=
           c("base::Wichmann-Hill",
             "base::Marsaglia-Multicarry",
             "base::Super-Duper",
             "base::Mersenne-Twister")[chain])
    # Let JAGS initialize other parameters automatically
  }
  
  parameters = switch(
    model_file,
    'metab_bayes_simple.txt'=c("GPP.daily", "ER.daily", "K600.daily", "err.obs.sigma"),
    'metab_bayes_procerr.txt'=c("GPP.daily", "ER.daily", "K600.daily", "err.obs.sigma", "err.proc.sigma")
  )
  
  n_cores = detectCores()
  if (!is.finite(n_cores)) { n_cores = 1 } 
  n_chains = max(3, min(max_cores , max(1, n_cores)))
  message(paste0("Found ",n_cores," cores; requesting ",n_chains," chains.\n"))
  
  runjags_out <- run.jags(
    method=c("rjags","parallel","snow")[2],
    model=system.file(paste0("models/", model_file), package="streamMetabolizer"),
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

