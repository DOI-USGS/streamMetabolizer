#' @include metab_model-class.R
NULL

#' Basic Bayesian metabolism model fitting function
#' 
#' Fits a model to estimate GPP and ER from input data on DO, 
#'   temperature, light, etc.
#'   
#' @author Alison Appling, Bob Hall
#' @param data data.frame with columns having the same names, units, and format 
#'   as the default. See \code{\link{mm_data}} for a full data description.
#' @param info Any metadata you would like to package within the metabolism 
#'   model.
#' @inheritParams runjags_bayes
#' @inheritParams mm_is_valid_day
#' @return A metab_bayes object containing the fitted model.
#' @examples
#' \dontrun{
#'  metab_bayes(data=data.frame(empty="shouldbreak"))
#' }
#' @export
#' @family metab_model
metab_bayes <- function(
  data=mm_data(local.time, DO.obs, DO.sat, depth, temp.water, light), # args for bayes_1ply
  info=NULL, # args for new("metab_bayes")
  #fit_distribs=c("GPP","ER","K600"), fit_autocor=c("GPP","ER","K600"), fit_resids=c("K600"), fit_function_K600=c("splineKvT","lmKvQ"),
  maxCores=4, adaptSteps=10, burnInSteps=40, numSavedSteps=400, thinSteps=1, # args for runjags_bayes
  tests=c('full_day', 'even_timesteps', 'complete_data'), day_start=-1.5, day_end=30 # args for mm_is_valid_day, mm_model_by_ply
) {
  
  # Check data for correct column names & units
  data <- mm_validate_data(data, "metab_bayes")
  
  # model the data, splitting into overlapping 31.5-hr 'plys' for each date
  bayes.all <- mm_model_by_ply(
    data, bayes_1ply, # for mm_model_by_ply
    day_start=day_start, day_end=day_end, # for mm_model_by_ply and mm_is_valid_day
    tests=tests, # for mm_is_valid_day
    maxCores=maxCores, adaptSteps=adaptSteps, burnInSteps=burnInSteps, numSavedSteps=numSavedSteps, thinSteps=thinSteps) # for bayes_1ply
  
  # Package and return results
  new("metab_bayes", 
      info=info,
      fit=bayes.all,
      # metab_bayes always uses the equivalent of calc_DO_mod
      args=list(calc_DO_fun=calc_DO_mod, day_start=day_start, day_end=day_end,
                maxCores=maxCores, adaptSteps=adaptSteps, burnInSteps=burnInSteps, 
                numSavedSteps=numSavedSteps, thinSteps=thinSteps),
      data=data,
      pkg_version=as.character(packageVersion("streamMetabolizer")))
}


#### helpers ####

#' Make daily metabolism estimates from input parameters
#' 
#' Called from metab_bayes().
#' 
#' @param data_ply data.frame of the form \code{mm_data(local.time, DO.obs, DO.sat,
#'   depth, temp.water, light)} and containing data for just one estimation-day
#'   (this may be >24 hours but only yields estimates for one 24-hour period)
#' @inheritParams runjags_bayes
#' @param ... additional args passed from mm_model_by_ply and ignored here
#' @return data.frame of estimates and \code{\link[stats]{nlm}} model 
#'   diagnostics
#' @keywords internal
bayes_1ply <- function(data_ply, maxCores=4, adaptSteps=1000, burnInSteps=4000, numSavedSteps=40000, thinSteps=1, ...) {
  
  # Provide ability to skip a poorly-formatted day for calculating 
  # metabolism, without breaking the whole loop. Just collect 
  # problems/errors as a list of strings and proceed. Also collect warnings.
  stop_strs <- mm_is_valid_day(data_ply, need_complete=c("DO.obs","DO.sat","depth","temp.water","light"))
  warn_strs <- character(0)
  
  # Calculate metabolism by Bayesian MCMC
  if(length(stop_strs) == 0) {
    bayes.1d <- withCallingHandlers(
      tryCatch({
        # first: try to run the bayes fitting function
        data.list <- prepjags_bayes(data_ply)
        runjags_bayes(dataList=data.list, maxCores=maxCores, adaptSteps=adaptSteps, 
                             burnInSteps=burnInSteps, numSavedSteps=numSavedSteps, thinSteps=thinSteps)
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
  date <- names(which.max(table(as.Date(data_ply$local.time))))
  timestep.days <- suppressWarnings(mean(as.numeric(diff(v(data_ply$local.time)), units="days"), na.rm=TRUE))
  
  # Format the data for Jags
  dataList = list(
    
    # Daily
    n = nrow(data_ply),
    
    # Every timestep
    frac.GPP = data_ply$light/sum(data_ply$light[strftime(data_ply$local.time,"%Y-%m-%d")==date]),
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
    dataList <- dataList[-which(names(dataList)=="DO.obs")]
  }
  
  dataList
}

#' Actually run JAGS on a formatted data ply
#' 
#' Seems to need to import rjags but does not, for now, because I can't get
#' rjags to install on the Condor cluster. Including an import rjags line here
#' allowed runjags to do its job last time I tried.
#' 
#' @param dataList a formatted list of inputs to the JAGS model
#' @param maxCores the maximum number of cores to apply to this run
#' @param adaptSteps the number of steps to use in adapting the model
#' @param burnInSteps the number of steps to run and ignore before starting to 
#'   collect MCMC 'data'
#' @param numSavedSteps the number of MCMC steps to save
#' @param thinSteps the number of steps to move before saving another step
#' @return a data.frame of outputs
#' @import runjags
#' @import parallel
#' @export
runjags_bayes <- function(dataList, maxCores=4, adaptSteps=1000, burnInSteps=4000, numSavedSteps=40000, thinSteps=1) {
  
  initsFun <- function(chain) {
    list(.RNG.name=
           c("base::Wichmann-Hill",
             "base::Marsaglia-Multicarry",
             "base::Super-Duper",
             "base::Mersenne-Twister")[chain])
    # Let JAGS initialize other parameters automatically
  }
  
  parameters = c("GPP.daily", "ER.daily", "K600.daily", "DO.err.tau")
  
  nCores = detectCores()
  if (!is.finite(nCores)) { nCores = 1 } 
  nChains = max(3, min(maxCores , max(1, nCores)))
  message(paste0("Found ",nCores," cores; requesting ",nChains," chains.\n"))
  
  runJagsOut <- run.jags(
    method=c("rjags","parallel","snow")[2],
    model=system.file("extdata/metab_bayes_jags.txt", package="streamMetabolizer"),
    monitor=parameters,
    data=dataList,
    inits=initsFun,
    n.chains=nChains,
    adapt=adaptSteps,
    burnin=burnInSteps,
    sample=ceiling(numSavedSteps/nChains),
    thin=thinSteps,
    summarise=TRUE,
    plots=FALSE)
  
  # format output into a data.frame
  jagsOut <- as.data.frame(c(
    mean=as.list(runJagsOut$summary$statistics[,'Mean']), 
    sd=as.list(runJagsOut$summary$statistics[,'SD'])))
  
  return(jagsOut)
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

