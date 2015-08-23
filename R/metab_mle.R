#' @include metab_model-class.R
NULL

#' Maximum likelihood metabolism model fitting function
#' 
#' Fits a model to estimate GPP and ER from input data on DO, temperature, 
#' light, etc.
#' 
#' @param data data.frame with columns having the same names, units, and format 
#'   as the default. See \code{\link{mm_data}} for a full data description.
#' @param info Any metadata you would like to package within the metabolism 
#'   model.
#' @inheritParams mm_is_valid_day
#' @param data_daily Optional. Leave missing if K600 is to be fitted, or give a
#'   three-column data.frame with columns as in default if K600 is to be taken
#'   as given
#' @param calc_DO_fun the function to use to build DO estimates from GPP, ER, 
#'   etc. default is calc_DO_mod, but could also be calc_DO_mod_by_diff
#' @return A metab_mle object containing the fitted model.
#'   
#' @import dplyr
#' @author Alison Appling, Jordan Read; modeled on LakeMetabolizer
#' @examples
#' \dontrun{
#'  metab_mle(data=data.frame(empty="shouldbreak"))
#' }
#' @export
#' @family metab_model
metab_mle <- function(
  data=mm_data(local.time, DO.obs, DO.sat, depth, temp.water, light), data_daily=mm_data(local.date, K600, optional='all'), calc_DO_fun=calc_DO_mod, # args for mle_1ply
  info=NULL, # args for new('metab_mle')
  tests=c('full_day', 'even_timesteps', 'complete_data'), day_start=-1.5, day_end=30 # args for mm_is_valid_day, mm_model_by_ply
) {
  
  # Check data for correct column names & units
  data <- mm_validate_data(data, "metab_mle")
  
  # model the data, splitting into overlapping 31.5-hr 'plys' for each date
  mle_all <- mm_model_by_ply(
    mle_1ply, data, if(!missing(daily.data)) daily.data else NULL, # for mm_model_by_ply
    day_start=day_start, day_end=day_end, # for mm_model_by_ply and mm_is_valid_day
    tests=tests, # for mm_is_valid_day
    K600=K600, calc_DO_fun=calc_DO_fun) # for mle_1ply
  
  # Package and return results
  new("metab_mle", 
      info=info,
      fit=mle_all,
      args=list(K600=K600, calc_DO_fun=calc_DO_fun, tests=tests, day_start=day_start, day_end=day_end), # keep in order passed to function
      data=c(list(data=data), if(!missing(data_daily)) list(data_daily=data_daily)),
      pkg_version=as.character(packageVersion("streamMetabolizer")))
}


#### helpers ####

#' Make daily metabolism estimates from input parameters
#' 
#' Called from metab_mle().
#' 
#' @param data_ply data.frame of the form \code{mm_data(local.time, DO.obs, DO.sat,
#'   depth, temp.water, light)} and containing data for just one estimation-day
#'   (this may be >24 hours but only yields estimates for one 24-hour period)
#' @param calc_DO_fun the function to use to build DO estimates from GPP, ER,
#'   etc. default is calc_DO_mod, but could also be calc_DO_mod_by_diff
#' @param ... additional args passed from mm_model_by_ply and ignored here
#' @return data.frame of estimates and \code{\link[stats]{nlm}} model 
#'   diagnostics
#' @keywords internal
mle_1ply <- function(data_ply, K600=NULL, calc_DO_fun=calc_DO_mod, 
                     tests=c('full_day', 'even_timesteps', 'complete_data'), day_start=-1.5, day_end=30, date, ...) {
  
  # Provide ability to skip a poorly-formatted day for calculating 
  # metabolism, without breaking the whole loop. Just collect 
  # problems/errors as a list of strings and proceed. Also collect warnings.
  timestep.days <- suppressWarnings(mean(as.numeric(diff(v(data_ply$local.time)), units="days"), na.rm=TRUE))
  validity <- mm_is_valid_day(
    data_ply, # data split by mm_model_by_ply
    tests=tests, day_start=day_start, day_end=day_end, # args passed from metab_mle
    timestep_days=timestep.days, need_complete=c("DO.obs","DO.sat","depth","temp.water","light")) # args supplied here
  stop_strs <- if(isTRUE(validity)) character(0) else validity
  warn_strs <- character(0)

  # Calculate metabolism by non linear minimization of an MLE function
  if(length(stop_strs) == 0) {
    nlm.args <- c(
      list(
        negloglik_1ply,
        hessian=TRUE,
        p = c(GPP=3, ER=-5, K600=5)[if(is.null(K600)) 1:3 else 1:2]
      ),
      if(!is.null(K600)) {
        list(K600=K600[match(date, as.character(K600$date)),"K600"])
      } else {
        list()
      },
      as.list(
        data_ply[c("DO.obs","DO.sat","depth","temp.water")]
      ),
      list(
        frac.GPP = data_ply$light/sum(data_ply[strftime(data_ply$local.time,"%Y-%m-%d")==as.character(date),'light']),
        frac.ER = timestep.days,
        frac.D = timestep.days,
        calc_DO_fun = calc_DO_fun
      ))
    
    mle.1d <- withCallingHandlers(
      tryCatch({
        # first: try to run the MLE fitting function
        mle.1d <- do.call(nlm, nlm.args)
        # if we were successful, also compute the confidence interval as in 
        # http://www.stat.umn.edu/geyer/5931/mle/mle.pdf section 2.3, which says: 
        # 'Inverse Fisher information gives the asymptotic variance matrix of the 
        # MLE. From it, we can construct asymptotic confidence intervals.' See also 
        # http://stats.stackexchange.com/questions/27033/in-r-given-an-output-from-optim-with-a-hessian-matrix-how-to-calculate-paramet
        inv_fish <- solve(mle.1d$hessian) 
        mle.1d$sd <- sqrt(diag(inv_fish))
        mle.1d
      }, error=function(err) {
        # on error: give up, remembering error. dummy values provided below
        stop_strs <<- c(stop_strs, err$message)
        NA
      }), warning=function(war) {
        # on warning: record the warning and run nlm again
        warn_strs <<- c(warn_strs, war$message)
        invokeRestart("muffleWarning")
      })
    
  } 
  
  # Return, reporting any results, warnings, and errors. if the model fitting
  # failed, use dummy data to fill in the output with NAs.
  if(length(stop_strs) > 0) {
    data.frame(
      GPP=NA, ER=NA, K600=NA,
      GPP.sd=NA, ER.sd=NA, K600.sd=NA,
      GPP.grad=NA, ER.grad=NA, K600.grad=NA,
      minimum=NA, code=NA, iterations=NA, 
      warnings=paste0(warn_strs, collapse="; "), 
      errors=paste0(stop_strs, collapse="; "),
      stringsAsFactors=FALSE)
  } else {
    data.frame(
      GPP=mle.1d$estimate[1], ER=mle.1d$estimate[2], K600=if(is.null(K600)) mle.1d$estimate[3] else nlm.args$K600,
      GPP.sd=mle.1d$sd[1], ER.sd=mle.1d$sd[2], K600.sd=if(is.null(K600)) mle.1d$sd[3] else NA,
      GPP.grad=mle.1d$gradient[1], ER.grad=mle.1d$gradient[2], K600.grad=if(is.null(K600)) mle.1d$gradient[3] else NA,
      minimum=mle.1d$minimum, code=mle.1d$code, iterations=mle.1d$iterations, 
      warnings=paste0(warn_strs, collapse="; "), 
      errors=paste0(stop_strs, collapse="; "),
      stringsAsFactors=FALSE)
  }
}

#' Return the likelihood value for a given set of parameters and observations
#' 
#' Called from mle_1ply(). From ?nlm, this function should be "the function to
#' be minimized, returning a single numeric value. This should be a function 
#' with first argument a vector of the length of p followed by any other 
#' arguments specified by the ... argument."
#' 
#' @param params a vector of length 2, where the first element is GPP and the 
#'   second element is ER (both mg/L/d)
#' @param K600.daily optional - if K600 is to be fit, leave NULL. If K600 is to
#'   be supplied, this is where to put it.
#' @param calc_DO_fun the function to use to build DO estimates from GPP, ER, 
#'   etc. default is calc_DO_mod, but could also be calc_DO_mod_by_diff
#' @keywords internal
negloglik_1ply <- function(params, K600.daily, DO.obs, DO.sat, depth, temp.water, 
                           frac.GPP, frac.ER, frac.D, calc_DO_fun) {
  
  # Count how many DO observations/predictions there should be
  n <- length(DO.obs)
  
  # Parse params vector (passed from nlm) and produce DO.mod estimates. It's
  # important that these arguments are named because various calc_DO_funs take 
  # slightly different subsets of these arguments
  DO.mod <- calc_DO_fun(
    GPP.daily=params[1], ER.daily=params[2], K600.daily=if(missing(K600.daily)) params[3] else K600.daily,
    DO.obs=DO.obs, DO.sat=DO.sat, depth=depth, temp.water=temp.water, 
    frac.GPP=frac.GPP, frac.ER=frac.ER, frac.D=frac.D, DO.mod.1=DO.obs[1], n=n)
  
  # calculate & return the negative log likelihood of DO.mod values relative
  # to DO.obs values. equivalent to Bob's original code & formula at
  # http://www.statlect.com/normal_distribution_maximum_likelihood.htm
  diffs.sq <- (DO.obs-DO.mod)^2 
  sigma.sq <- sum(diffs.sq)/n
  (n/2)*log(sigma.sq) + (n/2)*log(2*pi) + (1/(2*sigma.sq))*sum(diffs.sq)
}



#### metab_mle class ####

#' Metabolism model fitted by maximum likelihood estimation
#' 
#' \code{metab_mle} models use non-linear minimization of the negative log
#' likelihood to fit values of GPP, ER, and K for a given DO curve.
#' 
#' @exportClass metab_mle
#' @family metab.model.classes
setClass(
  "metab_mle", 
  contains="metab_model"
)

