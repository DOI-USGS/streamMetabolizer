#' @include metab_model-class.R
NULL

#' Maximum likelihood metabolism model fitting function
#' 
#' Uses maximum likelihood to fit a model to estimate GPP and ER from input data
#' on DO, temperature, light, etc.
#' 
#' @author Alison Appling, Jordan Read, Luke Winslow
#'   
#' @inheritParams metab
#' @return A metab_mle object containing the fitted model. This object can be 
#'   inspected with the functions in the \code{\link{metab_model_interface}}.
#'   
#' @examples
#' dat <- data_metab('3')
#' # PRK
#' mm <- metab_mle(data=dat)
#' predict_metab(mm)
#' 
#' # PR with fixed K on two days
#' dat_daily <- data.frame(date=as.Date(c("2012-09-18","2012-09-20")), K600=35)
#' mm <- metab_mle(data=dat, data_daily=dat_daily)
#' predict_metab(mm)
#' \dontrun{
#' plot_DO_preds(predict_DO(mm))
#' }
#' @export
#' @family metab_model
metab_mle <- function(
  specs=specs(mm_name('mle')),
  data=mm_data(solar.time, DO.obs, DO.sat, depth, temp.water, light), 
  data_daily=mm_data(date, K600, optional='all'), 
  info=NULL
) {
  
  if(missing(specs)) {
    # if specs is left to the default, it gets confused about whether specs() is
    # the argument or the function. tell it which:
    specs <- streamMetabolizer::specs(mm_name('mle'))
  }
  fitting_time <- system.time({
    # Check data for correct column names & units
    dat_list <- mm_validate_data(if(missing(data)) NULL else data, if(missing(data_daily)) NULL else data_daily, "metab_mle")
    data <- v(dat_list[['data']])
    data_daily <- v(dat_list[['data_daily']])
    
    # model the data, splitting into overlapping 31.5-hr 'plys' for each date
    mle_all <- mm_model_by_ply(
      mle_1ply, data=data, data_daily=data_daily, # for mm_model_by_ply
      day_start=specs$day_start, day_end=specs$day_end, day_tests=specs$day_tests, # for mm_model_by_ply
      specs=specs) # for mle_1ply and negloglik_1ply
  })
  
  # Package results
  mm <- metab_model(
    model_class="metab_mle",
    info=info,
    fit=mle_all,
    fitting_time=fitting_time,
    specs=specs,
    data=dat_list[['data']], # keep the units if given
    data_daily=dat_list[['data_daily']])
  
  # Update data with DO predictions
  mm@data <- predict_DO(mm)
  
  # Return
  mm
}


#### helpers ####

#' Make daily metabolism estimates from input parameters
#' 
#' Called from metab_mle().
#' 
#' @inheritParams mm_model_by_ply_prototype
#' @inheritParams metab
#' @return data.frame of estimates and \code{\link[stats]{nlm}} model 
#'   diagnostics
#' @importFrom stats nlm
#' @keywords internal
mle_1ply <- function(
  data_ply, data_daily_ply, ply_date, ply_validity, timestep_days, ..., # inheritParams mm_model_by_ply_prototype
  specs=specs('m_np_oi_pm_km.nlm') # inheritParams metab
) {
  
  # Provide ability to skip a poorly-formatted day for calculating 
  # metabolism, without breaking the whole loop. Just collect 
  # problems/errors as a list of strings and proceed. Also collect warnings.
  stop_strs <- if(isTRUE(ply_validity)) character(0) else ply_validity
  warn_strs <- character(0)

  # Collect K600 if it's available
  K600 <- if(is.null(data_daily_ply)) {
    NULL 
  } else {
    if(nrow(data_daily_ply)==0 || is.na(data_daily_ply$K600)) {
      # Daily K600 has generally been supplied but isn't available for this 
      # particular day. What do we do? Do we (A) estimate all three parameters 
      # because we can (K600 <- NULL), or (B) omit this day because the user
      # will be expecting consistency in methods across days (K600 <- {stop_strs
      # <- c(stop_strs, "data_daily$K600 is NA"); NA})? Going with (C) give
      # warning and then estimate all three parameters.
      warn_strs <- c(warn_strs, "data_daily$K600 is NA so fitting by MLE")
      NULL #(C)
    } else {
      data_daily_ply$K600
    }
  }
  
  # Calculate metabolism by non linear minimization of an MLE function
  if(length(stop_strs) == 0) {
    nlm.args <- c(
      list(
        f = negloglik_1ply,
        p = c(GPP=specs$GPP_init, ER=specs$ER_init, K600=specs$K600_init)[if(is.null(K600)) 1:3 else 1:2],
        hessian = TRUE,
        K600.daily=K600
      ),
      as.list(
        data_ply[c("DO.obs","DO.sat","depth","temp.water")]
      ),
      list(
        frac.GPP = data_ply$light/sum(data_ply$light[as.character(data_ply$solar.time,"%Y-%m-%d")==as.character(ply_date)]),
        frac.ER = timestep_days,
        frac.D = timestep_days,
        calc_DO_fun = specs$calc_DO_fun,
        ODE_method = specs$ODE_method
      ))
    
    mle.1d <- withCallingHandlers(
      tryCatch({
        # first: try to run the MLE fitting function
        mle.1d <- do.call(nlm, nlm.args)
        # if we were successful, also compute the confidence interval as in 
        # http://www.stat.umn.edu/geyer/5931/mle/mle.pdf section 2.3, which
        # says: 'Inverse Fisher information gives the asymptotic variance matrix
        # of the MLE. From it, we can construct asymptotic confidence
        # intervals.' See also 
        # http://stats.stackexchange.com/questions/27033/in-r-given-an-output-from-optim-with-a-hessian-matrix-how-to-calculate-paramet
        # This will be a gross underestimate of true uncertainty because the
        # errors are autocorrelated. Still, it's a start.
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
      warnings=paste0(unique(warn_strs), collapse="; "), 
      errors=paste0(unique(stop_strs), collapse="; "),
      stringsAsFactors=FALSE)
  } else {
    data.frame(
      GPP=mle.1d$estimate[1], ER=mle.1d$estimate[2], K600=if(is.null(K600)) mle.1d$estimate[3] else K600,
      GPP.sd=mle.1d$sd[1], ER.sd=mle.1d$sd[2], K600.sd=if(is.null(K600)) mle.1d$sd[3] else NA,
      GPP.grad=mle.1d$gradient[1], ER.grad=mle.1d$gradient[2], K600.grad=if(is.null(K600)) mle.1d$gradient[3] else NA,
      minimum=mle.1d$minimum, code=mle.1d$code, iterations=mle.1d$iterations, 
      warnings=paste0(unique(warn_strs), collapse="; "), 
      errors=paste0(unique(stop_strs), collapse="; "),
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
#' @inheritParams calc_DO_mod_w_fixed_error
#' @keywords internal
negloglik_1ply <- function(params, K600.daily, DO.obs, DO.sat, depth, temp.water, 
                           frac.GPP, frac.ER, frac.D, calc_DO_fun, ODE_method) {
  
  # Count how many DO observations/predictions there should be
  n <- length(DO.obs)
  
  # Parse params vector (passed from nlm) and produce DO.mod estimates. It's
  # important that these arguments are named because various calc_DO_funs take 
  # slightly different subsets of these arguments
  DO.mod <- calc_DO_fun(
    GPP.daily=params[1], ER.daily=params[2], K600.daily=if(is.null(K600.daily)) params[3] else K600.daily,
    DO.obs=DO.obs, DO.sat=DO.sat, depth=depth, temp.water=temp.water, 
    frac.GPP=frac.GPP, frac.ER=frac.ER, frac.D=frac.D, DO.mod.1=DO.obs[1], n=n, ODE_method=ODE_method)
  
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

