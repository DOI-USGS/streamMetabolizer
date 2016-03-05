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
#' # set the date in several formats
#' start.chron <- chron::chron(dates="08/23/12", times="22:00:00")
#' end.chron <- chron::chron(dates="08/25/12", times="06:00:00")
#' start.posix <- as.POSIXct(format(start.chron, "%Y-%m-%d %H:%M:%S"), tz="UTC")
#' end.posix <- as.POSIXct(format(end.chron, "%Y-%m-%d %H:%M:%S"), tz="UTC")
#' mid.date <- as.Date(start.posix + (end.posix - start.posix)/2, tz=lubridate::tz(start.posix))
#' start.numeric <- as.numeric(start.posix - as.POSIXct(format(mid.date, "%Y-%m-%d 00:00:00"),
#'    tz="UTC"), units='hours')
#' end.numeric <- as.numeric(end.posix - as.POSIXct(format(mid.date, "%Y-%m-%d 00:00:00"),
#'   tz="UTC"), units='hours')
#' 
#' # get, format, & subset data
#' vfrench <- streamMetabolizer:::load_french_creek(attach.units=FALSE)
#' vfrenchshort <- vfrench[vfrench$solar.time >= start.posix & vfrench$solar.time <= end.posix, ]
#' 
#' # PRK
#' get_fit(mm <- metab_mle(data=vfrenchshort, day_start=start.numeric, 
#'   day_end=end.numeric))[2,c("GPP","ER","K600","minimum")]
#' plot_DO_preds(predict_DO(mm))
#' streamMetabolizer:::load_french_creek_std_mle(vfrenchshort, estimate='PRK')
#' 
#' # PR
#' get_fit(mm <- metab_mle(data=vfrenchshort, data_daily=data.frame(date=mid.date, K600=35), 
#'   day_start=start.numeric, day_end=end.numeric))[2,c("GPP","ER","K600","minimum")]
#' get_fitting_time(mm)
#' plot_DO_preds(predict_DO(mm))
#' streamMetabolizer:::load_french_creek_std_mle(vfrenchshort, estimate='PR', K=35)
#' 
#' \dontrun{
#'   metab_mle(data=data.frame(empty="shouldbreak"))
#'  
#'   # PRK and PR with process error
#'   get_fit(mm <- metab_mle(data=vfrenchshort, 
#'     model_specs=specs('m_np_pi_pm_km.nlm'), 
#'     day_start=start.numeric, day_end=end.numeric))[2,c("GPP","ER","K600","minimum")]
#'   plot_DO_preds(predict_DO(mm))
#'   get_fit(mm <- metab_mle(data=vfrenchshort, data_daily=data.frame(date=mid.date, K600=35), 
#'     model_specs=specs('m_np_pi_pm_km.nlm'), 
#'     day_start=start.numeric, day_end=end.numeric))[2,c("GPP","ER","K600","minimum")]
#'   plot_DO_preds(predict_DO(mm))
#' }
#' @export
#' @family metab_model
metab_mle <- function(
  model_specs=specs(mm_name('mle')),
  data=mm_data(solar.time, DO.obs, DO.sat, depth, temp.water, light), 
  data_daily=mm_data(date, K600, optional='all'), 
  info=NULL,
  day_start=4, day_end=27.99, tests=c('full_day', 'even_timesteps', 'complete_data')
) {
  
  fitting_time <- system.time({
    # Check data for correct column names & units
    dat_list <- mm_validate_data(data, if(missing(data_daily)) NULL else data_daily, "metab_mle")
    
    # model the data, splitting into overlapping 31.5-hr 'plys' for each date
    mle_all <- mm_model_by_ply(
      mle_1ply, data=dat_list[['data']], data_daily=dat_list[['data_daily']], # for mm_model_by_ply
      day_start=day_start, day_end=day_end, # for mm_model_by_ply and mm_is_valid_day
      tests=tests, # for mm_is_valid_day
      model_specs=model_specs) # for mle_1ply and negloglik_1ply
  })
  
  # Package results
  mm <- metab_model(
    model_class="metab_mle",
    info=info,
    fit=mle_all,
    fitting_time=fitting_time,
    args=list(model_specs=model_specs, day_start=day_start, day_end=day_end, tests=tests), # keep in order passed to function
    data=dat_list[['data']],
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
#' @inheritParams mm_is_valid_day
#' @inheritParams metab
#' @return data.frame of estimates and \code{\link[stats]{nlm}} model 
#'   diagnostics
#' @importFrom stats nlm
#' @keywords internal
mle_1ply <- function(
  data_ply, data_daily_ply, day_start, day_end, ply_date, # inheritParams mm_model_by_ply_prototype
  tests=c('full_day', 'even_timesteps', 'complete_data'), # inheritParams mm_is_valid_day
  model_specs=specs('m_np_oi_pm_km.nlm')
) {
  
  # Provide ability to skip a poorly-formatted day for calculating 
  # metabolism, without breaking the whole loop. Just collect 
  # problems/errors as a list of strings and proceed. Also collect warnings.
  timestep.days <- suppressWarnings(mean(as.numeric(diff(v(data_ply$solar.time)), units="days"), na.rm=TRUE))
  validity <- mm_is_valid_day(
    data_ply, # data split by mm_model_by_ply
    tests=tests, day_start=day_start, day_end=day_end, # args passed from metab_mle
    timestep_days=timestep.days) # arg supplied here to avoid calculating twice
  stop_strs <- if(isTRUE(validity)) character(0) else validity
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
        p = c(GPP=model_specs$GPP_init, ER=model_specs$ER_init, K600=model_specs$K600_init)[if(is.null(K600)) 1:3 else 1:2],
        hessian = TRUE,
        K600.daily=K600
      ),
      as.list(
        data_ply[c("DO.obs","DO.sat","depth","temp.water")]
      ),
      list(
        frac.GPP = data_ply$light/sum(data_ply$light[as.character(data_ply$solar.time,"%Y-%m-%d")==as.character(ply_date)]),
        frac.ER = timestep.days,
        frac.D = timestep.days,
        calc_DO_fun = model_specs$calc_DO_fun,
        ODE_method = model_specs$ODE_method
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
      warnings=paste0(warn_strs, collapse="; "), 
      errors=paste0(stop_strs, collapse="; "),
      stringsAsFactors=FALSE)
  } else {
    data.frame(
      GPP=mle.1d$estimate[1], ER=mle.1d$estimate[2], K600=if(is.null(K600)) mle.1d$estimate[3] else K600,
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

