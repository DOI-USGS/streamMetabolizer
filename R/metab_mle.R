#' @include metab_model-class.R
NULL

#' Maximum likelihood metabolism model fitting function
#' 
#' Uses maximum likelihood to fit a model to estimate GPP and ER from input data
#' on DO, temperature, light, etc.
#' 
#' @inheritParams metab
#' @return A metab_mle object containing the fitted model. This object can be 
#'   inspected with the functions in the \code{\link{metab_model_interface}}. 
#'   The \code{code} column in \code{get_fit(mm)} is defined in the Value
#'   subsection of \code{?nlm}.
#' @examples
#' dat <- data_metab('3','30')
#' # PRK
#' mm <- metab_mle(data=dat)
#' predict_metab(mm)
#' 
#' # PR with fixed K on two days
#' dat_daily <- data.frame(date=as.Date(c("2012-09-18","2012-09-20")), K600.daily=35)
#' metab_mle(data=dat, data_daily=dat_daily)
#' 
#' # PRK with day-specific inits on some days
#' dat_daily <- data.frame(date=as.Date("2012-09-19"), 
#'   init.GPP.daily=4, init.K600.daily=60)
#' metab_mle(data=dat, data_daily=dat_daily)
#' 
#' # Nonlinear GPP or ER equations
#' metab_mle(specs(mm_name('mle', GPP_fun='satlight')), data=dat)
#' metab_mle(specs(mm_name('mle', ER_fun='q10temp')), data=dat)
#' \dontrun{
#' plot_DO_preds(predict_DO(mm))
#' }
#' @export
#' @family metab_model
metab_mle <- function(
  specs=specs(mm_name('mle')),
  data=mm_data(solar.time, DO.obs, DO.sat, depth, temp.water, light),
  data_daily=mm_data(date, K600.daily, init.GPP.daily, init.ER.daily, init.K600.daily, optional='all'),
  info=NULL
) {
  
  if(missing(specs)) {
    # if specs is left to the default, it gets confused about whether specs() is
    # the argument or the function. tell it which:
    specs <- streamMetabolizer::specs(mm_name('mle'))
  }
  if(mm_parse_name(specs$model_name)$ode_method %in% c('lsoda','lsodes','lsodar'))
    warning("we've seen bad results with ODE methods 'lsoda', 'lsodes', and 'lsodar'. Use at your own risk")
  
  fitting_time <- system.time({
    # Check data for correct column names & units
    dat_list <- mm_validate_data(if(missing(data)) NULL else data, if(missing(data_daily)) NULL else data_daily, "metab_mle")
    data <- v(dat_list[['data']])
    data_daily <- v(dat_list[['data_daily']])
    
    # model the data, splitting into overlapping 31.5-hr 'plys' for each date
    mle_all <- mm_model_by_ply(
      mle_1ply, data=data, data_daily=data_daily, # for mm_model_by_ply
      day_start=specs$day_start, day_end=specs$day_end, day_tests=specs$day_tests, # for mm_model_by_ply
      specs=specs) # for mle_1ply and create_calc_dDOdt
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
  
  # Update data with DO & metab predictions
  mm@data <- predict_DO(mm)
  mm@metab_daily <- predict_metab(mm)
  
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
  specs=specs('m_np_oi_tr_km.nlm') # inheritParams metab
) {
  
  # Provide ability to skip a poorly-formatted day for calculating
  # metabolism, without breaking the whole loop. Just collect
  # problems/errors as a list of strings and proceed. Also collect warnings.
  stop_strs <- if(isTRUE(ply_validity)) character(0) else ply_validity
  warn_strs <- character(0)
  
  if(length(stop_strs) == 0) {
    # Collect K600 and date-specific initial values if they're available. If a
    # value is named in data_daily_ply but not available (nrow(data_daily)==0 ||
    # value==NA), use the default: for xx_init values this is specs$xx_init, for
    # K600 this is NULL (fit by MLE)
    init.vals <- unlist(specs[grepl('^init.', names(specs))])
    K600 <- NULL
    . <- '.dplyr.var'
    if(!is.null(data_daily_ply)) {
      daily.inits <- names(init.vals) %>% { .[. %in% names(data_daily_ply)] }
      if(nrow(data_daily_ply)==0) {
        if(length(daily.inits) > 0) {
          warn_strs <- c(warn_strs, paste0(
            "nrow(data_daily)==0 so using specs for ", paste0(daily.inits, collapse=', ')))
        }
        if(exists('K600', data_daily_ply)) {
          warn_strs <- c(warn_strs, "nrow(data_daily)==0 so fitting K600 by MLE")
        }
      } else { # nrow(data_daily_ply)==1 because we're in mle_1ply, so it can only be 0 or 1
        for(init in daily.inits) {
          if(is.na(data_daily_ply[[init]])) 
            warn_strs <- c(warn_strs, paste0("data_daily$", init, "==NA so using specs"))
          else
            init.vals[[init]] <- data_daily_ply[[init]]
        }
        if(exists('K600', data_daily_ply)) {
          if(is.na(data_daily_ply$K600)) {
            warn_strs <- c(warn_strs, "data_daily$K600==NA so fitting by MLE")
          } else {
            K600 <- data_daily_ply$K600
          }
        }
      }
    }
    fix_K600 <- !is.null(K600)
    if(fix_K600) init.vals[['init.K600.daily']] <- NULL # if K600.daily will be fixed, we don't need to fit it
    
    # name the init.vals as their parameters
    init.vals <- setNames(unlist(init.vals), substring(names(init.vals), 6))
    
    # parse the model_name
    features <- mm_parse_name(specs$model_name)
    
    # create the series of nested functions to compute the negative log 
    # likelihood (NLL) for a trio of values for GPP.daily, ER.daily, and
    # K600.daily
    dDOdt <- create_calc_dDOdt(
      data_ply, ode_method=features$ode_method, GPP_fun=features$GPP_fun,
      ER_fun=features$ER_fun, deficit_src=features$deficit_src)
    DO <- create_calc_DO(dDOdt)
    # to fit DO.mod.1, err_obs_iid_sigma, and/or err_proc_iid_sigma, add these to par.names in create_calc_NLL and nlm.args$p
    NLL <- create_calc_NLL(DO, err_obs_iid=features$err_obs_iid, err_proc_iid=features$err_proc_iid)
    if(fix_K600) environment(NLL)$par.names %<>% { .[. != 'K600.daily'] } # remove the K600 parameter if we're fixing K600
    
    # package nlm arguments in a list. estimate of fscale is based on comparison
    # of dat=data_metab('10','5') and dat=data_metab('10','30') in mm <- 
    # metab_mle(specs(mm_name('mle', GPP_fun='satlight')), data=dat), for which
    # typical minima are -300 and -70, and obs/day are 288 and 48, respectively
    nlm.args <- c(
      list(
        f = NLL,
        p = init.vals,
        typsize = init.vals,
        fscale = -(1/timestep_days),
        hessian = TRUE),
      if(fix_K600) list(K600.daily = K600)
    )
    
    # Find the best metabolism values by non-linear minimization of the likelihood
    # (NLL) of the output from a DO-prediction function, which is the ODE solution
    # of predictions from a dDO/dt-prediction function of a trial set of
    # metabolism values
    mle.1d <- withCallingHandlers(
      tryCatch({
        # try to run the MLE fitting function
        # iter <<- NA # only useful if create_calc_NLL has plotting uncommented
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
        
        # return from tryCatch
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
  val.names <- names(init.vals) # trust nlm to return the parameters in the same order we passed them in
  stat.names <- c('estimate','sd','gradient')
  valstat.names <- paste0(rep(val.names, each=length(stat.names)), rep(c('','.sd','.grad'), times=length(val.names)))
  goodness.names <- c('minimum','iterations','code')
  if(length(stop_strs) > 0) {
    valstat.cols <- as.list(rep(as.numeric(NA), length(valstat.names))) %>%
      as.data.frame() %>%
      setNames(valstat.names)
    goodness.cols <- as.list(rep(as.numeric(NA), length(goodness.names))) %>%
      as.data.frame() %>%
      setNames(goodness.names) %>%
      mutate(code.str=as.character(NA))
  } else {
    valstat.cols <- c(do.call(rbind, mle.1d[stat.names])) %>%
      as.list() %>% as.data.frame() %>%
      setNames(valstat.names)
    goodness.cols <- as.data.frame(mle.1d[goodness.names]) %>%
      # code interpretations are from the Value section of the ?nlm page
      mutate(code.str=c(
        'relative gradient is close to zero, current iterate is probably solution',
        'successive iterates within tolerance, current iterate is probably solution',
        'last global step failed to locate a point lower than estimate. Either estimate is an approximate local minimum of the function or steptol is too small',
        'iteration limit exceeded',
        'maximum step size stepmax exceeded five consecutive times. Either the function is unbounded below, becomes asymptotic to a finite value from above in some direction or stepmax is too small'
      )[code])
    if(goodness.cols$code > 2) warn_strs <- c(warn_strs, goodness.cols$code.str)
  }
  err.cols <- data.frame(
    warnings=paste0(unique(warn_strs), collapse="; "),
    errors=paste0(unique(stop_strs), collapse="; "),
    stringsAsFactors=FALSE)
  data.frame(valstat.cols, goodness.cols, err.cols)
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

