#' Helper to predict_DO.metab_model
#' 
#' Usually assigned to model_fun within mm_model_by_ply, called from there
#' 
#' @inheritParams mm_model_by_ply_prototype
#' @param model_name the coded model name that will determine the GPP_fun,
#'   ER_fun, deficit_src, and ode_method to use in prediction
#' @return a data.frame of predictions
#' @importFrom stats rnorm
mm_predict_DO_1ply <- function(
  data_ply, data_daily_ply, day_start, day_end, ply_date, ..., 
  model_name) {
  
  # the daily metabolism-relevant parameter estimates are in data_daily_ply.
  # skip today (return DO.mod=NAs) if they're missing. otherwise, proceed to
  # predict DO
  if(nrow(data_daily_ply)==0 || any(is.na(data_daily_ply))) {
    return(data.frame(data_ply, DO.mod=rep(NA, nrow(data_ply))))
  }

  # prepare arguments to DO prediction function. enable 3 possible sources for 
  # DO.mod.1: (1) fitted params or (2) data_daily (both passed via get_params to
  # data_daily_ply) or (3) data_ply$DO.obs[1]
  metab.pars <- data_daily_ply %>%
    select(-date) %>%
    bind_cols(if(!exists('DO.mod.1', data_daily_ply)) data.frame(DO.mod.1 = data_ply$DO.obs[1]) else NULL)
  
  # get info on the model structure
  features <- mm_parse_name(model_name)
  
  # identify any observation or process error to be added, and in what combos
  if(features$type == 'sim') {
    data_ply$DO.obs <- NULL # remove the old one; we're gonna replace it. we've already used it for metab.pars above if needed
    n <- nrow(data_ply)
    err.obs <- as.numeric(stats::filter(rnorm(n, 0, data_daily_ply$err.obs.sigma), filter=data_daily_ply$err.obs.phi, method="recursive"))
    err.proc <- as.numeric(stats::filter(rnorm(n, 0, data_daily_ply$err.proc.sigma), filter=data_daily_ply$err.proc.phi, method="recursive"))      
    errs <- list(
      list(colname='DO.pure', err.obs=0, err.proc=0), # the DO implied by the daily metab params (error-free)
      list(colname='DO.mod', err.obs=0, err.proc=err.proc), # the 'true' DO (with proc err)
      list(colname='DO.obs', err.obs=err.obs, err.proc=err.proc)) # the 'observed' DO (with obs err)
  } else {
    errs <- list(
      list(colname='DO.mod', err.obs=0, err.proc=0))
  }
  
  for(e in errs) {
    # prepare DO prediction function
    dDOdt <- create_calc_dDOdt(
      data_ply, ode_method=features$ode_method, GPP_fun=features$GPP_fun,
      ER_fun=features$ER_fun, deficit_src=features$deficit_src,
      err.proc=e$err.proc)
    DO <- create_calc_DO(dDOdt, ode_method=features$ode_method, err.obs=e$err.obs)
    
    # call DO prediction function and assign to named column
    data_ply[e$colname] <- DO(metab.pars)
  }
  
  # return the data with modeled DO attached
  data_ply
}
