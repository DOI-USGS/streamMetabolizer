#' Helper to predict_DO.metab_model
#' 
#' Usually assigned to model_fun within mm_model_by_ply, called from there
#' 
#' @inheritParams mm_model_by_ply_prototype
#' @param model_name the coded model name that will determine the GPP_fun,
#'   ER_fun, deficit_src, and ode_method to use in prediction
#' @return a data.frame of predictions
mm_predict_DO_1ply <- function(
  data_ply, data_daily_ply, day_start, day_end, ply_date, ..., 
  model_name) {
  
  # the daily metabolism-relevant parameter estimates are in data_daily_ply.
  # skip today (return DO.mod=NAs) if they're missing. otherwise, proceed to
  # predict DO
  if(nrow(data_daily_ply)==0 || any(is.na(data_daily_ply))) {
    return(data.frame(data_ply, DO.mod=rep(NA, nrow(data_ply))))
  }

  # identify any observation and/or process error to be added (only known
  # application is for simulating data)
  err.obs <- if(exists('err.obs', data_ply)) data_ply$err.obs else 0
  err.proc <- if(exists('err.proc', data_ply)) data_ply$err.proc else 0
  
  # prepare DO prediction function
  features <- mm_parse_name(model_name)
  dDOdt <- create_calc_dDOdt(
    data_ply, ode_method=features$ode_method, GPP_fun=features$GPP_fun,
    ER_fun=features$ER_fun, deficit_src=features$deficit_src, # or maybe 'DO_mod' instead for prediction??
    err.proc=err.proc)
  DO <- create_calc_DO(dDOdt, ode_method=features$ode_method, err.obs=err.obs)
  
  # prepare arguments to DO prediction function. enable 3 possible sources for 
  # DO.mod.1: (1) fitted params or (2) data_daily (both passed via get_params to
  # data_daily_ply) or (3) data_ply$DO.obs[1]
  metab.pars <- data_daily_ply %>%
    select(-date) %>%
    bind_cols(if(!exists('DO.mod.1', data_daily_ply)) data.frame(DO.mod.1 = data_ply$DO.obs[1]) else NULL)
  
  # call DO prediction function
  DO.mod <- DO(metab.pars)
  
  # return the data with modeled DO attached
  data.frame(data_ply, DO.mod=DO.mod)
}
