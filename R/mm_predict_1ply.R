#' Helper to predict_DO.metab_model
#' 
#' Usually assigned to model_fun within mm_model_by_ply, called from there
#' 
#' @inheritParams mm_model_by_ply_prototype
#' @param model_name the coded model name that will determine the GPP_fun,
#'   ER_fun, deficit_src, and ode_method to use in prediction
#' @return a data.frame of predictions
mm_predict_1ply <- function(
  data_ply, data_daily_ply, day_start, day_end, ply_date, ..., 
  model_name) {
  
  # The daily metabolism estimates are in data_daily_ply. Skip today (return
  # DO.mod=NAs) if they're missing. Otherwise, proceed to predict DO
  if(nrow(data_daily_ply)==0 || any(is.na(data_daily_ply[c("GPP","ER","K600")]))) {
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
  
  # call DO prediction function
  DO.mod <- data_daily_ply %>%
    rename(GPP.daily=GPP, ER.daily=ER, K600.daily=K600) %>%
    select(-date) %>%
    DO()
  
  # return the data with modeled DO attached
  data.frame(data_ply, DO.mod=DO.mod)
}
