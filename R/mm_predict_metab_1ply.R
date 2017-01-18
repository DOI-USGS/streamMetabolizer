#' Helper to predict_metab.metab_model
#' 
#' Usually assigned to model_fun within mm_model_by_ply, called from there
#' 
#' @inheritParams mm_model_by_ply_prototype
#' @param model_name the coded model name that will determine the GPP_fun,
#'   ER_fun, deficit_src, and ode_method to use in prediction
#' @return a data.frame of predictions
#' @import dplyr
mm_predict_metab_1ply <- function(
  data_ply, data_daily_ply, day_start, day_end, ply_date, ply_validity, ..., 
  model_name) {
  
  # record ply validity failures if present (only gets tested for sim models)
  stop_strs <- if(isTRUE(ply_validity)) character(0) else ply_validity
  warn_strs <- character(0)
  
  # skip today if we're missing metabolism estimates and/or input data (return a
  # near-empty or empty data.frame, respectively, and don't bother reporting on 
  # ply validity). order the checks such that skip_for gets set to the most
  # important reason we're skipping the date (if we are)
  skip_for <- c()
  if(nrow(data_daily_ply) == 0) skip_for <- 'empty'
  else if(data_daily_ply %>% select(-date) %>% sapply(is.na) %>% all()) skip_for <- 'all_NA'
  else if(length(stop_strs) > 0) skip_for <- 'validity'
  if(length(skip_for) > 0) {
    na.df <- data.frame(
      GPP=NA_real_, GPP.lower=NA_real_, GPP.upper=NA_real_,
      ER=NA_real_, ER.lower=NA_real_, ER.upper=NA_real_, 
      warnings=NA_character_, errors=NA_character_, 
      stringsAsFactors=FALSE)
    out.df <- switch(
      skip_for,
      'validity'=na.df %>% mutate(
        warnings=paste0(unique(warn_strs), collapse="; "),
        errors=paste0(unique(stop_strs), collapse="; ")),
      'all_NA'=na.df,
      'empty'=filter(na.df, FALSE))
    return(out.df)
  }
  
  # prepare metab prediction functions
  features <- mm_parse_name(model_name)
  env.dDOdt <- 
    create_calc_dDOdt(
      data_ply, ode_method=features$ode_method, GPP_fun=features$GPP_fun,
      ER_fun=features$ER_fun, deficit_src=features$deficit_src,
      err.proc=0) %>%
    environment()
  t <- env.dDOdt$data$t
  
  # call prediction functions and/or compute CIs for GPP, ER, & D
  preds <- bind_cols(lapply(c('GPP','ER'), function(met) { # lapply can also apply to 'D', but leaving out for now
    precalc.names <- paste0(met, '.daily', c('','.sd'))
    if(all(precalc.names %in% names(data_daily_ply))) {
      # warning: this shortcut won't work if we permit time periods other than
      # 24 hours!!
      met.preds <- mm_sd_to_ci(data_daily_ply[precalc.names])
    } else {
      met.preds <- tryCatch({
        data.frame(
          mean(mapply(
            function(t, DO.mod.t) {
              switch(
                met,
                GPP=, ER=env.dDOdt[[met]](t=t, metab.pars=data_daily_ply),
                D=env.dDOdt[[met]](t=t, metab.pars=data_daily_ply, DO.mod.t=DO.mod.t) * env.dDOdt$data$depth
              )
            }, t=t, DO.mod.t=data_ply$DO.mod)),
          lower = NA, # could use delta method to derive CIs here
          upper = NA)},
        error=function(e) NA)
    }
    met.preds %>% setNames(paste0(met, c("",".lower",".upper")))
  }))
  
  # return the modeled daily mean metabolism and reaeration rates
  err.cols <- data.frame(
    warnings=paste0(unique(warn_strs), collapse="; "),
    errors=paste0(unique(stop_strs), collapse="; "), # so far there will never be error strings; just a placeholder
    stringsAsFactors=FALSE)
  data.frame(preds, err.cols)
}
