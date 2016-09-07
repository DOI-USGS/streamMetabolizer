#' Helper to predict_metab.metab_model
#' 
#' Usually assigned to model_fun within mm_model_by_ply, called from there
#' 
#' @inheritParams mm_model_by_ply_prototype
#' @param model_name the coded model name that will determine the GPP_fun,
#'   ER_fun, deficit_src, and ode_method to use in prediction
#' @return a data.frame of predictions
mm_predict_metab_1ply <- function(
  data_ply, data_daily_ply, day_start, day_end, ply_date, ply_validity, ..., 
  model_name) {
  
  # record date test failures as warnings (flipped from usual model-fitting
  # practice of recording these as errors)
  warn_strs <- if(isTRUE(ply_validity)) character(0) else ply_validity
  stop_strs <- character(0)
  
  # skip today if we're missing metabolism estimates and/or input data (return a
  # near-empty or empty data.frame, respectively)
  if(nrow(data_daily_ply) == 0) {
    na_num <- as.numeric(NA)
    na_df <- data.frame(
      date=ply_date, 
      GPP=na_num, GPP.lower=na_num, GPP.upper=na_num,
      ER=na_num, ER.lower=na_num, ER.upper=na_num)
    if(nrow(data_ply) > 0) {
      # in the more common case, we just need a row of NAs to reflect that
      # the daily metabolism-relevant parameter estimates are missing
      return(na_df)
    } else {
      # if mm_model_by_ply comes up with 0 data rows total, it calls this 
      # function once with empty data_daily_ply and data_ply to figure out which
      # column names to return. Give those names here as a rowless data.frame
      return(na_df[c(),])
    }
  }
  
  # prepare metab prediction functions
  features <- mm_parse_name(model_name)
  env.dDOdt <- 
    create_calc_dDOdt(
      data_ply, ode_method=features$ode_method, GPP_fun=features$GPP_fun,
      ER_fun=features$ER_fun, deficit_src=features$deficit_src, # or maybe 'DO_mod' instead for prediction??
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
