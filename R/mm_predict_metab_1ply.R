#' Helper to predict_metab.metab_model
#' 
#' Usually assigned to model_fun within mm_model_by_ply, called from there
#' 
#' @inheritParams mm_model_by_ply_prototype
#' @param model_name the coded model name that will determine the GPP_fun,
#'   ER_fun, deficit_src, and ode_method to use in prediction
#' @return a data.frame of predictions
mm_predict_metab_1ply <- function(
  data_ply, data_daily_ply, day_start, day_end, ply_date, ..., 
  model_name) {
  
  # the daily metabolism-relevant parameter estimates are in data_daily_ply.
  # skip today (return DO.mod=NAs) if they're missing. otherwise, proceed to
  # predict DO
  if(nrow(data_daily_ply)==0) {
    return(data.frame(data_ply, GPP=NA))
  }
  
  # prepare metab prediction function
  features <- mm_parse_name(model_name)
  dDOdt <- create_calc_dDOdt(
    data_ply, ode_method=features$ode_method, GPP_fun=features$GPP_fun,
    ER_fun=features$ER_fun, deficit_src=features$deficit_src, # or maybe 'DO_mod' instead for prediction??
    err.proc=0)
  t <- environment(dDOdt)$data$t
  
  # call prediction functions and/or compute CIs for GPP, ER, & D
  preds <- bind_cols(lapply(c('GPP','ER'), function(met) { # lapply can also apply to 'D', but leaving out for now
    precalc.names <- paste0(met, '.daily', c('','.sd'))
    if(all(precalc.names %in% names(data_daily_ply))) {
      met.preds <- mm_sd_to_ci(data_daily_ply[precalc.names])
    } else {
      met.preds <- tryCatch({
        data.frame(
          mean(mapply(
            function(t, DO.mod.t) {
              switch(
                met,
                GPP=, ER=environment(dDOdt)[[met]](t=t, metab.pars=data_daily_ply),
                D=environment(dDOdt)[[met]](t=t, metab.pars=data_daily_ply, DO.mod.t=DO.mod.t) * environment(dDOdt)$data$depth
              )
            }, t=t, DO.mod.t=data_ply$DO.mod)),
          lower = NA, # could use delta method to derive CIs here
          upper = NA)},
        error=function(e) NA)
    }
    met.preds %>% setNames(paste0(met, c("",".lower",".upper")))
  }))
  
  # return the modeled daily mean metabolism and reaeration rates
  preds
}
