#' Helper to predict_DO.metab_model
#' 
#' Usually assigned to model_fun within mm_model_by_ply, called from there
#' 
#' @inheritParams mm_model_by_ply_prototype
#' @param calc_DO_fun the function to use to build DO estimates from GPP, ER, 
#'   etc. default is calc_DO_mod, but could also be calc_DO_mod_by_diff
#' @param calc_DO_args a list of other arguments passed to calc_DO_fun
#' @return a data.frame of predictions
mm_predict_1ply <- function(
  data_ply, data_daily_ply, day_start, day_end, ply_date, timestep_days, ..., 
  calc_DO_fun, calc_DO_args=list()) {
  
  # The daily metabolism estimates are in data_daily_ply. Skip today (return
  # DO.mod=NAs) if they're missing. Otherwise, proceed to predict DO
  if(nrow(data_daily_ply)==0 || any(is.na(data_daily_ply[c("GPP","ER","K600")]))) {
    return(data.frame(data_ply, DO.mod=rep(NA, nrow(data_ply))))
  }
  
  # prepare arguments, appending calc_DO_args if there are any
  tot.GPP <- sum(data_ply$light[as.character(data_ply$solar.time,"%Y-%m-%d")==ply_date])
  frac.GPP <- if(tot.GPP > 0) data_ply$light/tot.GPP else 1/nrow(data_ply)
  all_args <- c(
    list(GPP.daily=data_daily_ply$GPP, ER.daily=data_daily_ply$ER, K600.daily=data_daily_ply$K600, 
         DO.obs=data_ply$DO.obs, DO.sat=data_ply$DO.sat, depth=data_ply$depth, temp.water=data_ply$temp.water, 
         frac.GPP=frac.GPP, frac.ER=timestep_days, frac.D=timestep_days, 
         DO.mod.1=if('DO.obs' %in% names(data_ply)) data_ply$DO.obs[1] else data_daily_ply$DO.mod.1,
         n=nrow(data_ply)),
    calc_DO_args)
  
  # call calc_DO_fun with the prepared arguments
  DO.mod <- do.call(calc_DO_fun, all_args)
  
  # return the data with modeled DO attached
  data.frame(data_ply, DO.mod=DO.mod)
}
