#' Helper to predict_DO.metab_model
#' 
#' Usually assigned to model_fun within mm_model_by_ply, called from there
#' 
#' @inheritParams mm_model_by_ply_prototype
#' @param calc_DO_fun the function to use to build DO estimates from GPP, ER, 
#'   etc. default is calc_DO_mod, but could also be calc_DO_mod_by_diff
#' @return a data.frame of predictions
mm_predict_1ply <- function(data_ply, data_daily_ply, day_start, day_end, local_date, calc_DO_fun) {
  
  # get the daily metabolism estimates, and skip today (return DO.mod=NAs) if
  # they're missing
  metab_est <- data_daily_ply
  if(nrow(metab_est)==0 || is.na(metab_est$GPP)) {
    return(data.frame(data_ply, DO.mod=rep(NA, nrow(data_ply))))
  }
  
  # if we have metab estimates, use them to predict DO
  . <- local.time <- ".dplyr.var"
  data_ply %>%
    do(with(., {
      
      # prepare auxiliary data
      n <- length(local.time)
      timestep.days <- suppressWarnings(mean(as.numeric(diff(local.time), units="days"), na.rm=TRUE))
      frac.GPP <- light/sum(light[as.character(local.time,"%Y-%m-%d")==local_date])
      
      # produce DO.mod estimates for today's GPP and ER
      DO.mod <- calc_DO_fun(
        GPP.daily=metab_est$GPP, 
        ER.daily=metab_est$ER, 
        K600.daily=metab_est$K600, 
        DO.obs=DO.obs, DO.sat=DO.sat, depth=depth, temp.water=temp.water, 
        frac.GPP=frac.GPP, frac.ER=timestep.days, frac.D=timestep.days, DO.mod.1=DO.obs[1], n=n)
      
      data.frame(., DO.mod=DO.mod)
    }))
}
