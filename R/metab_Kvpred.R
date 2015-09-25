#' @include metab_model-class.R
NULL

#' Regress initial daily K estimates against predictors
#' 
#' Takes daily estimates of K, usually from nighttime regression, and regress 
#' against predictors such as discharge. Returns a metab_model that only 
#' predicts daily K, nothing else.
#' 
#' @author Alison Appling
#' @inheritParams metab_model_prototype
#' @param method the name of an interpolation or regression method relating K to the predictor[s] of choice.
#' @return A metab_Kvpred object containing the fitted model.
#' @examples
#' \dontrun{
#'  metab_Kvpred(data=data.frame(empty="shouldbreak"))
#' }
#' @export
#' @family metab_model
metab_Kvpred <- function(
  data=mm_data(NULL), data_daily=mm_data(local.date, K600, discharge.daily, velocity.daily), info=NULL, # inheritParams metab_model_prototype
  method="KvQ_regression"
) {
  
  fitting_time <- system.time({
    # Check data for correct column names & units
    dat_list <- mm_validate_data(if(missing(data)) NULL else data, data_daily, "metab_Kvpred")
    data <- dat_list[['data']]
    data_daily <- dat_list[['data_daily']]
    
    # model the data all days at a time, after first filtering out bad days
    data_daily_filtered <- data_daily ## need to finish
    warning("need to actually filter data_daily")
    Kvpred_all <- Kvpred_allply(
      data_daily_filtered)
  })
  
  # Package and return results
  metab_model(
    "metab_Kvpred", 
    info=info,
    fit=Kvpred_all,
    fitting_time=fitting_time,
    args=list(),
    data=data,
    data_daily=data_daily)
}


#### helpers ####

#' Make daily K estimates from preliminary daily K estimates
#' 
#' Called from metab_Kvpred().
#' 
#' @param data_all data.frame of the form \code{mm_data(local.time, K600, DO.obs,
#'   DO.sat, depth, temp.water, light)} and containing data for just one
#'   estimation-day (this may be >24 hours but only yields estimates for one
#'   24-hour period)
#' @param method the method to use
#' @return list of a data.frame of estimates and model diagnostics
#' @keywords internal
Kvpred_allply <- function(data_all, method) {
  stop("sorry; not yet implemented")
  # Return, reporting any results, warnings, and errors
  #   data.frame(local.date=NA,
  #              GPP=bayes.1d$mean.GPP.daily, GPP.sd=bayes.1d$sd.GPP.daily, 
  #              ER=bayes.1d$mean.ER.daily, ER.sd=bayes.1d$sd.ER.daily, 
  #              K600=bayes.1d$mean.K600.daily, K600.sd=bayes.1d$sd.K600.daily,
  #              warnings=paste0(warn_strs, collapse="; "), 
  #              errors=paste0(stop_strs, collapse="; "),
  #              stringsAsFactors=FALSE)
}


#### helpers to the helper ####



#### metab_Kvpred class ####

#' Interpolation model of daily K for metabolism
#' 
#' \code{metab_Kvpred} models use initial daily estimates of K, along with 
#' predictors such as Q (discharge) or U (velocity) or T (time) to leverage all
#' available data to reach better, less variable daily estimates of K
#' 
#' @exportClass metab_Kvpred
#' @family metab.model.classes
setClass(
  "metab_Kvpred", 
  contains="metab_model"
)

