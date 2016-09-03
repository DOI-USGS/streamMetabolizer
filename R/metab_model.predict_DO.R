#' @include metab_model-class.R
NULL

#' Make dissolved oxygen predictions from a fitted metab_model.
#' 
#' Makes fine-scale predictions of dissolved oxygen using fitted coefficients, 
#' etc. from the metabolism model.
#' 
#' @inheritParams predict_DO
#' @return A data.frame of predictions, as for the generic 
#'   \code{\link{predict_DO}}.
#' @import dplyr
#' @importFrom lubridate tz
#' @importFrom unitted u v get_units
#' @export
#' @family predict_DO
predict_DO.metab_model <- function(metab_model, date_start=NA, date_end=NA, 
                                   day_start=get_specs(metab_model)$day_start, day_end=get_specs(metab_model)$day_end,
                                   ..., attach.units=FALSE, use_saved=TRUE) {
  
  # pull args from the model
  specs <- get_specs(metab_model)
  
  # get the input data; filter if requested
  data <- get_data(metab_model) %>% 
    v() %>% # units are discarded here. someday we'll get them worked through the whole system.
    mm_filter_dates(date_start=date_start, date_end=date_end, day_start=day_start, day_end=day_end)
  
  # if allowed and available, use previously stored values for DO.mod rather than re-predicting them now
  if(isTRUE(use_saved) && !is.null(data) && "DO.mod" %in% names(data)) {
    return(data)
  }
  
  # get the metabolism estimates; filter as we did for data
  metab_ests <- get_params(metab_model, date_start=date_start, date_end=date_end, uncertainty='none', messages=FALSE)
  
  # if we lack the data or params to predict, return NULL now
  if(is.null(data) || nrow(data) == 0 || is.null(metab_ests)) 
    return(NULL)
  
  # re-process the input data with the metabolism estimates to predict DO
  preds <- mm_model_by_ply(
    mm_predict_DO_1ply, data=data, data_daily=metab_ests, # for mm_model_by_ply
    day_start=day_start, day_end=day_end, day_tests=c(), timestep_days=FALSE, # for mm_model_by_ply
    model_name=specs$model_name) %>% # for mm_predict_DO_1ply
    mm_filter_dates(date_start=date_start, date_end=date_end) # trim off the extra
  
  if(attach.units) preds <- u(preds, get_units(mm_data())[gsub('DO.mod','DO.obs',names(preds))])
  preds
}