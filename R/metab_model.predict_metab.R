#' @include metab_model-class.R
NULL

#' Make metabolism predictions from a fitted metab_model.
#' 
#' Makes daily predictions of GPP, ER, and K600 with upper and lower bounds
#' reflecting a 95% CI.
#' 
#' @inheritParams predict_metab
#' @param attach.units logical. Should units be attached to the output?
#' @return A data.frame of predictions, as for the generic 
#'   \code{\link{predict_metab}}.
#' @importFrom stats qnorm setNames
#' @import dplyr
#' @importFrom unitted u v get_units
#' @export
#' @family predict_metab
predict_metab.metab_model <- function(metab_model, date_start=NA, date_end=NA, 
                                      day_start=get_specs(metab_model)$day_start, day_end=day_start+24,
                                      ..., attach.units=FALSE, use_saved=TRUE) {
  
  if(isTRUE(use_saved) && !is.null(metab_model@metab_daily)) {
    # if allowed and available, use previously stored values rather than
    # re-predicting them now
    preds <- metab_model@metab_daily
    
  } else {
    # otherwise predict them now
    
    # get the metabolism parameters; filter if requested
    metab_ests <- get_params(metab_model, date_start=date_start, date_end=date_end, uncertainty='sd', messages=FALSE)
    
    # return now if there's no estimation to be done
    if(is.null(metab_ests)) return(NULL)
    
    # pull args from the model
    specs <- get_specs(metab_model)
    if(missing(day_end) && specs$day_end != day_end)
      warning("default day_end (day_start + 24 = ", day_end, ") overrides get_specs(mm)$day_end (", specs$day_end, ")")
    if((day_end - day_start) != 24) {
      warning("predictions are means of non-24-hour periods because (day_end - day_start) != 24")
    }
    
    # get the instantaneous data, including DO.mod; filter if requested
    data <- predict_DO(
      metab_model, date_start=date_start, date_end=date_end, day_start=day_start, day_end=day_end,
      attach.units=FALSE, use_saved=TRUE)
    
    # re-process the input data with the metabolism estimates to predict daily mean metabolism
    preds <- mm_model_by_ply(
      mm_predict_metab_1ply, data=data, data_daily=metab_ests, # for mm_model_by_ply
      day_start=day_start, day_end=day_end, day_tests=c("full day", "even_timesteps", "complete_data"), timestep_days=FALSE, # for mm_model_by_ply
      model_name=specs$model_name) # for mm_predict_DO_1ply
    
    # attach warnings and errors if available
    preds.cols <- append(names(preds), 'msgs.fit', after=which(names(preds) == 'warnings') - 1)
    fit <- get_fit(metab_model)
    if(!is.null(fit) && all(exists(c('date','warnings','errors'), fit))) {
      messages <- fit %>%
        select(date, warnings, errors) %>%
        compress_msgs('msgs.fit')
      preds <- full_join(preds, messages, by='date', copy=TRUE)
    } else {
      preds <- mutate(preds, msgs.fit=NA)
    }
    preds <- select_(preds, .dots=preds.cols)
  }
  
  if(attach.units) {
    pred.units <- get_units(mm_data())[sapply(names(preds), function(x) strsplit(x, '\\.')[[1]][1], USE.NAMES=FALSE)]
    preds <- u(preds, pred.units)
  }
  preds
}