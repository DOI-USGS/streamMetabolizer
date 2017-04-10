#' @include metab_model-class.R
NULL


#' @describeIn predict_metab This implementation is shared by many model types
#' @export
#' @importFrom stats qnorm setNames
#' @import dplyr
#' @importFrom unitted u v get_units
predict_metab.metab_model <- function(
  metab_model, date_start=NA, date_end=NA, 
  day_start=get_specs(metab_model)$day_start, day_end=min(day_start+24, get_specs(metab_model)$day_end),
  ..., attach.units=FALSE, use_saved=TRUE) {
  
  if(isTRUE(use_saved) && !is.null(metab_model@metab_daily)) {
    # if allowed and available, use previously stored values rather than
    # re-predicting them now
    if((!missing(day_start) && day_start != get_specs(metab_model)$day_start) || 
       (!missing(day_end) && day_end != min(day_start+24, get_specs(metab_model)$day_end))) 
      warning("using saved daily metabolism values and so ignoring new day_start and/or day_end values")
    preds <- metab_model@metab_daily %>%
      mm_filter_dates(date_start=date_start, date_end=date_end)
    
  } else {
    # otherwise predict them now
    
    # get the metabolism parameters; filter if requested
    metab_ests <- get_params(metab_model, date_start=date_start, date_end=date_end, uncertainty='sd', messages=FALSE)
    
    # return now if there's no estimation to be done
    if(is.null(metab_ests)) return(NULL)
    
    # pull args from the model
    specs <- get_specs(metab_model)
    
    # consider the appropriateness of day_start and day_end
    if(specs$day_start != day_start || specs$day_end != day_end) {
      if(day_start < specs$day_start) stop("day_start may not be earlier than the day_start stored in metab_model@specs")
      if(day_end > specs$day_end) stop("day_end may not be later than the day_end stored in metab_model@specs")
      message(paste(
        "daily metabolism predictions are for the period from", day_start, "to", day_end, "hours on each date\n",
        "(differs from the model-fitting range of", specs$day_start, "to", specs$day_end, "hours)"))
    }
    if(day_end - day_start > 24) {
      # give error mostly because mm_model_by_ply can't currently handle this, 
      # but also because it's hard to interpret a mean metabolism for a period 
      # other than 24 hours, given that light and temperature and DO deficits
      # all vary systematically with time of day
      stop("day_end - day_start must not exceed 24 hours for metabolism prediction")
    } else if((day_end - day_start) < 24) {
      if(mm_parse_name(specs$model_name)$type != 'night') {
        # for most models, stop, because the GPP_fun=='linlight' gives the wrong
        # answers for <24-hour periods (normalizes by a different light value
        # such that the period-specific average is always either the daily
        # average or 0, regardless of which period is selected). and who knows
        # what new GPP or ER functions might break similarly, so stay on the
        # safe side by requiring what most people will want anyway (a
        # 24-hour-period prediction)
        stop("day_end - day_start < 24 hours; this is unacceptable except for metab_night")
      } else {
        # but metab_night may not have 24-hour periods available, and 
        # mm_model_by_ply CAN handle <24-hour periods, so make an exception for 
        # metab_night. and actually don't even give a warning because it's
        # really fine - metab_night will give the correct answers (within its
        # abilities) regardless of the period of time specified here
        message("for metab_night, GPP estimates are 0 because they're for nighttime only ")
      }
    }
    
    # get the instantaneous data, including DO.mod (which we need for predicting
    # D); filter to a 24-hour period
    data <- predict_DO(metab_model, date_start=date_start, date_end=date_end, attach.units=FALSE, use_saved=TRUE) %>%
      mm_filter_hours(day_start=day_start, day_end=day_end)
    
    # re-process the input data with the metabolism estimates to predict daily
    # mean metabolism. no need to test days again unless we didn't do it during
    # fitting (which only occurs for sim)
    day_tests_pred <- if(mm_parse_name(specs$model_name)$type == 'sim') specs$day_tests else c() 
    req_ts_pred <- if(mm_parse_name(specs$model_name)$type == 'sim') specs$required_timestep else NA 
    preds <- mm_model_by_ply(
      mm_predict_metab_1ply, data=data, data_daily=metab_ests, # for mm_model_by_ply
      day_start=day_start, day_end=day_end, day_tests=day_tests_pred, required_timestep=req_ts_pred, timestep_days=FALSE, # for mm_model_by_ply
      model_name=specs$model_name) # for mm_predict_DO_1ply
    
    # attach warnings and errors if available
    warnings <- errors <- '.dplyr.var'
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
