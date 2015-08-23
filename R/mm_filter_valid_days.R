#' Remove entries in data
#' 
#' @param data
#' @inheritParams mm_is_valid_day
#' @return list of data and data_daily with same structure as inputs but with invalid days removed
#' @export
mm_filter_valid_days <- function(
  data, data_daily,
  tests=c('full_day', 'even_timesteps', 'complete_data'), 
  day_start=6, day_end=30, timestep_days=NA, need_complete=names(data)) {
  
  # filter the instantaneous data using mm_is_valid_day
  filter_fun <- function(data, data_daily, day_start, day_end, date, tests, timestep_days, need_complete) {
    stop_strs <- mm_is_valid_day(
      day=data, tests=tests,
      day_start=day_start, day_end=day_end,
      timestep_days=timestep_days, need_complete=need_complete)
    if(isTRUE(stop_strs)) {
      data 
      } else {
        removed <<- c(removed, list(data.frame(local.date=date, errors=paste0(stop_strs, collapse="; "))))
        NULL
      }
  }
  # actually run the filtering function over all days, recording days that were removed
  removed <- list()
  data_filtered <- mm_model_by_ply(
    model_fun=filter_fun, data=data, data_daily=data_daily, 
    day_start=day_start, day_end=day_end, tests=tests, timestep_days=timestep_days, need_complete=need_complete)
  removed <- do.call(rbind, removed)
  
  # filter the daily data to match & retur
  if(!missing(data_daily)) {
    data_daily_filtered <- data_daily[data_daily$local.date %in% unique(data_filtered$date),]
  } else {
    data_daily_filtered <- NULL
  }
  
  # return
  list(data=data_filtered, data_daily=data_daily_filtered, removed=removed)
}