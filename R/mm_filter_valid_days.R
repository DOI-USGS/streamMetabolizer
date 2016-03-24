#' Remove entries in data
#' 
#' Filter out any data rows that don't pass the specified tests for completeness
#' and regularity
#' 
#' @param data data.frame of instantaneous observations, to be filtered to only 
#'   those points on days that pass the specified tests in mm_is_valid_day
#' @param data_daily data.frame of daily estimates/statistics, to be filtered in
#'   accordance with the filtering of data
#' @inheritParams mm_model_by_ply
#' @return list of data and data_daily with same structure as inputs but with 
#'   invalid days removed, plus a third data.frame of dates that were removed
#' @examples
#' dat <- data_metab(res='30', num_days='10', flaws='missing middle')
#' datfilt <- mm_filter_valid_days(dat)
#' datfilt$removed
#' c(nrow(dat), nrow(datfilt$data))
#' @export
mm_filter_valid_days <- function(
  data, data_daily=NULL, # redefine from metab
  day_start=4, day_end=27.99, day_tests=c('full_day', 'even_timesteps', 'complete_data'), timestep_days=TRUE # inheritParams mm_model_by_ply
) {
  
  # function to filter the instantaneous data using validity, record dates that
  # are removed and the reasons for removal in a parent variable named removed
  filter_fun <- function(data_ply, ply_date, ply_validity, ...) {
    if(isTRUE(ply_validity)) { # day is valid
      data_ply
    } else {
      removed <<- c(removed, list(data.frame(date=ply_date, errors=paste0(ply_validity, collapse="; "))))
      NULL
    }
  }

  # run the filtering function over all days, recording days that were removed
  removed <- list()
  data_filtered <- mm_model_by_ply(
    model_fun=filter_fun, data=data, data_daily=data_daily, 
    day_start=day_start, day_end=day_end, day_tests=day_tests)
  removed <- do.call(rbind, removed) # removed was populated by <<- calls within filter_fun
  
  # filter the daily data to match & return
  if(!is.null(data_daily)) {
    daily_unmatched <- as.Date(setdiff(
      as.character(data_daily$date), 
      c(unique(format(data$solar.time, "%Y-%m-%d")), as.character(removed$date))))
    daily_removed <- data.frame(
      date=daily_unmatched, 
      errors=rep("date in data_daily but not data", length(daily_unmatched)))
    removed <- rbind(removed, daily_removed)
    removed <- removed[order(removed$date),]
    rownames(removed) <- NULL
    data_daily_filtered <- data_daily[data_daily$date %in% unique(data_filtered$date),]
  } else {
    data_daily_filtered <- NULL
  }
  
  # return
  list(data=data_filtered, data_daily=data_daily_filtered, removed=removed)
}
