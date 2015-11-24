#' Filter unit or daily data by inclusive start & end dates
#' 
#' @param data either instantaneous/unit or daily data, having columns for local.time or local.date, respectively, to filter
#' @inheritParams predict_DO
#' @importFrom lubridate tz floor_date
#' @keywords internal
#' @examples 
#' tm <- Sys.time()
#' dt <- Sys.Date()
#' udat <- data.frame(local.time=tm + as.difftime(1:100, units='hours'), value=1:100)
#' ddat <- data.frame(local.date=dt + as.difftime(1:100, units='days'), value=1:100)
#' streamMetabolizer:::mm_filter_dates(udat, tm)
#' streamMetabolizer:::mm_filter_dates(udat, date_start=dt+as.difftime(1, units="days"))
mm_filter_dates <- function(data, date_start=NA, date_end=NA, day_start=4, day_end=27.99, date_format="%Y-%m-%d") {
  
  if(is.null(data) || nrow(data) == 0 || (is.character(data) && data == "generic metab_model class; no actual fit")) return(data)
  date_col <- unlist(sapply(c('local.time','local.date'), grep, names(data), fixed=TRUE, value=TRUE, USE.NAMES = FALSE))[1]
  data <- data[!is.na(data[[date_col]]),]
  # format dates
  tidy_date <- function(date) { 
    local.tz <- if(date_col=='local.time') tz(data$local.time) else 'UTC'
    loc.date <- as.Date(date, format=date_format, tz=local.tz)
    switch(
      date_col,
      'local.time' = as.POSIXct(format(loc.date, date_format), format=date_format, tz=local.tz),
      'local.date' = loc.date
    )
  }
  date_start <- tidy_date(date_start)
  date_end <- tidy_date(date_end)
  # filter
  switch(
    date_col,
    'local.time' = {
      # format dates
      data %>% dplyr::filter(
        if(!is.na(date_start)) local.time >= (date_start + as.difftime(day_start, units="hours")) else TRUE,
        if(!is.na(date_end)) local.time <= (date_end + as.difftime(day_end, units="hours")) else TRUE) 
    }, 
    'local.date' = {
      data %>% dplyr::filter(
        if(!is.na(date_start)) local.date >= date_start else TRUE,
        if(!is.na(date_end)) local.date <= date_end else TRUE)
    }
  )
}