#' Filter unit or daily data by inclusive start & end dates
#' 
#' @param data either instantaneous/unit or daily data, having columns for
#'   solar.time or date, respectively, to filter
#' @inheritParams predict_DO
#' @importFrom lubridate tz floor_date
#' @keywords internal
#' @examples 
#' tm <- as.POSIXct("2017-10-02 00:00:00 UTC")
#' dt <- as.Date("2017-10-02")
#' udat <- data.frame(solar.time=tm + as.difftime(1:100, units='hours'), value=1:100)
#' udat1 <- streamMetabolizer:::mm_filter_dates(udat)
#' udat2 <- streamMetabolizer:::mm_filter_dates(udat, date_start=dt, date_end=dt)
#' udat3 <- streamMetabolizer:::mm_filter_dates(udat, date_start=dt, date_end=dt, 
#'   day_start=12, day_end=14)
#' c(nrow(udat), nrow(udat1), nrow(udat2), nrow(udat3))
#' ddat <- data.frame(date=dt + as.difftime(1:100, units='days'), value=1:100)
#' ddat1 <- streamMetabolizer:::mm_filter_dates(ddat)
#' ddat2 <- streamMetabolizer:::mm_filter_dates(ddat, date_start=dt+10, date_end=dt+20)
#' c(nrow(ddat), nrow(ddat1), nrow(ddat2))
mm_filter_dates <- function(data, date_start=NA, date_end=NA, day_start=4, day_end=28, date_format="%Y-%m-%d") {
  
  if(is.null(data) || nrow(data) == 0 || (is.character(data) && data == "generic metab_model class; no actual fit")) return(data)
  date_col <- unlist(sapply(c('solar.time','date'), grep, names(data), fixed=TRUE, value=TRUE, USE.NAMES = FALSE))[1]
  data <- data[!is.na(data[[date_col]]),]
  # format dates
  tidy_date <- function(date) { 
    local.tz <- if(date_col=='solar.time') tz(data$solar.time) else 'UTC'
    loc.date <- as.Date(date, format=date_format, tz=local.tz)
    switch(
      date_col,
      'solar.time' = as.POSIXct(format(loc.date, date_format), format=date_format, tz=local.tz),
      'date' = loc.date
    )
  }
  date_start <- tidy_date(date_start)
  date_end <- tidy_date(date_end)
  # filter
  switch(
    date_col,
    'solar.time' = {
      # format dates
      data %>% dplyr::filter(
        if(!is.na(date_start)) solar.time >= (date_start + as.difftime(day_start, units="hours")) else TRUE,
        if(!is.na(date_end)) solar.time <= (date_end + as.difftime(day_end, units="hours")) else TRUE) 
    }, 
    'date' = {
      data %>% dplyr::filter(
        if(!is.na(date_start)) date >= date_start else TRUE,
        if(!is.na(date_end)) date <= date_end else TRUE)
    }
  )
}