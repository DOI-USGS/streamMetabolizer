#' Limit to a specific time range on each date
#' 
#' Within each date (as labeled by the 'date' column of \code{data}, select the 
#' values of solar.time that are within the time range specified by day_start
#' and day_end). This function only removes rows and cannot add them; to add
#' overlap starting from a continuous time series, see
#' \code{\link{mm_model_by_ply}}.
#' 
#' @param data a data.frame containing date and solar.time columns (POSIXct)
#' @param day_start the start time of each day, inclusive, in hours
#' @param day_end the end time of each day, exclusive, in hours
#' @import dplyr
#' @keywords internal
mm_filter_hours <- function(data, day_start, day_end) {
  data %>%
    mutate(d0 = with_tz(as.POSIXct(data$date), 'UTC'),
           ds = d0 + as.difftime(day_start, units='hours'),
           de = d0 + as.difftime(day_end, units='hours')) %>%
    filter(solar.time >= ds, solar.time < de) %>%
    select(-d0, -ds, -de)
}
