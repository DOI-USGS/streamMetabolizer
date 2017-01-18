#' Calculate solar.time from local.time
#' 
#' Calculate the appropriate solar.time column for input to metab(). The input 
#' must be POSIXct clock time and should have the correct timezone information 
#' embedded in the object, whether the tz is UTC, local time with daylight
#' savings, or local standard time. The output is always mean solar time (not
#' apparent; see \code{convert_UTC_to_solaritme}).
#' 
#' @inheritParams convert_localtime_to_UTC
#' @inheritParams convert_UTC_to_solartime
#' @import dplyr
#' @examples 
#' local.time <- as.POSIXct('2016-05-27 12:00:00', tz='America/New_York')
#' solar.time <- calc_solar_time(local.time, longitude=3)
#' @export
calc_solar_time <- function(local.time, longitude) {
  
  # quick check to try to warn users away from using non-solar time
  stated_tz <- attr(local.time, 'tz')
  if(stated_tz == 'UTC') {
    implied_tz <- tryCatch(lookup_timezone(latitude=51.5, longitude), error=function(e) NULL) # london: 51.5N, 0W
    if(!is.null(implied_tz) && !implied_tz$std_offset %in% -1:1) # give some buffer. UTC longitudes vary with latitude
      warning('found non-UTC longitude for UTC timezone. Are you sure you passed in a local time?')
  }
  
  local.time %>%
    convert_localtime_to_UTC() %>%
    convert_UTC_to_solartime(longitude=longitude, time.type='mean solar')
}
