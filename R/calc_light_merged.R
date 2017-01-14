#' Merge modeled and observed PAR into a single timeseries
#' 
#' Merge two time series (one observed, one modeled) of photosynthetically 
#' active radiation (PAR) for a series of date-times. You can also think about 
#' this as a way to smoothly interpolate a time series of observations.
#' 
#' @param PAR.obs a 2-column data.frame with columns solar.time and light, as in
#'   argument default, containing the full time series of observed light (should
#'   be at a lower temporal resolution than \code{PAR.mod})
#' @param solar.time a vector of mean solar times for which the light should be 
#'   modeled and merged with the values in PAR.obs
#' @inheritParams calc_light
#' @param max.gap difftime or NA. If difftime, the maximum gap between a light 
#'   observation and a time point in solar.time, beyond which no value will be 
#'   given for light at that solar.time. If NA, all values will be modeled, even
#'   if they are many days away from a light observation.
#' @param attach.units logical. Should the returned vector be a unitted object?
#' @import dplyr
#' @export
#' @examples 
#' \dontrun{
#' library(mda.streams)
#' library(dplyr)
#' library(lubridate)
#' library(ggplot2)
#' library(unitted)
#' coords <- get_site_coords('nwis_08062500')
#' PAR.obs <- get_ts('par_calcSw', 'nwis_08062500') %>%
#'   mutate(solar.time = convert_UTC_to_solartime(DateTime, coords$lon)) %>%
#'     select(solar.time, light=par) %>%
#'     subset(solar.time %within% interval(ymd("2007-10-31"), ymd("2007-11-03"), tz=tz(solar.time)) &
#'           !(solar.time %within% interval(as.POSIXct("2007-11-01 10:00", tz=tz(solar.time)), 
#'                                          as.POSIXct("2007-11-01 12:00", tz=tz(solar.time)))))
#' PAR.mod <- data.frame(solar.time=seq(as.POSIXct('2007-11-01', tz=tz(PAR.obs$solar.time)),
#'     as.POSIXct('2007-11-04', tz=tz(PAR.obs$solar.time)), by=as.difftime(10, units='mins'))) %>%
#'   mutate(light=calc_light(solar.time, latitude=coords$lat, longitude=coords$lon))
#' PAR.merged <- calc_light_merged(PAR.obs, PAR.mod$solar.time, 
#'   latitude=coords$lat, longitude=coords$lon, max.gap=as.difftime(12, units='hours'))
#' ggplot(bind_rows(mutate(v(PAR.obs), type='obs'), mutate(v(PAR.mod), type='mod'), 
#'                  mutate(v(PAR.merged), type='merged')),
#'   aes(x=solar.time, y=light, color=type)) + geom_line() + geom_point()
#' ggplot(bind_rows(mutate(v(PAR.obs), type='obs'), mutate(v(PAR.mod), type='mod'), 
#'                  mutate(v(PAR.merged), type='merged')) %>% 
#'        filter(solar.time %within% interval(ymd("2007-11-01"), ymd("2007-11-02"), 
#'                                            tz=tz(PAR.obs$solar.time))), 
#'   aes(x=solar.time, y=light, color=type)) + geom_line() + geom_point()
#' }
#' @export
calc_light_merged <- function(
  PAR.obs=mm_data(solar.time, light), 
  solar.time, latitude, longitude, max.PAR=max(PAR.obs$light),
  max.gap=as.difftime(3, units="hours"), attach.units=is.unitted(PAR.obs)) {
  
  . <- is.mod <- obs <- mod <- resid.int <- merged <- '.dplyr.var'
  
  # join the tses, noting which solar.times apply to which ts
  PAR.merged <- u(data.frame(solar.time, is.mod=TRUE, is.obs=NA)) %>%
    full_join(rename(PAR.obs, obs=light), by='solar.time') %>%
    {.[order(.$solar.time),]} %>%
    mutate(
      is.mod = ifelse(is.na(is.mod), FALSE, is.mod),
      is.obs = !is.na(obs))
  
  # if requested, remove rows with no supporting observations nearby in time
  if(!is.na(max.gap)) {
    # calculate distance from nearest observation and remove 'mod' rows that are
    # too far from an obs
    # figure out what the time gap is between a given x datetime and the
    # nearest y datetime
    x_date_num <- as.numeric(PAR.merged$solar.time)
    y_date_num <- as.numeric(PAR.merged$solar.time[PAR.merged$is.obs])
    prev_y <- approx(x=y_date_num, y=y_date_num, xout=x_date_num, method="constant", f=0, rule=2)
    next_y <- approx(x=y_date_num, y=y_date_num, xout=x_date_num, method="constant", f=1, rule=2) 
    min_gap <- pmin(abs(prev_y$y-prev_y$x), abs(next_y$y-next_y$x))
    PAR.merged <- PAR.merged[min_gap <= as.numeric(max.gap, units="secs"), ]
  }
  
  # create a 'merged' time series where modeled values are adjusted to flow
  # through observed values
  PAR.merged <- PAR.merged %>%
    mutate(
      # model light
      mod = calc_light(solar.time, latitude, longitude, max.PAR),
      # compute the residuals (as proportions, dealing with 0 specially)
      resid.prop = ifelse(mod==u(0, 'umol m^-2 s^-1') | (is.na(obs) & mod==u(0, 'umol m^-2 s^-1')), 1, obs/mod))
  # interpolate the residuals to match up with every modeled light value. pipes fail with this approx call, so use boring notation
  PAR.merged$resid.int <- approx(x=PAR.merged$solar.time, y=PAR.merged$resid.prop, xout=PAR.merged$solar.time, rule=2)$y
  PAR.merged <- PAR.merged %>%
    # do the correction from mod scale to obs scale
    mutate(merged = mod * resid.int) %>%
    # collect just the rows and cols we want
    subset(is.mod) %>%
    select(solar.time, light=merged)
  
  # return
  PAR.merged
}

