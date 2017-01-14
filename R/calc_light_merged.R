#' Merge modeled and observed PAR into a single timeseries
#' 
#' Merge two time series of photosynthetically active radiation (PAR) for a 
#' series of date-times.
#' 
#' @inheritParams calc_light
#' @param PAR.obs a 2-column data.frame with columns solar.time and light, as in
#'   argument default, containing the full time series of observed light (should
#'   be at a lower temporal resolution than \code{PAR.mod})
#' @param PAR.mod a 2-column data.frame with columns solar.time and light, as in
#'   argument default, containing the full time series of modeled light
#' @param extrapolate logical. Should the output include periods of time where 
#'   only modeled light is available? Periods where only observed light is
#'   available (i.e., those not included in solar.time) are always excluded.
#' @param attach.units logical. Should the returned vector be a unitted object?
#' @import dplyr
#'   
#' @examples 
#' \dontrun{
#' library(mda.streams)
#' library(dplyr)
#' library(lubridate)
#' coords <- get_site_coords('nwis_08062500')
#' PAR.obs <- get_ts('par_calcSw', 'nwis_08062500') %>%
#'   mutate(solar.time = convert_UTC_to_solartime(DateTime, coords$lon)) %>%
#'     select(solar.time, light=par) %>%
#'     subset(solar.time >= as.POSIXct('2007-10-31', tz=tz(solar.time)) &
#'            solar.time < as.POSIXct('2007-11-03', tz=tz(solar.time)))
#' PAR.mod <- data.frame(solar.time=seq(as.POSIXct('2007-11-01', tz=tz(solar.time)),
#'     as.POSIXct('2007-11-04', tz=tz(solar.time)), by=as.difftime(10, units='mins'))) %>%
#'   mutate(light=calc_light(solar.time, latitude=coords$lat, longitude=coords$lon))
#' PAR.merged <- calc_light_merged(PAR.obs, PAR.mod$solar.time, latitude=coords$lat, longitude=coords$lon)
#' ggplot(bind_rows(mutate(PAR.obs, type='obs'), mutate(PAR.mod, type='mod'), mutate(PAR.merged, type='merged')), 
#'   aes(x=solar.time, y=light, color=type)) + geom_line() + geom_point()
#' }
#' @export
calc_light_merged <- function(
  PAR.obs=mm_data(solar.time, light),
  solar.time, 
  latitude, longitude, max.PAR=2326, coef.SW.to.PAR=formals(convert_SW_to_PAR)$coef,
  extrapolate=TRUE, attach.units=is.unitted(PAR.obs)) {
  
  # join the tses, noting which solar.times apply to which ts
  PAR.merged <- u(data.frame(solar.time, is.mod=TRUE, is.obs=NA)) %>%
    full_join(rename(PAR.obs, obs=light), by='solar.time') %>%
    {.[order(.$solar.time),]} %>%
    mutate(
      is.mod = ifelse(is.na(is.mod), FALSE, is.mod),
      is.obs = !is.na(obs))
  
  # if requested, remove rows with no supporting observations nearby in time
  if(!isTRUE(extrapolate)) {
    # calculate distance from nearest observation and remove 'mod' rows that are
    # too far from an obs
  }
  
  # create a 'merged' time series where modeled values are adjusted to flow
  # through observed values
  PAR.merged <- PAR.merged %>%
    mutate(
      mod = calc_light(solar.time, latitude, longitude, max.PAR, coef.SW.to.PAR),
      resid.prop = ifelse(mod==0 | (is.na(obs) & mod==0), 1, obs/mod))
  PAR.merged$resid.int <- approx(x=PAR.merged$solar.time, y=PAR.merged$resid.prop, xout=PAR.merged$solar.time, rule=2)$y # pipes fail here
  PAR.merged <- PAR.merged %>%
    mutate(
      #resid.calmed = ifelse(obs==0, 0, resid.int)
      merged = mod * resid.int) %>%
    subset(is.mod) %>%
    select(solar.time, light=merged)

}

