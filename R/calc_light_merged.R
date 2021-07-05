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
#' @param max.PAR the maximum PAR, as in calc_light. if NA, this function does
#'   its best to guess a max.PAR that will make modeled light pretty similar to
#'   cloud-free days of observed light
#' @param max.gap difftime or NA. If difftime, the maximum gap between a light
#'   observation and a time point in solar.time, beyond which no value will be
#'   given for light at that solar.time. If NA, all values will be modeled, even
#'   if they are many days away from a light observation.
#' @param attach.units (deprecated, effectively FALSE in future) logical. Should
#'   the returned vector be a unitted object?
#' @import dplyr
#' @importFrom unitted u v
#' @importFrom lifecycle deprecated is_present
#' @export
#' @examples
#' \dontrun{
#' library(dplyr)
#' library(ggplot2)
#' library(unitted)
#' timebounds <- as.POSIXct(c('2008-03-12 00:00', '2008-03-12 23:59'), tz='UTC')
#' coords <- list(lat=32.4, lon=-96.5)
#' PAR.obs <- tibble::tibble(
#'   solar.time=seq(timebounds[1], timebounds[2], by=as.difftime(3, units='hours')),
#'   light=c(0, 0, 85.9, 1160.5, 1539.0, 933.9, 0, 0)
#' ) %>% as.data.frame()
#' PAR.mod <- tibble::tibble(
#'   solar.time=seq(timebounds[1], timebounds[2], by=as.difftime(0.25, units='hours')),
#'   light=calc_light(solar.time, latitude=coords$lat, longitude=coords$lon)
#' ) %>% as.data.frame()
#' PAR.merged <- calc_light_merged(PAR.obs, PAR.mod$solar.time,
#'   latitude=coords$lat, longitude=coords$lon, max.gap=as.difftime(20, units='hours'))
#' ggplot(bind_rows(mutate(v(PAR.obs), type='obs'), mutate(v(PAR.mod), type='mod'),
#'                  mutate(v(PAR.merged), type='merged')) %>%
#'        mutate(type=ordered(type, levels=c('obs','mod','merged'))),
#'   aes(x=solar.time, y=light, color=type)) + geom_line() + geom_point() + theme_bw()
#' }
#' @export
calc_light_merged <- function(
  PAR.obs=mm_data(solar.time, light),
  solar.time, latitude, longitude, max.PAR=NA,
  max.gap=as.difftime(3, units="hours"), attach.units=deprecated()) {

  # check units-related arguments. old default was attach.units=is.unitted(PAR.obs)
  if (lifecycle::is_present(attach.units)) {
    unitted_deprecate_warn("calc_light_merged(attach.units)")
  } else if (is.unitted(PAR.obs)) {
    unitted_deprecate_warn("calc_light_merged(PAR.obs.unitted)")
    message('in calc_light_merged, setting attach.units=TRUE because is.unitted(PAR.obs)')
    attach.units <- TRUE
  } else {
    attach.units <- FALSE
  }

  # ensure units are correct and present within this function
  arg_units <- list(
    PAR.obs=get_units(mm_data(solar.time, light)),
    solar.time='',
    latitude='degN',
    longitude='degE')
  for(argname in names(arg_units)) {
    arg <- get(argname)
    unit <- arg_units[[argname]]
    if(is.unitted(arg)) {
      stopifnot(all(get_units(arg) == unit))
    } else {
      assign(argname, u(arg, unit))
    }
  }

  . <- is.mod <- obs <- mod <- resid.abs.int <- resid.prop.int <- merged <- '.dplyr.var'

  # set smart default for max.PAR to make the ts pretty
  if(is.na(max.PAR)) {
    # figure out what & when the max insolation is in the input data
    date.max.obs <- PAR.obs[which.max(PAR.obs$light), 'solar.time'] # this is mean solar time; will also need apparent solar
    max.obs <- PAR.obs[which.max(PAR.obs$light), 'light']

    # given some known value of old.max.PAR, what would the model predict for
    # solar insolation at the date when obs is at its maximum?
    old.max.PAR <- u(2000, 'umol m^-2 s^-1') # doesn't matter what it is
    mod.at.max.obs <- calc_light(date.max.obs, latitude, longitude, old.max.PAR)

    # find the max.PAR such that modeled and observed PAR match closely at the
    # time+date of observed max par. assumes linear relationship between PAR and
    # Sw as currently in LakeMetabolizer::sw.to.par.base as of 1/17/2017, and
    # hopes that max.obs is on a representatively sunny day
    max.PAR <- old.max.PAR * (max.obs / mod.at.max.obs)
  }

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
  # through observed values.

  # compute the modeled values and then the residuals as both differences and
  # proportions, dealing with 0 specially
  PAR.merged <- PAR.merged %>%
    mutate(
      mod = calc_light(solar.time, latitude, longitude, max.PAR),
      resid.abs = obs - mod,
      resid.prop = ifelse(mod==u(0, 'umol m^-2 s^-1'), NA, obs / mod))

  # interpolate the residuals to match up with every modeled light value. pipes
  # fail with this approx call, so use boring notation
  PAR.obsonly <- PAR.merged[!is.na(PAR.merged$obs), ]
  PAR.merged$resid.abs.int <- approx(x=PAR.obsonly$solar.time, y=v(PAR.obsonly$resid.abs), xout=v(PAR.merged$solar.time), rule=2)$y
  PAR.merged$resid.prop.int <- approx(x=PAR.obsonly$solar.time, y=PAR.obsonly$resid.prop, xout=PAR.merged$solar.time, rule=2)$y

  # do the correction from mod scale to obs scale
  PAR.merged <- PAR.merged %>%
    mutate(
      merged = u(ifelse(resid.prop.int <= 1, mod * resid.prop.int, pmax(0, v(mod) + resid.abs.int)), get_units(mod)))

  # collect just the rows and cols we want
  PAR.merged <- PAR.merged %>%
    subset(is.mod) %>%
    select(solar.time, light=merged)

  if(!attach.units) PAR.merged <- v(PAR.merged)

  # return
  PAR.merged
}
