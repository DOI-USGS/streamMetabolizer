#' Calculate modeled light from solar.time
#'
#' Calculate photosynthetically active radiation (PAR) for a series of
#' date-times and site coordinates.
#'
#' @param solar.time mean solar time, as required for input to metabolism
#'   models. See \code{\link{mm_data}} and \code{\link{calc_solar_time}}.
#' @inheritParams calc_zenith_angle
#' @inheritParams convert_solartime_to_UTC
#' @param max.PAR numeric or unitted_numeric: the PAR (umol m^-2 s^-1) that each
#'   day should reach at peak light
#' @inheritParams calc_solar_insolation
#' @import dplyr
#' @importFrom lifecycle deprecated is_present
#' @examples
#' solar.time <- lubridate::force_tz(as.POSIXct('2016-09-27 12:00'), 'UTC')
#' calc_light(solar.time, 40, -120)
#' library(unitted)
#' calc_light(u(solar.time), u(40, 'degN'), u(-120, 'degE'), u(2315, 'umol m^-2 s^-1'))
#' @export
calc_light <- function(
  solar.time, latitude, longitude, max.PAR=u(2326, 'umol m^-2 s^-1'),
  attach.units=deprecated()) {

  # check units-related arguments. old default was attach.units=is.unitted(solar.time)
  if (lifecycle::is_present(attach.units)) {
    # only warn if it's TRUE
    if(isTRUE(attach.units)) unitted_deprecate_warn("calc_light(attach.units)")
  } else if (is.unitted(solar.time)) {
    unitted_deprecate_warn("calc_light(solar.time.unitted)")
    message('in calc_light, setting attach.units=TRUE because is.unitted(solar.time)')
    attach.units <- TRUE
  } else {
    attach.units <- FALSE
  }

  coef.SW.to.PAR <- formals(convert_SW_to_PAR)$coef # shouldn't really matter what is b/c we convert out and back
  app.solar.time <- v(solar.time) %>%
    convert_solartime_to_UTC(longitude=longitude, time.type='mean solar') %>%
    convert_UTC_to_solartime(longitude=longitude, time.type='apparent solar')
  sw <- calc_solar_insolation(
    app.solar.time, latitude=latitude,
    max.insolation=convert_PAR_to_SW(max.PAR, coef=1/coef.SW.to.PAR),
    format=c("degrees", "radians"), attach.units=attach.units)
  par <- convert_SW_to_PAR(sw, coef=coef.SW.to.PAR)

  par
}
