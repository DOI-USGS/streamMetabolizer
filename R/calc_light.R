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
#' @param coef.SW.to.PAR conversion coefficient; see \code{code} in 
#'   \code{\link{convert_SW_to_PAR}}
#' @inheritParams calc_solar_insolation
#' @import dplyr
#' @examples 
#' solar.time <- lubridate::force_tz(as.POSIXct('2016-09-27 12:00'), 'UTC')
#' calc_light(solar.time, 40, -120)
#' @export
calc_light <- function(solar.time, latitude, longitude, max.PAR=2326, 
                       coef.SW.to.PAR=formals(convert_SW_to_PAR)$coef, attach.units=is.unitted(solar.time)) {
  
  app.solar.time <- solar.time %>%
    convert_solartime_to_UTC(longitude=longitude, time.type='mean solar') %>%
    convert_UTC_to_solartime(longitude=longitude, time.type='apparent solar')
  sw <- calc_solar_insolation(
    app.solar.time, latitude=latitude, 
    max.insolation=convert_PAR_to_SW(max.PAR, coef=1/coef.SW.to.PAR), 
    format=c("degrees", "radians"), attach.units=attach.units)
  par <- convert_SW_to_PAR(sw, coef=coef.SW.to.PAR)
  
  par
}
