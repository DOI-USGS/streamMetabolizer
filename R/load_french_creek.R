#' Load a short dataset from French Creek
#'
#' @import dplyr
#' @importFrom unitted u
#' @importFrom lifecycle deprecated
#' @importFrom utils read.csv
#' @importFrom lubridate with_tz
#' @param attach.units logical. Should units be attached to the data.frame?
#' @return a data.frame, unitted if attach.units==TRUE
load_french_creek <- function(attach.units=deprecated()) {
  # check arguments
  if (lifecycle::is_present(attach.units)) {
    # Signal the deprecation to the user
    deprecate_warn(
      "0.12.0", "streamMetabolizer::load_french_creek(attach.units)",
      details = "In the future, streamMetabolizer will stop using the unitted package entirely.")
  }

  # load the file
  file.name <- system.file("extdata", "french.csv", package="streamMetabolizer") # data from French Creek, Hotchkiss and Hall, In press, Ecology
  french <- read.csv(file.name, stringsAsFactors=FALSE, header=TRUE)

  . <- oxy <- temp <- station <- solar.time <- '.dplyr.var'

  # subset to 'low' site, remove NA oxys (1658) and remaining duplicates (n=1),  and ensure order by date
  french <- french %>%
    filter(station == 'low') %>% # subset to data from only one station (also, low has the cleanest data) (already subset in extdata)
    filter(!is.na(oxy)) %>%
    distinct() %>%
    arrange(solar.time)

  # rename DO.obs, temp.water
  french <- rename(french, DO.obs=oxy, temp.water=temp)

  # datetime
  tz_french <- lubridate::tz(convert_UTC_to_localtime(as.POSIXct("2012-09-10 00:00:00", tz="UTC"), latitude=41.33, longitude=-106.3, time.type="standard"))
  french$local.time <- lubridate::with_tz(as.POSIXct(paste(french$date, french$time), format="%m/%d/%Y %H:%M:%S", tz="America/Denver"), tz_french) # original is in MDT
  french$utc.time <- convert_localtime_to_UTC(french$local.time)
  french$solar.time <- convert_UTC_to_solartime(french$utc.time, longitude=-106.3, time.type='mean solar')

  # DO at sat
  pressure_air_mb <- 523 * 1.33322368 # 523 mmHg * 1.33322368 mb mmHg^-1 -> 697.27 mb, applicable at ~10000 ft
  french$DO.sat <- calc_DO_sat(temp.water=french$temp.water, pressure.air=pressure_air_mb)

  # depth
  french$depth <- 0.16 # meters

  # light
  french$app.solar.time <- convert_UTC_to_solartime(french$utc.time, longitude=-106.3, time.type='apparent solar')
  french$light <- convert_PAR_to_SW(2326) %>%
    calc_solar_insolation(app.solar.time=french$app.solar.time, latitude=41.33, max.insolation=.) %>%
    convert_SW_to_PAR()

  # set columns
  french <- french %>%
    select(c(solar.time, DO.obs, DO.sat, depth, temp.water, light))

  # add units if requested
  if (lifecycle::is_present(attach.units) && isTRUE(attach.units)) {
    french <- unitted::u(
      french,
      c(solar.time="", DO.obs="mgO2 L^-1", DO.sat="mgO2 L^-1", depth="m", temp.water="degC", light="umol m^-2 s^-1"))
  }

  return(french)
}
