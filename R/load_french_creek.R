#' Load a short dataset from French Creek
#' 
#' @importFrom unitted u
#' @param attach.units logical. Should units be attached to the data.frame?
#' @return a data.frame, unitted if attach.units==TRUE
load_french_creek <- function(attach.units=TRUE) {
  # load the file
  file.name <- system.file("extdata", "french_creek_data.csv", package="streamMetabolizer")
  french <- read.csv(file.name, stringsAsFactors=FALSE, skip=2)
  french_creek_units <- sapply(read.csv(file.name, stringsAsFactors=FALSE, nrow=1), function(col) as.character(col))
  french <- if(attach.units) u(french, french_creek_units) else french
  
  # pre-calculate DO at sat
  french$DO.sat <- calc_DO_at_sat(temp.water=french$temp.water, pressure.air=u(1000, "mb"))
  
  # fix up the time & light assuming this is the french creek near buffalo
  french$date.time <- as.POSIXct(strptime(french$date.time, format="%m/%d/%Y %H:%M"), tz="MST")
  french$utc.time <- convert_localtime_to_GMT(french$date.time)
  french$local.time <- force_tz(convert_GMT_to_localtime(french$utc.time, longitude=-106.753099, latitude=44.362594), "UTC")
  french$solar.time <- convert_GMT_to_solartime(french$utc.time, longitude=-106.753099, time.type='apparent solar')
  french$light <- convert_SW_to_PAR(calc_solar_insolation(solar.time=v(french$solar.time), latitude=44.362594, attach.units=TRUE))
  
  # return w/ proper units & columns
  french <- if(attach.units) french else v(french)
  french[c("local.time","DO.obs","DO.sat","depth","temp.water","light")]
}