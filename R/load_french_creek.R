#' Load a short dataset from French Creek
#' 
#' @import dplyr
#' @importFrom unitted u rename_.unitted_data.frame
#' @importFrom utils read.csv
#' @importFrom lubridate with_tz
#' @param attach.units logical. Should units be attached to the data.frame?
#' @return a data.frame, unitted if attach.units==TRUE
load_french_creek <- function(attach.units=TRUE) {
  # load the file
  file.name <- system.file("extdata", "french.csv", package="streamMetabolizer") # data from French Creek, Hotchkiss and Hall, In press, Ecology
  french <- read.csv(file.name, stringsAsFactors=FALSE, header=TRUE) 
  french <- french[french$station=="low",] # subset to data from only one station (also, low has the cleanest data) (already subset in extdata)
  french_creek_units <- c(station=NA, siteno=NA, sonde=NA, date=NA, time=NA, temp='degC', oxy='mgO2 L^-1')
  french <- unitted::u(french, french_creek_units)
  
  # remove NA oxys (1658) and remaining duplicates (n=1)
  french <- unique(french[!is.na(french$oxy),])
  
  # rename DO.obs, temp.water
  french <- unitted::rename_.unitted_data.frame(french, DO.obs='oxy', temp.water='temp')
  
  # datetime
  tz_french <- lubridate::tz(convert_UTC_to_localtime(as.POSIXct("2012-09-10 00:00:00", tz="UTC"), latitude=41.33, longitude=-106.3, time.type="standard"))
  french$local.time <- lubridate::with_tz(as.POSIXct(paste(french$date, french$time), format="%m/%d/%Y %H:%M:%S", tz="America/Denver"), tz_french) # original is in MDT
  french$utc.time <- convert_localtime_to_UTC(french$local.time)
  french$solar.time <- convert_UTC_to_solartime(french$utc.time, longitude=-106.3, time.type='mean solar')
  
  # DO at sat
  french$DO.sat <- calc_DO_at_sat(temp.water=french$temp.water, pressure.air=unitted::u(523, "mmHg")*unitted::u(1.33322368, "mb mmHg^-1")) # ~10000 ft, 523 mmHg -> 697.27 mb
    
  # depth
  french$depth <- unitted::u(0.16, "m")
  
  # light
  french$app.solar.time <- convert_UTC_to_solartime(french$utc.time, longitude=-106.3, time.type='apparent solar')
  french$light <- convert_SW_to_PAR(calc_solar_insolation(app.solar.time=unitted::v(french$app.solar.time), latitude=41.33, max.insolation=convert_PAR_to_SW(2326), attach.units=TRUE))
  
  # return w/ proper units & columns
  french <- if(attach.units) french else v(french)
  french[c("solar.time","DO.obs","DO.sat","depth","temp.water","light")]
}