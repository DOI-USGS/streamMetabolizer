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
  french_creek_units <- c(station=NA, siteno=NA, sonde=NA, date=NA, time=NA, temp='degC', oxy='mg L^-1')
  french <- u(french, french_creek_units)
  
  # remove NA oxys (1658) and remaining duplicates (n=1)
  french <- unique(french[!is.na(french$oxy),])
  
  # rename DO.obs, temp.water
  french <- unitted::rename_.unitted_data.frame(french, DO.obs='oxy', temp.water='temp')
  
  # datetime
  french$local.time <- with_tz(as.POSIXct(paste(french$date, french$time), format="%m/%d/%Y %H:%M:%S", tz="America/Denver"), "MST") # it's in MDT
  
  # DO at sat
  french$DO.sat <- calc_DO_at_sat(temp.water=french$temp.water, pressure.air=u(523, "mmHg")*u(1.33322368, "mb mmHg^-1")) # ~10000 ft, 523 mmHg -> 697.27 mb
    
  # depth
  french$depth <- u(0.16, "m")
  
  # light
  french$utc.time <- convert_localtime_to_GMT(french$local.time)
  french$solar.time <- convert_GMT_to_solartime(french$utc.time, longitude=-106.3, time.type='apparent solar')
  french$light <- convert_SW_to_PAR(calc_solar_insolation(solar.time=v(french$solar.time), latitude=41.33, max.insolation=convert_PAR_to_SW(2326), attach.units=TRUE))
  
  # return w/ proper units & columns
  french <- if(attach.units) french else v(french)
  french[c("local.time","DO.obs","DO.sat","depth","temp.water","light")]
}