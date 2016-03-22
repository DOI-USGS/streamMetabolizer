#' Load a short dataset from Spring Creek
#' 
#' @import dplyr
#' @importFrom unitted u rename_.unitted_data.frame
#' @importFrom utils read.csv
#' @importFrom lubridate with_tz
#' @param attach.units logical. Should units be attached to the data.frame?
#' @return a data.frame, unitted if attach.units==TRUE
load_spring_creek <- function(attach.units=TRUE) {
  # load the file
  file.name <- system.file("extdata", "spring14.csv", package="streamMetabolizer") # data from Spring Creek, Laramie, WY
  utc.time <- oxy <- temp <- solar.time <- app.solar.time <- ".dplyr.var"
  spring <- read.csv(file.name, stringsAsFactors=FALSE, header=TRUE) %>% 
    transmute(
      utc.time = as.POSIXct(time, origin="1970-01-01", tz="UTC"),
      local.time = with_tz(utc.time, "America/Denver"),
      DO.obs = u(oxy, 'mgO2 L^-1'),
      temp.water = u(temp, 'degC')) %>% 
    u() %>%
    mutate(
      DO.sat = calc_DO_at_sat(temp.water=temp.water, pressure.air=u(595, "mmHg")*unitted::u(1.33322368, "mb mmHg^-1")),
      depth = u(rep(0.18, length(temp.water)), "m"),
      solar.time = convert_UTC_to_solartime(utc.time, longitude=-105.6, time.type='mean solar'),
      app.solar.time = convert_UTC_to_solartime(utc.time, longitude=-105.6, time.type='apparent solar'),
      light = convert_SW_to_PAR(calc_solar_insolation(app.solar.time=v(app.solar.time), latitude=41.33, max.insolation=convert_PAR_to_SW(2326), attach.units=TRUE))
    )
    
  # return w/ proper units & columns
  spring <- if(attach.units) spring else v(spring)
  spring[c("solar.time","DO.obs","DO.sat","depth","temp.water","light")]
}