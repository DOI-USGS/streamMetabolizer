#' @title Convert DateTime from GMT to local solar time
#'   
#' @param date.time date-time values in POSIXct format.
#' @param time.type character. "apparent solar", i.e. true solar time, is noon 
#'   when the sun is at its zenith. "mean solar" approximates apparent solar 
#'   time but with noons exactly 24 hours apart.
#' @return a POSIXct object that says it's in tz="GMT" but that's actually in
#'   solar time, with noon being very close to solar noon
#' @importFrom lubridate tz with_tz
#' @export
#' @references Yard, Bennett, Mietz, Coggins, Stevens, Hueftle, and Blinn. 2005.
#'   Influence of topographic complexity on solar insolation estimates for the 
#'   Colorado River, Grand Canyon, AZ. Ecological Modelling.
convert_GMT_to_solartime <- function(date.time, longitude, time.type=c("apparent solar", "mean solar"), ...){
  time.type <- match.arg(time.type)
  
  # format checking - require tz=GMT and expected units
  if(class(v(date.time))[1] != "POSIXct") stop("expecting date.time as a POSIXct object")
  if(tz(date.time) != "GMT") date.time <- with_tz(date.time, tzone="GMT")
  if(!is.unitted(date.time)) date.time <- u(date.time)
  if(is.unitted(longitude)) {
    if(get_units(longitude) == "degE") longitude <- u(-longitude, "degW")
    verify_units(longitude, "degW")
  } else {
    longitude <- u(longitude, "degW")
  }
  
  # calculate mean.solar time, which approximates solar noon at clock noon to
  # within ~20 minutes
  longitude.GMT <- u(0, "degW")
  time.adjustment <- u(3.989, "mins degW^-1")*(longitude.GMT-longitude)
  mean.solar <- date.time + as.difftime(time.adjustment, units=get_units(time.adjustment))
  
  # either return mean solar time or adjust to true (apparent) solar time
  if(time.type=="mean solar") {
    return(mean.solar)
  } else { # always "apparent solar"
    # Use the equation of time to compute the discrepancy between apparent and
    # mean solar time. E is in minutes.
    jday <- convert_date_to_doyhr(mean.solar) - 1 # subtract 1 for jdays between 0 and 364
    E <- u(9.87*sin(radi((2*360*(jday-81))/365)) - 
             7.53*cos(radi((360*(jday-81))/365)) - 
             1.5*sin(radi((360*(jday-81))/365)), 
           units="mins") # Equation of time as in Yard et al. 2005
    return(mean.solar + as.difftime(E, units=get_units(E)))
  }
}

#' Convert time from GMT to local time.
#' 
#' Convert time from GMT to local time, either standard or with daylight
#' savings. Recommended for post-analysis visualization only; most functions in
#' streamMetabolizer use times in GMT.
convert_GMT_to_localtime <- function(date, time.type=c("standard local", "daylight local")) {
  stop("not yet implemented")
}