#' @title Convert DateTime from GMT to local solar time
#'   
#' @param date.time date-time values in POSIXct format.
#' @param time.type character. "apparent solar", i.e. true solar time, is noon 
#'   when the sun is at its zenith. "mean solar" approximates apparent solar 
#'   time but with noons exactly 24 hours apart.
#' @return a POSIXct object that says it's in tz="GMT" but that's actually in
#'   solar time, with noon being very close to solar noon
#' @importFrom lubridate tz with_tz
#' @importFrom unitted u v is.unitted
#' @export
#' @references Yard, Bennett, Mietz, Coggins, Stevens, Hueftle, and Blinn. 2005.
#'   Influence of topographic complexity on solar insolation estimates for the 
#'   Colorado River, Grand Canyon, AZ. Ecological Modelling.
convert_GMT_to_solartime <- function(date.time, longitude, time.type=c("apparent solar", "mean solar")){
  time.type <- match.arg(time.type)
  
  # format checking - require tz=GMT and expected units
  if(is.unitted(date.time)) date.time <- v(date.time)
  if(class(date.time)[1] != "POSIXct") stop("expecting date.time as a POSIXct object")
  if(!(tz(date.time) %in% c("GMT","Etc/GMT-0","Etc/GMT+0"))) stop("expecting tz=GMT")
  # alternative to above: date.time <- with_tz(date.time, tzone="GMT") # hidden feature, or bad/weak error checking?
  if(is.unitted(longitude)) {
    if(get_units(longitude) == "degW") longitude <- u(-1*longitude, "degE")
    verify_units(longitude, "degE")
  } else {
    longitude <- u(longitude, "degE")
  }
  
  # calculate mean.solar time, which approximates solar noon at clock noon to
  # within ~20 minutes
  longitude.GMT <- u(0, "degE")
  time.adjustment <- u(3.989, "mins degE^-1")*(longitude.GMT + longitude)
  mean.solar <- date.time + as.difftime(time.adjustment, units=get_units(time.adjustment))
  
  # either return mean solar time or adjust to true (apparent) solar time
  if(time.type=="mean solar") {
    return(mean.solar)
  } else { # always "apparent solar"
    # Use the equation of time to compute the discrepancy between apparent and
    # mean solar time. E is in minutes.
    jday <- convert_date_to_doyhr(mean.solar) - 1 # subtract 1 for jdays between 0 and 364
    E <- u(9.87*sin(to_radians((2*360*(jday-81))/365)) - 
             7.53*cos(to_radians((360*(jday-81))/365)) - 
             1.5*sin(to_radians((360*(jday-81))/365)), 
           units="mins") # Equation of time as in Yard et al. 2005
    return(mean.solar + as.difftime(E, units=get_units(E)))
  }
}

#' @title Convert DateTime from local solar time to GMT
#'   
#' @param date.time date-time values in POSIXct format. Timezone will be
#'   assumed/coerced to GMT.
#' @param time.type character. "apparent solar", i.e. true solar time, is noon 
#'   when the sun is at its zenith. "mean solar" approximates apparent solar 
#'   time but with noons exactly 24 hours apart.
#' @return a POSIXct object in GMT
#' @importFrom lubridate force_tz
#' @importFrom unitted u v is.unitted
#' @export
#' @references Yard, Bennett, Mietz, Coggins, Stevens, Hueftle, and Blinn. 2005.
#'   Influence of topographic complexity on solar insolation estimates for the 
#'   Colorado River, Grand Canyon, AZ. Ecological Modelling.
convert_solartime_to_GMT <- function(date.time, longitude, time.type=c("apparent solar", "mean solar")) {
  pretend.gmt <- lubridate::force_tz(date.time, "GMT")
  conversion <- pretend.gmt - convert_GMT_to_solartime(pretend.gmt, longitude, time.type)
  return(date.time + conversion)
}
