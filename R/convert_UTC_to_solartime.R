#' Convert DateTime from UTC to local solar time
#' 
#' Convert DateTime from UTC to local solar time, which may be either apparent 
#' solar (perfect match between noon and solar zenith) or mean solar (exactly 24
#' hours between solar noons).
#' 
#' @param date.time date-time values in POSIXct format and UTC timezone.
#' @param longitude numeric, in degrees, either positive and unitted ("degE" or 
#'   "degW") or with sign indicating direction (positive = East)
#' @param time.type character. "apparent solar", i.e. true solar time, is noon 
#'   when the sun is at its zenith. "mean solar" approximates apparent solar 
#'   time but with noons exactly 24 hours apart. Elsewhere in this package,
#'   variables named "solar.time" are mean solar time, whereas "app.solar.time"
#'   is apparent solar and "any.solar.time" is either.
#' @return a POSIXct object that says it's in tz="UTC" but that's actually in 
#'   solar time, with noon being very close to solar noon
#' @importFrom lubridate tz with_tz
#' @importFrom unitted u v is.unitted
#' @export
#' @references Yard, Bennett, Mietz, Coggins, Stevens, Hueftle, and Blinn. 2005.
#'   Influence of topographic complexity on solar insolation estimates for the 
#'   Colorado River, Grand Canyon, AZ. Ecological Modelling.
convert_UTC_to_solartime <- function(date.time, longitude, time.type=c("apparent solar", "mean solar")){
  time.type <- match.arg(time.type)
  
  # format checking - require tz=UTC and expected units
  if(is.unitted(date.time)) date.time <- v(date.time)
  if(class(date.time)[1] != "POSIXct") stop("expecting date.time as a POSIXct object")
  if(!(tz(date.time) %in% c("GMT","Etc/GMT-0","Etc/GMT+0","UTC"))) stop("expecting tz=UTC")
  # alternative to above: date.time <- with_tz(date.time, tzone="UTC") # hidden feature, or bad/weak error checking?
  if(is.unitted(longitude)) {
    if(get_units(longitude) == "degW") longitude <- u(-1*longitude, "degE")
    verify_units(longitude, "degE")
  } else {
    longitude <- u(longitude, "degE")
  }
  
  # calculate mean.solar time, which approximates solar noon at clock noon to
  # within ~20 minutes
  longitude.UTC <- u(0, "degE")
  time.adjustment <- u(3.989, "mins degE^-1")*(longitude.UTC + longitude)
  mean.solar <- date.time + as.difftime(time.adjustment, units=get_units(time.adjustment))
  
  # either return mean solar time or adjust to true (apparent) solar time
  if(time.type=="mean solar") {
    out <- mean.solar
  } else { # always "apparent solar"
    # Use the equation of time to compute the discrepancy between apparent and
    # mean solar time. E is in minutes.
    jday <- convert_date_to_doyhr(mean.solar) - 1 # subtract 1 for jdays between 0 and 364
    E <- u(9.87*sin(to_radians((2*360*(jday-81))/365)) - 
             7.53*cos(to_radians((360*(jday-81))/365)) - 
             1.5*sin(to_radians((360*(jday-81))/365)), 
           units="mins") # Equation of time as in Yard et al. 2005
    out <- mean.solar + as.difftime(E, units=get_units(E))
  }
  
  if(is.unitted(date.time) && is.unitted(longitude)) out <- u(out)
  return(out)
}

#' Convert DateTime from local solar time to UTC
#' 
#' Convert DateTime to UTC from local solar time, which may be either apparent 
#' solar (perfect match between noon and solar zenith) or mean solar (exactly 24
#' hours between solar noons).
#' 
#' @param any.solar.time either apparent or mean solar time (specified by 
#'   time.type); date-time values in POSIXct format. Timezone must be UTC.
#' @param longitude numeric, in degrees, either positive and unitted ("degE" or 
#'   "degW") or with sign indicating direction (positive = East), describing 
#'   location of the site
#' @param time.type character indicating whether any.solar.time values are in 
#'   apparent or mean solar time. "apparent solar", i.e. true solar time, is 
#'   noon when the sun is at its zenith. "mean solar" approximates apparent 
#'   solar time but with noons exactly 24 hours apart.
#' @return a POSIXct object in UTC
#' @export
#' @references Yard, Bennett, Mietz, Coggins, Stevens, Hueftle, and Blinn. 2005.
#'   Influence of topographic complexity on solar insolation estimates for the 
#'   Colorado River, Grand Canyon, AZ. Ecological Modelling.
convert_solartime_to_UTC <- function(any.solar.time, longitude, time.type=c("apparent solar", "mean solar")) {
  conversion <- any.solar.time - convert_UTC_to_solartime(any.solar.time, longitude, time.type)
  return(any.solar.time + conversion)
}
