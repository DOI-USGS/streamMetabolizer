#' Convert time from GMT to local time.
#' 
#' Convert time from GMT to local time, either standard or with daylight 
#' savings. Recommended for post-analysis visualization only; most functions in 
#' streamMetabolizer use times in GMT.
#' 
#' @param date.time POSIXct object the date and time in GMT
#' @param latitude numeric, in degrees, either positive and unitted ("degN" or 
#'   "degS") or with sign indicating direction (positive = North)
#' @param longitude numeric, in degrees, either positive and unitted ("degE" or 
#'   "degW") or with sign indicating direction (positive = East)
#' @param time.type character. The type of time zone desired - either standard 
#'   time without any daylight savings time or daylight time where daylight
#'   savings is on during the appropriate days
#' @importFrom lubridate with_tz
#' @importFrom XML xmlParse
#' @importFrom RCurl getURL
#' @importFrom unitted u v
#' @references 
#' http://stackoverflow.com/questions/23414340/convert-to-local-time-zone-using-latitude-and-longitude
#' @export
convert_GMT_to_localtime <- function(date.time, latitude, longitude, time.type=c("standard local", "daylight local")) {
  
  # format checking - require expected time.type, tz=GMT, and expected units
  time.type <- match.arg(time.type)
  if(is.unitted(date.time)) date.time <- v(date.time)
  if(class(date.time)[1] != "POSIXct") stop("expecting date.time as a POSIXct object")
  if(!(tz(date.time) %in% c("GMT","Etc/GMT-0","Etc/GMT+0","UTC"))) stop("expecting tz=GMT")
  # alternative to above: date.time <- with_tz(date.time, tzone="GMT") # hidden feature, or bad/weak error checking?
  with_units <- (is.unitted(longitude) | is.unitted(latitude))
  if(with_units) {
    if(is.unitted(longitude)) {
      if(get_units(longitude) == "degW") longitude <- u(-1*longitude, "degE")
      verify_units(longitude, "degE")
      longitude <- v(longitude)
    } else {
      stop("if either longitude or latitude is unitted, both should be unitted. longitude is not.")
    }
    if(is.unitted(latitude)) {
      if(get_units(latitude) == "degS") latitude <- u(-1*latitude, "degN")
      verify_units(latitude, "degN")
      latitude <- v(latitude)
    } else {
      stop("if either longitude or latitude is unitted, both should be unitted. latitude is not.")
    }    
  }
  
  #' Use Google API to determine local time zone
  #' 
  #' Some parameter definitions below are copied directly from the API webpage.
  #' 
  #' @param latitude degrees latitude (positive for north) of the location to 
  #'   look up.
  #' @param longitude degrees longitude (positive for east) of the location to 
  #'   look up.
  #' @param timestamp POSIXct representation of a time - determines daylight
  #'   savings offset, if any
  #' @references https://developers.google.com/maps/documentation/timezone/
  google_timezone <- function(latitude, longitude, timestamp) {
    api.url <- sprintf("https://maps.googleapis.com/maps/api/timezone/xml?location=%s,%s&timestamp=%d&sensor=false", 
                       latitude, 
                       longitude, 
                       as.numeric(as.POSIXct(timestamp, origin="1970-01-01 00:00:00 GMT")))
    api.out <- RCurl::getURL(api.url, .opts = list(ssl.verifypeer = FALSE))           
    out.parsed <- XML::xmlParse(api.out)
    return(list(
      tz = out.parsed[["string(//time_zone_id)"]],
      dst_offset = u(as.numeric(out.parsed[["string(//dst_offset)"]])/3600, "hours"),
      std_offset = u(as.numeric(out.parsed[["string(//raw_offset)"]])/3600, "hours")
    ))
  }
  # find a single wintertime date to pass to Google (i.e., definitely not daylight savings time)
  date.time.std <- if(latitude >= 0) as.POSIXct("2015-01-01 00:00:00", tz="GMT") else as.POSIXct("2015-07-01 00:00:00", tz="GMT")
  tz.info <- google_timezone(latitude, longitude, date.time.std)
  
  # check the output for validity
  if(tz.info$tz == "" || is.na(tz.info$std_offset)) stop("sorry, could not find time zone for specified lat/long")
  
  # return in specified format
  if(time.type == "daylight local") {
    lubridate::with_tz(date.time, tz.info$tz)
  } else {
    # "POSIX has positive signs west of Greenwich" - http://opensource.apple.com/source/system_cmds/system_cmds-230/zic.tproj/datfiles/etcetera
    std.tz <- sprintf("Etc/GMT%s%d", if(tz.info$std_offset > u(0, "hours")) "-" else "+", abs(as.numeric(v(tz.info$std_offset))))
    if(std.tz %in% c("Etc/GMT+0", "Etc/GMT-0")) std.tz <- "GMT"
    lubridate::with_tz(date.time, std.tz)
  }  
}

#' Convert time from local time to GMT.
#' 
#' Convert time from local time (either standard or with daylight savings) to 
#' GMT.
#' 
#' @param local.time POSIXct date+time of interest, already in local time as
#'   specified by the tz attribute
#' @importFrom lubridate with_tz
#' @references 
#' http://stackoverflow.com/questions/23414340/convert-to-local-time-zone-using-latitude-and-longitude
#' @export
convert_localtime_to_GMT <- function(local.time) {
  return(with_tz(local.time, "GMT"))
}