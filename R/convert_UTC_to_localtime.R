#' Convert time from UTC to local time.
#' 
#' Convert time from UTC to local time, either standard or with daylight 
#' savings. Recommended for post-analysis visualization only; most functions in 
#' streamMetabolizer use times in UTC.
#' 
#' @param date.time POSIXct object the date and time in UTC
#' @param latitude numeric, in degrees, either positive and unitted ("degN" or 
#'   "degS") or with sign indicating direction (positive = North)
#' @param longitude numeric, in degrees, either positive and unitted ("degE" or 
#'   "degW") or with sign indicating direction (positive = East)
#' @param time.type character. The type of time zone desired - either standard 
#'   time without any daylight savings time or daylight time where daylight
#'   savings is on during the appropriate days
#' @importFrom lubridate with_tz
#' @importFrom unitted u v
#' @references 
#' http://stackoverflow.com/questions/23414340/convert-to-local-time-zone-using-latitude-and-longitude
#' @export
convert_UTC_to_localtime <- function(date.time, latitude, longitude, time.type=c("standard local", "daylight local")) {
  
  # format checking - require expected time.type, tz=UTC, and expected units
  time.type <- match.arg(time.type)
  if(is.unitted(date.time)) date.time <- v(date.time)
  if(class(date.time)[1] != "POSIXct") stop("expecting date.time as a POSIXct object")
  if(!(tz(date.time) %in% c("GMT","Etc/GMT-0","Etc/GMT+0","UTC"))) stop("expecting tz=UTC")
  # alternative to above: date.time <- with_tz(date.time, tzone="UTC") # hidden feature, or bad/weak error checking?
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
  
  # ask the cache or Google for a time offset. The Google API limits requests to
  # 5/sec, so caching helps limit the number of requests we have to make.
  lookup_key <- sprintf("%.10f,%.10f", latitude, longitude)
  tz_info <- pkg.env$tz_lookups[[lookup_key]]
  if(is.null(tz_info) || tz_info$retry > 0) { # give lookup a second chance if it didn't come out right
    retry <- tz_info$retry # save the retry info
    tz_info <- lookup_google_timezone(latitude, longitude) # ask Google
    tz_info$retry <- if(tz_info$tz == "" || is.na(tz_info$std_offset)) retry - 1 else 0 # check the output for validity
    pkg.env$tz_lookups[[lookup_key]] <- tz_info # update tz_lookups
    if(tz_info$retry > 0) stop("sorry, could not find time zone for specified lat/long") # if there was a problem, throw the error now
  }
  
  # return in specified format
  if(time.type == "daylight local") {
    lubridate::with_tz(date.time, tz_info$tz)
  } else {
    # "POSIX has positive signs west of Greenwich" - http://opensource.apple.com/source/system_cmds/system_cmds-230/zic.tproj/datfiles/etcetera
    std.tz <- sprintf("Etc/GMT%s%d", if(tz_info$std_offset > u(0, "hours")) "-" else "+", abs(as.numeric(v(tz_info$std_offset))))
    if(std.tz %in% c("Etc/GMT+0", "Etc/GMT-0")) std.tz <- "UTC"
    lubridate::with_tz(date.time, std.tz)
  }  
}

#' Convert time from local time to UTC.
#' 
#' Convert time from local time (either standard or with daylight savings) to 
#' UTC.
#' 
#' @param local.time POSIXct date+time of interest, already in local time as
#'   specified by the tz attribute
#' @importFrom lubridate with_tz
#' @references 
#' http://stackoverflow.com/questions/23414340/convert-to-local-time-zone-using-latitude-and-longitude
#' @export
convert_localtime_to_UTC <- function(local.time) {
  return(with_tz(local.time, "UTC"))
}

