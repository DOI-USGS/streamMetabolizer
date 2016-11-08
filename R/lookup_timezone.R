#' Determine the local time zone from the coordinates
#' 
#' Uses the Google API (and/or the package cache) to determine the local
#' timezone name, offset, and DST offset of a site
#' 
#' @param latitude degrees latitude (positive for north) of the location to look
#'   up.
#' @param longitude degrees longitude (positive for east) of the location to 
#'   look up.
#' @export
#' @examples 
#' lookup_timezone(41.33, -106.3)
lookup_timezone <- function(latitude, longitude) {
  # ask the cache for a time offset. The Google API limits requests to 5/sec, so
  # caching helps limit the number of requests we have to make
  lookup_key <- sprintf("%.10f,%.10f", latitude, longitude)
  tz_info <- pkg.env$tz_lookups[[lookup_key]]
  
  # if it's not in the cache, go to Google. give lookup another chance if it
  # didn't come out right and we haven't tried too many times already
  if(is.null(tz_info) || tz_info$retry > 0) {
    retry <- if(is.null(tz_info)) 3 else tz_info$retry # set/extract the retry info
    tz_info <- tryCatch(
      lookup_google_timezone(latitude, longitude),
      error=function(e) list(tz="", dst_offset=NA, std_offset=NA))
    failure <- tz_info$tz == "" || is.na(tz_info$std_offset)
    tz_info$retry <- if(failure) retry - 1 else 0 # check the output for validity
    pkg.env$tz_lookups[[lookup_key]] <- tz_info # update tz_lookups
  }
  # failes here if we've failed this time or given up completely
  if(tz_info$tz == "") stop("sorry, could not find time zone for specified lat/long")
  tz_info
}

#' Use Google API to determine local time zone
#' 
#' This function uses two packages, \code{RCurl} and \code{XML}, that are not
#' required for the \code{streamMetabolizer} package as a whole. If these are
#' not already installed, run \code{install.packages(c('RCurl','XML'))} before
#' calling \code{lookup_google_timezone}.
#' 
#' Some parameter definitions below are copied directly from the API webpage.
#' 
#' @inheritParams lookup_timezone
#' @param timestamp POSIXct representation of a time - determines daylight 
#'   savings offset, if any. the default is Jan 1 for northern latitudes and 
#'   July 1 for southern latitudes, i.e., a time surely not during daylight 
#'   savings time.
#' @references https://developers.google.com/maps/documentation/timezone/
#' @export
lookup_google_timezone <- function(
  latitude, longitude, 
  timestamp=if(latitude >= 0) as.POSIXct("2015-01-01 00:00:00", tz="UTC") else as.POSIXct("2015-07-01 00:00:00", tz="UTC")) {
  
  called_as_internal <- all(c(':::','streamMetabolizer') %in% as.character(sys.call()[[1]])) ||
    any(sapply(sys.calls()[-sys.nframe()], function(sc) tail(as.character(sc[[1]], 1))) %in% ls(envir = asNamespace("streamMetabolizer")))
  if(!called_as_internal) {
    .Deprecated('lookup_timezone') # deprecation plan is to make this internal. everybody can and should use the lookup_timezone instead
  }
  
  
  # check for required packages specific to this function
  if(!requireNamespace("RCurl", quietly = TRUE)) {
    stop("the RCurl package must be installed to use this function")
  }
  if(!requireNamespace("XML", quietly = TRUE)) {
    stop("the XML package must be installed to use this function")
  }
  
  # ask google
  api.url <- sprintf("https://maps.googleapis.com/maps/api/timezone/xml?location=%s,%s&timestamp=%d", 
                     # &sensor=false - "The Google Maps API previously required that you include the sensor parameter to indicate whether your application used a sensor to determine the user's location. This parameter is no longer required."
                     latitude, 
                     longitude, 
                     as.numeric(as.POSIXct(timestamp, origin="1970-01-01 00:00:00 UTC")))
  api.out <- RCurl::getURL(api.url, .opts = list(ssl.verifypeer = FALSE))           
  out.parsed <- XML::xmlParse(api.out)
  return(list(
    tz = out.parsed[["string(//time_zone_id)"]],
    dst_offset = u(as.numeric(out.parsed[["string(//dst_offset)"]])/3600, "hours"),
    std_offset = u(as.numeric(out.parsed[["string(//raw_offset)"]])/3600, "hours")
  ))
}
