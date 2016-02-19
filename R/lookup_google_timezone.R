#' Use Google API to determine local time zone
#' 
#' This function uses two packages, \code{RCurl} and \code{XML}, that are not
#' required for the \code{streamMetabolizer} package as a whole. If these are
#' not already installed, run \code{install.packages(c('RCurl','XML'))} before
#' calling \code{lookup_google_timezone}.
#' 
#' Some parameter definitions below are copied directly from the API webpage.
#' 
#' @param latitude degrees latitude (positive for north) of the location to look
#'   up.
#' @param longitude degrees longitude (positive for east) of the location to 
#'   look up.
#' @param timestamp POSIXct representation of a time - determines daylight 
#'   savings offset, if any. the default is Jan 1 for northern latitudes and 
#'   July 1 for southern latitudes, i.e., a time surely not during daylight 
#'   savings time.
#' @references https://developers.google.com/maps/documentation/timezone/
#' @export
lookup_google_timezone <- function(
  latitude, longitude, 
  timestamp=if(latitude >= 0) as.POSIXct("2015-01-01 00:00:00", tz="UTC") else as.POSIXct("2015-07-01 00:00:00", tz="UTC")) {
  
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