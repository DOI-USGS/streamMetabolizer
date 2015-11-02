#' Use USGS API (USGS Elevation Point Query Service) to determine approximate 
#' local elevation
#' 
#' This is meant to supply an APPROXIMATE elevation, with no guarantees on 
#' precision or on the lifetime of the API service used by the function. This 
#' function uses two packages, \code{RCurl} and \code{XML}, that are not 
#' required for the \code{streamMetabolizer} package as a whole. If these are 
#' not already installed, run \code{install.packages(c('RCurl','XML'))} before 
#' calling \code{lookup_usgs_elevation}.
#' 
#' @param latitude degrees latitude (positive for north) of the location to look
#'   up.
#' @param longitude degrees longitude (positive for east) of the location to 
#'   look up.
#' @param units character, one of Meters or Feet, specifying the units in which
#'   to return the elevation
#' @references http://ned.usgs.gov/epqs/
#' @importFrom unitted u
#' @export
lookup_usgs_elevation <- function(
  latitude, longitude, units=c("Meters","Feet")) {
  
  # confirm that units are among the accepted values for ned.usgs.gov
  units <- match.arg(units)
  
  # check for required packages specific to this function
  if(!requireNamespace("RCurl", quietly = TRUE)) {
    stop("the RCurl package must be installed to use this function")
  }
  if(!requireNamespace("XML", quietly = TRUE)) {
    stop("the XML package must be installed to use this function")
  }
  
  # ask the USGS
  api.url <- sprintf("http://nationalmap.gov/epqs/pqs.php?x=%f&y=%f&units=%s&output=xml",
                     longitude, latitude, units)
  api.out <- RCurl::getURL(api.url, .opts = list(ssl.verifypeer = FALSE))           
  out.parsed <- XML::xmlParse(api.out)
  out.units <- switch(
    out.parsed[["string(//Units)"]],
    "Feet"="ft",
    "Meters"="m")
  return(list(
    data_source = out.parsed[["string(//Data_Source)"]],
    elevation = u(as.numeric(out.parsed[["string(//Elevation)"]]), out.units)
  ))
}
