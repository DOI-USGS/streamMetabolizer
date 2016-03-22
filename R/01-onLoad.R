#'@import methods
#'@keywords internal
.onAttach <- function(libname, pkgname) {
  packageStartupMessage("This information is preliminary or provisional and is subject to revision. It is being provided to meet the need for timely best science. The information has not received final approval by the U.S. Geological Survey (USGS) and is provided on the condition that neither the USGS nor the U.S. Government shall be held liable for any damages resulting from the authorized or unauthorized use of the information. Although this software program has been used by the USGS, no warranty, expressed or implied, is made by the USGS or the U.S. Government as to the accuracy and functioning of the program and related program material nor shall the fact of distribution constitute any such warranty, and no responsibility is assumed by the USGS in connection therewith.")
}
library(methods)

#' Define a package environment for storing data specific to a project during an
#' R session
#' 
#' @return the package environment
#' @keywords internal
define_pkg_env <- function() {
  pkg.env <- new.env()
  pkg.env$tz_lookups <- list(
    # populate with values that are used in test-convert.R and load_french_creek.R
    "51.4800000000,-0.0000000000"=list(tz="Europe/London", dst_offset=u(0,"hours"), std_offset=u(0,"hours"), retry=0),
    "41.0000000000,105.3000000000"=list(tz="Asia/Shanghai", dst_offset=u(0,"hours"), std_offset=u(8,"hours"), retry=0),
    "37.0000000000,-105.3000000000"=list(tz="America/Denver", dst_offset=u(0,"hours"), std_offset=u(-7,"hours"), retry=0),
    "34.0000000000,-80.0000000000"=list(tz="America/New_York", dst_offset=u(0,"hours"), std_offset=u(-5,"hours"), retry=0),
    "44.3625940000,-106.7530990000"=list(tz="America/Denver", dst_offset=u(0,"hours"), std_offset=u(-7,"hours"), retry=0),
    "41.3300000000,-106.3000000000"=list(tz="America/Denver", dst_offset=u(0,"hours"), std_offset=u(-7,"hours"), retry=0) # French Creek
  )
  return(pkg.env)
}
pkg.env <- define_pkg_env()
