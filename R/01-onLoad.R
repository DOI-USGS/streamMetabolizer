#' @import methods
#' @importFrom utils available.packages contrib.url
#' @keywords internal
.onAttach <- function(libname, pkgname) {
  packageStartupMessage("This package is in development. We are using it for our own early applications and welcome flexible, resilient new users who can help us make the package better. Details of the user interface and model implementations will change. Please give us feedback at https://github.com/USGS-R/streamMetabolizer/issues/new.\n")
  packageStartupMessage("This information is preliminary or provisional and is subject to revision. It is being provided to meet the need for timely best science. The information has not received final approval by the U.S. Geological Survey (USGS) and is provided on the condition that neither the USGS nor the U.S. Government shall be held liable for any damages resulting from the authorized or unauthorized use of the information. Although this software program has been used by the USGS, no warranty, expressed or implied, is made by the USGS or the U.S. Government as to the accuracy and functioning of the program and related program material nor shall the fact of distribution constitute any such warranty, and no responsibility is assumed by the USGS in connection therewith.")
  
  # Load deSolve because otherwise after a few model runs we're likely to get
  # the following error
  requireNamespace('deSolve', quietly=TRUE)
  ## Error in .Call("call_rkFixed", as.double(y), as.double(times), Func, Initfunc,  :
  ##   "call_rkFixed" not resolved from current namespace (deSolve)
  ## Error in .C("unlock_solver") :
  ##   "unlock_solver" not resolved from current namespace (deSolve)
  
  # Check whether this package is up to date on GRAN
  GRAN_update_code <- paste0(
    '  update.packages(oldPkgs=c("streamMetabolizer","unitted"),\n',
    '    dependencies=TRUE, repos=c("https://owi.usgs.gov/R", "https://cran.rstudio.com"))')
  tryCatch({
    GRAN_pkg <- available.packages(contrib.url("https://owi.usgs.gov/R"))
    GRAN_version <- package_version(GRAN_pkg[[pkgname, 'Version']])
    local_version <- packageVersion(pkgname)
    if(local_version < GRAN_version) {
      packageStartupMessage(
        'Time to update to ', pkgname, ' version ', GRAN_version, '! You have ', local_version, '. Get stable updates with\n',
        GRAN_update_code)
    }
  }, error=function(e) {
    packageStartupMessage("Can't check GRAN for new package versions just now. We'll try again next time.")
  })
  
  # Check whether this package is up to date on GitHub
  github_owner <- 'USGS-R'
  github_branch <- 'develop'
  github_pkg_ref <- paste0(github_owner,'/',pkgname,'@',github_branch)
  github_update_code <- paste0(
    '  devtools::install_github("',github_pkg_ref,'")')
  if(requireNamespace('devtools', quietly=TRUE)) {
    tryCatch({
      github_ref <- devtools:::github_resolve_ref(
        devtools::github_release(), 
        devtools:::parse_git_repo(github_pkg_ref))$ref
      github_version <- package_version(gsub('v', '', github_ref))
      if(local_version < github_version) {
        packageStartupMessage(
          'New development version of ', pkgname, ' (', github_version, ') is ready! You have ', local_version, '. Get dev updates with\n',
          github_update_code)
      }
    }, error=function(e) {
      packageStartupMessage("Can't check GitHub for new package versions just now. We'll try again next time.")
    })
  }
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
