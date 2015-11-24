#' Return the data types that may be used by metab_models using the 
#' metab_model_interface.
#' 
#' @description Produces a unitted data.frame with the column names, units, and
#' data format to be used by metab_models that comply strictly with the 
#' metab_model_interface. These are the columns that may be included:
#' 
#' \itemize{
#' 
#' \item{ \code{local.time} date-time values in solar mean time (local, 
#' longitude-specific, NON-DAYLIGHT-SAVINGS time), in POSIXct format with a 
#' nominal time zone of UTC.}
#' 
#' \item{ \code{DO.obs} dissolved oxygen concentration observations, \eqn{mg 
#' O[2] L^{-1}}{mg O2 / L}}
#' 
#' \item{ \code{DO.sat} dissolved oxygen concentrations if the water were at 
#' equilibrium saturation \eqn{mg O[2] L^{-1}}{mg O2 / L}. Calculate using 
#' \link{calc_DO_at_sat}}
#' 
#' \item{ \code{depth} stream depth, \eqn{m}{m}}.
#' 
#' \item{ \code{temp.water} water temperature, \eqn{degC}}.
#' 
#' \item{ \code{light} photosynthetically active radiation, \eqn{\mu mol\ m^{-2}
#' s^{-1}}{micro mols / m^2 / s}}
#' 
#' \item{ \code{local.date} dates of interest in Date format}
#' 
#' \item{ \code{DO.obs} dissolved oxygen concentration observations, \eqn{mg 
#' O[2] L^{-1}}{mg O2 / L}}
#'
#' \item{ \code{GPP} daily estimates of GPP}
#' 
#' \item{ \code{ER} daily estimates of ER}
#' 
#' \item{ \code{K600} daily estimates of K600}
#' 
#' \item{ \code{discharge.daily} daily mean river discharge}
#' 
#' \item{ \code{velocity.daily} daily mean river flow velocity}
#' 
#' }
#' 
#' @details Most models will require a subset of these data columns. Specialized
#' models may deviate from this format, but this is discouraged.
#' 
#' @param ... column names to select, as passed to \code{\link[dplyr]{select}}
#' @param optional one or more character strings listing the columns, if any, 
#'   that may be excluded. If 'all', the entire data.frame may be omitted. If 
#'   'none', the entire data.frame must be included as prototyped. If specific 
#'   column names are given, those columns may be omitted entirely or passed to 
#'   the metab_model() call as all NAs.
#' @return data data.frame with columns as in the description
#'   
#' @export
#' @importFrom unitted u
#' @importFrom lazyeval lazy_dots
#' @import dplyr
#' @examples
#' # all possible columns
#' mm_data()
#' 
#' # columns typical of instantaneous data
#' mm_data(local.time, DO.obs, DO.sat, depth, temp.water, light)
#' 
#' # columns typical of daily data
#' mm_data(local.date, K600, discharge, velocity)
mm_data <- function(..., optional='none') {
  dat <- u(data.frame(
    local.time=u(as.POSIXct("2050-03-14 15:10:00",tz="UTC"), NA), 
    DO.obs=    u(10.1,"mgO2 L^-1"), 
    DO.sat=    u(14.2,"mgO2 L^-1"), 
    depth=     u(0.5,"m"), 
    temp.water=u(21.8,"degC"), 
    light=     u(300.9,"umol m^-2 s^-1"), 
    local.date=u(as.Date("2050-03-14", tz="UTC"), NA), 
    DO.mod.1=  u(7.5,"mgO2 L^-1"),
    GPP=       u(5,"gO2 m^-2 d^-1"), 
    ER=        u(5,"gO2 m^-2 d^-1"), 
    K600=      u(5,"d^-1"), 
    K600.lower=u(4.5,"d^-1"), 
    K600.upper=u(5.6,"d^-1"), 
    discharge= u(9,"m^3 s^-1"), 
    velocity=  u(2,"m s^-1")))
  .dots <- lazy_dots(...)
  .nulldot <- length(.dots) == 1 && is.null(.dots[[1]]$expr)
  dat <- if(isTRUE(.nulldot)) {
    u(NULL)
  } else if(length(.dots) == 0) {
    dat 
  } else {
    select_(dat, .dots=.dots)
  }

  # add information about which columns, if any, are optional.
  optional <- if(missing(optional)) {
    if(isTRUE(.nulldot)) {
      'all'
    } else {
      'none' 
    }
  } else {
    opt <- match.arg(optional, choices=c('all','none',names(dat)), several.ok=TRUE)
    if(any(c('all','none') %in% optional) && length(optional) != 1) 
      stop("if optional is 'all' or 'none', it should be length 1")
    if(all(names(dat) %in% opt)) 
      opt <- 'all'
    opt
  }
  attr(dat, 'optional') <- optional
  
  # return
  dat
}

# Because metab_models will call mm_data(...) to define their default data, it
# makes sense to declare all the potential columns as global variables here;
# otherwise we'd need to do it before defining any of those functions.
globalVariables(names(mm_data()))