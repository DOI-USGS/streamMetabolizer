#' Return the data types that may be used by metab_models using the 
#' metab_model_interface.
#' 
#' Produces a unitted data.frame with the column names, units, and data format 
#' to be used by metab_models that comply strictly with the 
#' metab_model_interface.
#' 
#' Most models will require a subset of these data columns. Specialized models 
#' may deviate from this format, but this is discouraged.
#' 
#' @param ... column names to select, as passed to \code{\link[dplyr]{select}}
#' @return data data.frame with columns \itemize{
#'   
#'   \item{ \code{local.time} date-time values in local, NON-DAYLIGHT-SAVINGS
#'   time, in POSIXct format with the true local tz format.}
#'   
#'   \item{ \code{solar.time} date-time values in solar time, in POSIXct format 
#'   with a nominal time zone of UTC.}
#'   
#'   \item{ \code{DO.obs} dissolved oxygen concentration observations, \eqn{mg 
#'   O[2] L^{-1}}{mg O2 / L}}
#'   
#'   \item{ \code{DO.sat} dissolved oxygen concentrations if the water were at 
#'   equilibrium saturation \eqn{mg O[2] L^{-1}}{mg O2 / L}. Calculate using 
#'   \link{calc_DO_at_sat}}
#'   
#'   \item{ \code{depth} stream depth, \eqn{m}{m}}.
#'   
#'   \item{ \code{temp.water} water temperature, \eqn{degC}}.
#'   
#'   \item{ \code{light} photosynthetically active radiation, \eqn{\mu mol\ 
#'   m^{-2} s^{-1}}{micro mols / m^2 / s}}
#'   
#'   }
#'   
#' @export
#' @importFrom unitted u
#' @importFrom lazyeval lazy_dots
#' @import dplyr
#' @examples
#' mm_data()
#' mm_data(depth, light, local.time)
mm_data <- function(...) {
  dat <- u(data.frame(
    local.time=u(as.POSIXct("2050-03-14 15:10:00",tz="UTC"), NA), 
    solar.time=u(as.POSIXct("2050-03-14 15:9:27",tz="UTC"), NA), 
    DO.obs=    u(10.1,"mgO2 L^-1"), 
    DO.sat=    u(14.2,"mgO2 L^-1"), 
    depth=     u(0.5,"m"), 
    temp.water=u(21.8,"degC"), 
    light=     u(300.9,"umol m^-2 s^-1")))
  .dots = lazy_dots(...)
  if(length(.dots) == 0) dat else select_(dat, .dots=.dots)
}
# Because metab_models will call mm_data(...) to define their default data, it
# makes sense to declare all the potential columns as global variables here;
# otherwise we'd need to do it before defining any of those functions.
globalVariables(names(mm_data()))