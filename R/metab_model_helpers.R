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
#'   \item{ \code{date.time} date-time values in solar time, in POSIXct format
#'   with a nominal time zone of UTC.
#'   
#'   \item{ \code{DO.obs} dissolved oxygen concentration observations, \eqn{mg 
#'   O[2] L^{-1}}{mg O2 / L}}
#'   
#'   \item{ \code{DO.sat} dissolved oxygen concentrations if the water were at 
#'   equilibrium saturation \eqn{mg O[2] L^{-1}}{mg O2 / L}}. Calculate using 
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
#' @import dplyr
#' @examples
#' mm_data()
#' mm_data(depth, light, date.time)
mm_data <- function(...) {
  dat <- u(data.frame(
    date.time= u(as.POSIXct("2050-03-14 15:9:27",tz="UTC"), NA), 
    DO.obs=    u(10.1,"mgO2 L^-1"), 
    DO.sat=    u(14.2,"mgO2 L^-1"), 
    depth=     u(0.5,"m"), 
    temp.water=u(21.8,"degC"), 
    light=     u(300.9,"umol m^-2 s^-1")))
  .dots = lazyeval::lazy_dots(...)
  if(length(.dots) == 0) dat else select_(dat, .dots=.dots)
}

#' Evaluate whether the data argument is properly formatted.
#' 
#' Will most often be called from within a metab_model constructor.
#' 
#' @param data a data.frame that has been (or will be) passed into a metab_model
#'   constructor
#' @param metab_class character the class name of the metab_model constructor
#' @examples
#' \dontrun{
#' mm_validate_data(dplyr::select(mm_data(),-temp.water), "metab_mle")
#' }
#' @export
mm_validate_data <- function(data, metab_class) {
  
  # the expectation is set by the default data argument to the specific metabolism class
  expected.data <- formals(metab_class)$data %>% eval()
  
  # check for missing or extra columns
  missing.columns <- setdiff(names(expected.data), names(data))
  extra.columns <- setdiff(names(data), names(expected.data))
  if(length(missing.columns) > 0) {
    stop(paste0("data is missing these columns: ", paste0(missing.columns, collapse=", ")))
  }
  if(length(extra.columns) > 0) {
    stop(paste0("data should omit these extra columns: ", paste0(extra.columns, collapse=", ")))
  }
  
  # put the data columns in the same order as expected.data
  data <- data[names(expected.data)]
  
  # check for units mismatches. column names will already match exactly.
  mismatched.units <- which(get_units(expected.data) != get_units(data))
  if(length(mismatched.units) > 0) {
    data.units <- get_units(data)[mismatched.units]
    expected.units <- get_units(expected.data)[mismatched.units]
    stop(paste0("unexpected units: ", paste0(
      "(", 1:length(mismatched.units), ") ", 
      names(data.units), " = ", data.units, ", expected ", expected.units,
      collapse="; ")))
  }
  
  # return the data, which may have had its columns reordered during validation
  return(data)
}


# mm_validate_day <- function(data, cols) {
#   stop("only partly implemented")
#   # Provide ability to skip a poorly-formatted day for calculating 
#   # metabolism, without breaking the whole loop. Just collect 
#   # problems/errors as a list of strings and proceed. Also collect warnings.
#   stop_strs <- warn_strs <- list()
#   
#   ## Error checks:
#   # Require that the data consist of three consecutive days (10:30 pm on Day 1 to 6 am on Day 3)
#   if(!isTRUE(all.equal(diff(range(day$date.time)) %>% as.numeric(units="days"), 
#                        as.difftime(31.5, units="hours") %>% as.numeric(units="days"), 
#                        tol=as.difftime(31, units="mins") %>% as.numeric(units="days")))) {
#     stop_strs <- c(stop_strs, "incomplete time series")
#   }
#   # Require that on each day date.time has a ~single, ~consistent time step
#   timestep.days <- suppressWarnings(mean(as.numeric(diff(day$date.time), units="days"), na.rm=TRUE))
#   timestep.deviations <- suppressWarnings(diff(range(as.numeric(diff(day$date.time), units="days"), na.rm=TRUE)))
#   if(length(stop_strs) == 0 & is.finite(timestep.days) & is.finite(timestep.deviations)) {
#     # max-min timestep length can't be more than 1% of mean timestep length
#     if((timestep.deviations / timestep.days) > 0.001) { 
#       stop_strs <- c(stop_strs, "uneven timesteps")
#     }
#     # all timesteps per day must add up to 31.5 hrs (31.5/24 days), plus or minus 0.51 hrs
#     if(abs(timestep.days * length(day$date.time) - 31.5/24) > 0.51/24) { 
#       stop_strs <- c(stop_strs, paste0("sum(timesteps) != 31.5 hours"))
#     }
#   } else {
#     stop_strs <- c(stop_strs, "can't measure timesteps")
#   }
#   # Require complete data
#   if(any(is.na(day$DO.obs))) stop_strs <- c(stop_strs, "NAs in DO.obs")
#   if(any(is.na(day$DO.sat))) stop_strs <- c(stop_strs, "NAs in DO.sat")
#   if(any(is.na(day$depth))) stop_strs <- c(stop_strs, "NAs in depth")
#   if(any(is.na(day$temp.water))) stop_strs <- c(stop_strs, "NAs in temp.water")
#   if(any(is.na(day$light))) stop_strs <- c(stop_strs, "NAs in light")
#   
#   # Return
#   list(stop_strs=stop_strs, warn_strs=warn_strs)
# }
# 
