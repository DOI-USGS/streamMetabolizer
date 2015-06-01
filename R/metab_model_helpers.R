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
#' @examples
#' mm_data()
#' dplyr::select(mm_data(), depth, light)
mm_data <- function() {
  u(data.frame(
    date.time= u(as.POSIXct("2050-03-14 15:9:27",tz="UTC"), NA), 
    DO.obs=    u(10.1,"mgO2 L^-1"), 
    DO.sat=    u(14.2,"mgO2 L^-1"), 
    depth=     u(0.5,"m"), 
    temp.water=u(21.8,"degC"), 
    light=     u(300.9,"umol m^-2 s^-1")))
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