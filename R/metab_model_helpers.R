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
#' @importFrom lazyeval lazy_dots
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
  .dots = lazy_dots(...)
  if(length(.dots) == 0) dat else select_(dat, .dots=.dots)
}
# Because metab_models will call mm_data(...) to define their default data, it
# makes sense to declare all the potential columns as global variables here;
# otherwise we'd need to do it before defining any of those functions.
globalVariables(names(mm_data()))

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


#' Validate one day of data, returning a vector of error strings if needed
#' 
#' Provides ability to skip a poorly-formatted day for calculating metabolism, 
#' without breaking the whole loop. Rather than producing errors, quietly 
#' collects problems/errors as a list of strings for the calling function to 
#' handle.
#' 
#' Assumes that the data have already been validated as in 
#' \code{\link{mm_validate_data}}
#' 
#' @param day data for one day
#' @param tests list of tests to conduct
#' @param day_start argument for 'full_day' test: expected start of day data in 
#'   number of hours from the midnight that begins the modal date. For example, 
#'   day_start=-1.5 indicates that data describing 2006-06-26 should begin at 
#'   2006-06-25 22:30, or as near as possible given the timestep.
#' @param day_end argument for 'full_day' test: expected end of day data in 
#'   number of hours from the midnight that begins the modal date. For example, 
#'   day_end=30 indicates that data describing 2006-06-26 should end at 
#'   2006-06-27 06:00, or as near as possible given the timestep.
#' @param need_complete character vector of the names of columns in day that
#'   must be complete (without NAs)
#' @return character vector of errors, or empty list
#' @export
mm_is_valid_day <- function(day, tests=c('full_day', 'even_timesteps', 'complete_data'), 
                            day_start=-1.5, day_end=30,
                            need_complete=names(day)) {
  
  # check input
  tests <- match.arg(tests, several.ok = TRUE)
  day_start <- as.difftime(day_start, units="hours")
  day_end <- as.difftime(day_end, units="hours")
  
  # initialize vectors
  stop_strs <- character(0)

  # estimate time steps - useful for a few tests
  timestep <- as.difftime(suppressWarnings(mean(as.numeric(diff(v(day$date.time)), units="days"), na.rm=TRUE)), units="days")
  if(!is.finite(timestep)) {
    stop_strs <- c(stop_strs, "can't measure timesteps")
  }

  # Require that the data span the full expected period (e.g., from 10:30pm on
  # preceding day to 6am on following day)
  if('full_day' %in% tests & is.finite(timestep)) {
    date_counts <- table(format(day$date.time, "%Y-%m-%d"))
    date_start <- as.POSIXct(paste0(names(date_counts)[which.max(date_counts)], " 00:00:00"), tz="UTC")
    similar_time <- function(a, b, tol=timestep*0.5) {
      isTRUE(all.equal(as.numeric(a, units="days"), as.numeric(b, units="days"), tol=as.numeric(tol, units="days")))
    }
    if(!similar_time(min(day$date.time)-date_start, day_start))
      stop_strs <- c(stop_strs, "data don't start when expected")
    if(!similar_time(max(day$date.time)-date_start, day_end))
      stop_strs <- c(stop_strs, "data don't end when expected")
  }
  
  # Require that on each day date.time has a ~single, ~consistent time step
  if('even_timesteps' %in% tests & is.finite(timestep)) {
    timestep.deviations <- suppressWarnings(diff(range(as.numeric(diff(v(day$date.time)), units="days"), na.rm=TRUE)))
    if(!is.finite(timestep.deviations)) {
      stop_strs <- c(stop_strs, "can't measure range of timestep lengths")
    } else {
      # max-min timestep length shouldn't be more than 2% of mean timestep length
      if((timestep.deviations / as.numeric(timestep, units="days")) > 0.002) { 
        stop_strs <- c(stop_strs, "uneven timesteps")
      }
    }
  }
  
  if('complete_data' %in% tests) {
    # Require complete data in all columns
    for(col in need_complete) {
      if(any(is.na(day[col]))) 
        stop_strs <- c(stop_strs, paste0("NAs in ", col))
    }
  }
  
  # Return
  stop_strs
}

