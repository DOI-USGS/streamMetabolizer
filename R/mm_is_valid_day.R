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
#' @inheritParams mm_model_by_ply_prototype
#' @inheritParams mm_model_by_ply
#' @param day_tests list of tests to conduct to determine whether each date 
#'   worth of data is valid for modeling. \code{full_day}: Do the data span the 
#'   full expected period (e.g., from 10:30pm on preceding day to 6am on 
#'   following day)? \code{even_timesteps}: are all of the timesteps within the 
#'   day the same length, to within a tolerance of 0.2\% of the timestep length? 
#'   \code{complete_data}: are all columns of input data available at every 
#'   timestep? \code{pos_discharge}: is discharge greater than 0 at every 
#'   timestep? A further test is implied if \code{required_timestep} is a non-NA
#'   numeric.
#' @param required_timestep NA or numeric (length 1). If numeric, the timestep 
#'   length in days that a date must have to pass the validity check (to within 
#'   a tolerance of 0.2\% of the value of \code{required_timestep})
#' @param ply_date the Date this data_ply is intended to match. May be NA
#' @param timestep_days the expected timestep length in fraction of a day; for 
#'   example, a 1-hour timestep is 1/24 is 0.0416667. This is calculated within 
#'   the function if timestep_days is NA. May be supplied as an argument to (1) 
#'   pre-calculate the value for efficiency, or (2) require a specific timestep.
#' @return character vector of errors if day is invalid, or TRUE if it's valid
#' @importFrom lubridate tz
#' @examples
#' mm_is_valid_day(data_metab('1'))
#' mm_is_valid_day(data_metab('1', flaws='missing middle'))
#' mm_is_valid_day(data_metab('1', flaws='missorted'))
#' mm_is_valid_day(data_metab('1', flaws='duplicated'))
#' mm_is_valid_day(data_metab('1', flaws=c('duplicated','missing end')))
#' mm_is_valid_day(data_metab('3'))
#' @export
mm_is_valid_day <- function(
  data_ply, # inheritParams mm_model_by_ply_prototype
  day_start=4, day_end=27.99, # inheritParams mm_model_by_ply
  day_tests=c('full_day', 'even_timesteps', 'complete_data', 'pos_discharge'), 
  required_timestep=NA,
  ply_date=as.Date(format(data_ply[max(1,nrow(data_ply)/2),'solar.time'], "%Y-%m-%d")),
  timestep_days=NA
) {
  
  # check input
  if(!missing(day_tests) && length(day_tests) == 0 && is.na(required_timestep[1])) return(TRUE)
  day_tests <- match.arg(day_tests, several.ok = TRUE)
  day_start <- as.difftime(day_start, units="hours")
  day_end <- as.difftime(day_end, units="hours")
  
  # initialize vectors
  stop_strs <- character(0)
  
  # find the mean timestep if it will be needed
  if(any(c('full_day','even_timesteps') %in% day_tests) || !isTRUE(is.na(required_timestep))) {
    timestep.days <- if(!is.null(timestep_days) && is.na(timestep_days)) {
      tryCatch(
        mm_get_timestep(data_ply$solar.time, format='mean'),
        error=function(e) {
          stop_strs <<- c(stop_strs, e$message)
          NA
        })
    } else {
      timestep_days
    }
    if(length(timestep.days) == 0 || !is.finite(timestep.days)) {
      stop_strs <- c(stop_strs, "no timesteps")
      timestep.days <- NA
    }
  }

  # Require that the data span the full expected period (e.g., from 10:30pm on
  # preceding day to 6am on following day)
  if('full_day' %in% day_tests & is.finite(timestep.days)) {
    # date_counts <- table(format(data_ply$solar.time, "%Y-%m-%d"))
    # ply_date <- names(date_counts)[which.max(date_counts)]
    date_start <- as.POSIXct(paste0(as.character(ply_date), " 00:00:00"), tz=lubridate::tz(v(data_ply$solar.time)))
    similar_time <- function(a, b, tol) {
      abs(as.numeric(a, units="days") - as.numeric(b, units="days")) < as.numeric(tol, units="days")
    }
    if(!similar_time(min(data_ply$solar.time)-date_start, day_start, tol=timestep.days))
      stop_strs <- c(stop_strs, "data don't start when expected")
    if(!similar_time(max(data_ply$solar.time)-date_start, day_end, tol=timestep.days))
      stop_strs <- c(stop_strs, "data don't end when expected")
  }
  
  # Require that on each day solar.time has a ~single, ~consistent time step
  # (with tolerance of 0.2% of mean timestep length)
  if('even_timesteps' %in% day_tests & is.finite(timestep.days)) {
    tryCatch(
      ts <- mm_get_timestep(data_ply$solar.time, format='unique', require_unique=TRUE, tol=0.002*timestep.days),
      error=function(e) {
        stop_strs <<- c(stop_strs, "uneven timesteps")
      })
  }
  
  # Require the actual average timestep to match the value specified in
  # required_timestep (within tolerance of 0.2% of required_timestep length)
  if(!isTRUE(is.na(required_timestep)) & is.finite(timestep.days)) {
    if(length(required_timestep) != 1 || !is.numeric(required_timestep)) {
      stop('expecting required_timestep as single numeric value (or NA)')
    }
    if(!isTRUE(all.equal(timestep.days, required_timestep, tolerance=0.002*required_timestep))) {
      stop_strs <- c(stop_strs, sprintf("timestep != %0.08f d", required_timestep))
    }
  }
  
  # Require complete data in all columns
  if('complete_data' %in% day_tests) {
    for(col in names(data_ply)) {
      if(any(is.na(data_ply[[col]]))) 
        stop_strs <- c(stop_strs, paste0("NAs in ", col))
    }
  }
  
  # Require discharge (if present) to be positive at all times
  if('pos_discharge' %in% day_tests) {
    if('discharge' %in% names(data_ply) && any(!is.na(data_ply$discharge))) {
      if(any(data_ply$discharge[which(!is.na(data_ply$discharge))] <= 0)) 
        stop_strs <- c(stop_strs, "discharge <= 0")
    }
  }
  
  # Return the stop strings if there was a problem; otherwise, return TRUE
  if(length(stop_strs) > 0) {
    stop_strs
  } else {
    TRUE
  }
}
