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
#' @inheritParams mm_model_by_ply_prototype
#' @param tests list of tests to conduct
#' @param timestep_days the expected timestep length in fraction of a day; for 
#'   example, a 1-hour timestep is 1/24 is 0.0416667. This is calculated within 
#'   the function if timestep_days is NA.
#' @return character vector of errors if day is invalid, or TRUE if it's valid
#' @importFrom lubridate tz
#' @export
mm_is_valid_day <- function(day, day_start=4, day_end=27.99, 
                            tests=c('full_day', 'even_timesteps', 'complete_data'), 
                            timestep_days=NA) {
  
  # check input
  tests <- match.arg(tests, several.ok = TRUE)
  day_start <- as.difftime(day_start, units="hours")
  day_end <- as.difftime(day_end, units="hours")
  
  # initialize vectors
  stop_strs <- character(0)
  
  # estimate time steps - useful for a few tests
  timestep.days <- if(is.na(timestep_days)) suppressWarnings(mean(as.numeric(diff(v(day$solar.time)), units="days"), na.rm=TRUE)) else timestep_days
  timestep <- as.difftime(timestep.days, units="days")
  if(!is.finite(timestep)) {
    stop_strs <- c(stop_strs, "can't measure timesteps")
  }

  # Require that the data span the full expected period (e.g., from 10:30pm on
  # preceding day to 6am on following day)
  if('full_day' %in% tests & is.finite(timestep)) {
    date_counts <- table(format(day$solar.time, "%Y-%m-%d"))
    ply_date <- names(date_counts)[which.max(date_counts)]
    date_start <- as.POSIXct(paste0(ply_date, " 00:00:00"), tz=lubridate::tz(v(day$solar.time)))
    similar_time <- function(a, b, tol=timestep) {
      abs(as.numeric(a, units="days") - as.numeric(b, units="days")) < as.numeric(tol, units="days")
    }
    if(!similar_time(min(day$solar.time)-date_start, day_start))
      stop_strs <- c(stop_strs, "data don't start when expected")
    if(!similar_time(max(day$solar.time)-date_start, day_end))
      stop_strs <- c(stop_strs, "data don't end when expected")
  }
  
  # Require that on each day solar.time has a ~single, ~consistent time step
  if('even_timesteps' %in% tests & is.finite(timestep)) {
    timestep.deviations <- suppressWarnings(diff(range(as.numeric(diff(v(day$solar.time)), units="days"), na.rm=TRUE)))
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
    for(col in names(day)) {
      if(any(is.na(day[col]))) 
        stop_strs <- c(stop_strs, paste0("NAs in ", col))
    }
  }
  
  # Return the stop strings with Date attribute if there was a problem; otherwise, return TRUE
  if(length(stop_strs) > 0) {
    stop_strs
  } else {
    TRUE
  }
}
