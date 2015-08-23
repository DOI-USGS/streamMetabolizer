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
#' @param day_start expected start of day data in number of hours from the
#'   midnight that begins the modal date. For example, day_start=-1.5 indicates
#'   that data describing 2006-06-26 should begin at 2006-06-25 22:30, or at the
#'   first observation time that occurs after that time if day_start doesn't
#'   fall exactly on an observation time.
#' @param day_end expected end of day data in number of hours from the midnight
#'   that begins the modal date. For example, day_end=30 indicates that data
#'   describing 2006-06-26 should end at 2006-06-27 06:00, or at the last
#'   observation time that occurs before that time if day_end doesn't fall
#'   exactly on an observation time.
#' @param timestep_days the expected timestep length in fraction of a day; for 
#'   example, a 1-hour timestep is 1/24 is 0.0416667. This is calculated within 
#'   the function if timestep_days is NA.
#' @param need_complete character vector of the names of columns in day that 
#'   must be complete (without NAs)
#' @return character vector of errors, or empty list
#' @export
mm_is_valid_day <- function(day, tests=c('full_day', 'even_timesteps', 'complete_data'), 
                            day_start=-1.5, day_end=30, timestep_days=NA,
                            need_complete=names(day)) {
  
  # check input
  tests <- match.arg(tests, several.ok = TRUE)
  day_start <- as.difftime(day_start, units="hours")
  day_end <- as.difftime(day_end, units="hours")
  
  # initialize vectors
  stop_strs <- character(0)
  
  # estimate time steps - useful for a few tests
  timestep.days <- if(is.na(timestep_days)) suppressWarnings(mean(as.numeric(diff(v(day$local.time)), units="days"), na.rm=TRUE)) else timestep_days
  timestep <- as.difftime(timestep.days, units="days")
  if(!is.finite(timestep)) {
    stop_strs <- c(stop_strs, "can't measure timesteps")
  }

  # Require that the data span the full expected period (e.g., from 10:30pm on
  # preceding day to 6am on following day)
  if('full_day' %in% tests & is.finite(timestep)) {
    date_counts <- table(format(day$local.time, "%Y-%m-%d"))
    local_date <- names(date_counts)[which.max(date_counts)]
    date_start <- as.POSIXct(paste0(local_date, " 00:00:00"), tz="UTC")
    similar_time <- function(a, b, tol=timestep) {
      abs(as.numeric(a, units="days") - as.numeric(b, units="days")) < as.numeric(tol, units="days")
    }
    if(!similar_time(min(day$local.time)-date_start, day_start))
      stop_strs <- c(stop_strs, "data don't start when expected")
    if(!similar_time(max(day$local.time)-date_start, day_end))
      stop_strs <- c(stop_strs, "data don't end when expected")
  }
  
  # Require that on each day local.time has a ~single, ~consistent time step
  if('even_timesteps' %in% tests & is.finite(timestep)) {
    timestep.deviations <- suppressWarnings(diff(range(as.numeric(diff(v(day$local.time)), units="days"), na.rm=TRUE)))
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
  
  # Return the stop strings with Date attribute if there was a problem; otherwise, return TRUE
  if(length(stop_strs) > 0) {
    stop_strs
  } else {
    TRUE
  }
}
