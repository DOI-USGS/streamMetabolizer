#' Evaluate whether the data argument is properly formatted.
#' 
#' Will most often be called from within a metab_model constructor.
#' 
#' @inheritParams metab
#' @param metab_class character the class name of the metab_model constructor
#' @param data_tests list of tests to conduct to determine whether the input
#'   data.frames are properly formatted to allow modeling to begin
#' @import dplyr
#' @importFrom lubridate is.POSIXct is.Date
#' @importFrom stats setNames
#' @examples
#' \dontrun{
#' mm_validate_data(dplyr::select(mm_data(),-temp.water), metab_class="metab_mle")
#' }
#' @export
mm_validate_data <- function(
  data=mm_data(NULL), data_daily=mm_data(NULL), #inheritParams metab
  metab_class, data_tests=c('missing_cols','extra_cols','na_times','units')
) {
  
  data_types <- setNames(c("data","data_daily"),c("data","data_daily"))
  dat_all <- lapply(data_types, function(data_type) {
    
    # pick out the data.frame for this loop
    dat <- get(data_type)
    
    # the data expectation is set by the default data argument to the specific metabolism class
    expected.data <- formals(metab_class)[[data_type]] %>% eval()
    optional.data <- attr(expected.data, 'optional')
    if('all' %in% optional.data) optional.data <- c('all', names(expected.data))
    
    # quick return if dat is NULL
    if(is.null(v(dat))) {
      if('all' %in% optional.data) {
        return(dat)
      } else {
        stop(paste0(data_type, " is NULL but required"), call.=FALSE)
      }
    }
    
    # check for missing or extra columns
    if('missing_cols' %in% data_tests) {
      missing.columns <- setdiff(names(expected.data), names(dat))
      missing.columns <- setdiff(missing.columns, optional.data) # optional cols don't count
      if(length(missing.columns) > 0) {
        stop(paste0(data_type, " is missing these columns: ", paste0(missing.columns, collapse=", ")), call.=FALSE)
      }
    }
    if('extra_cols' %in% data_tests) {
      extra.columns <- setdiff(names(dat), names(expected.data))
      if(length(extra.columns) > 0) {
        stop(paste0(data_type, " should omit these extra columns: ", paste0(extra.columns, collapse=", ")), call.=FALSE)
      }
    }
    
    # check for NA timestamps, better to run after missing_cols so the best 
    # error can be thrown if timecol is missing. check here, too, in case 
    # missing_cols was not among the data_tests or the metab_model data were 
    # specified without a timestamp column
    if('na_times' %in% data_tests) {
      # match against the known timestamp column names rather than a
      # substring grep for 'date'/'time', which would also match non-
      # timestamp columns such as 'travel.time'
      timecol <- intersect(c('solar.time','date'), names(dat))
      if(length(timecol) != 1) stop("in ", data_type, " found ", length(timecol), " possible timestamp columns", call.=FALSE)
      na.times <- which(is.na(dat[[timecol]]))
      if(length(na.times) > 0) {
        stop(paste0(data_type, " has NA date stamps in these rows: ", paste0(na.times, collapse=", ")), call.=FALSE)
      }
      if(timecol=='solar.time' && !lubridate::is.POSIXct(dat[[timecol]])) stop("expecting 'solar.time' to be of class 'POSIXct'", call.=FALSE)
      if(timecol=='solar.time' && !(lubridate::tz(dat[[timecol]]) %in% c('UTC','GMT'))) stop("expecting 'solar.time' to have timezone 'UTC'", call.=FALSE)
      if(timecol=='date' && !lubridate::is.Date(dat[[timecol]])) stop("expecting 'date' to be of class 'Date'", call.=FALSE)
    }
    
    # put the data columns in the same order as expected.data and eliminate any 
    # extra columns. accommodate (don't try to include) missing columns, which
    # will necessarily be optional if missing_cols was tested above
    keeper.columns <- names(expected.data)[names(expected.data) %in% names(dat)]
    dat <- dat[keeper.columns]
    expected.data <- expected.data[keeper.columns]
    
    # check for units mismatches. column names will already match exactly.
    if('units' %in% data_tests) {
      mismatched.units <- which(get_units(expected.data) != get_units(dat))
      if(length(mismatched.units) > 0) {
        data.units <- get_units(dat)[mismatched.units]
        expected.units <- get_units(expected.data)[mismatched.units]
        stop(paste0("unexpected units in ", data_type, ": ", paste0(
          "(", 1:length(mismatched.units), ") ",
          names(data.units), " = ", data.units, ", expected ", expected.units,
          collapse="; ")), call.=FALSE)
      }
    }

    # check travel.time positivity, if present. only two-station models supply
    # this column, so this is a no-op for other model types: the column filter
    # above has already dropped travel.time from any dat whose metab_class
    # doesn't declare it. A zero or negative travel time makes the upstream
    # lag meaningless, so this is a structural error.
    #
    # No upper bound is enforced here. An upper bound on travel time is a
    # question about individual days, not about the dataset's structure, and
    # the right response is to drop those days rather than reject the whole
    # dataset -- something validation's fail-fast contract can't express. That
    # ceiling now lives in mm_align_2s() (see mm_lag_2s.R).
    if('travel.time' %in% names(dat)) {
      travel.time <- v(dat$travel.time)
      if(any(travel.time <= 0)) {
        stop('travel.time must be > 0', call.=FALSE)
      }
    }

    # return the data, whose columns may be reordered/filtered
    dat
  })

  # return the data.frames, which may have had their columns reordered during validation and are packaged as a list
  return(dat_all)
}


#' Two-station-specific data validation
#'
#' Checks the lead-in existence requirement specific to
#' \code{\link{metab_bayes_2s}}: at least one row must have enough preceding
#' upstream observations to cover its own travel time, given the (median)
#' timestep of \code{data$solar.time}. This is a structural question -- is
#' there any data at all to lag from? -- and so belongs with validation's
#' other fail-fast checks.
#'
#' Rows that individually lack lead-in are not an error; they are dropped as
#' lead-in rows by \code{mm_align_2s} (see \code{mm_lag_2s.R}), which also
#' owns the per-day travel-time ceiling. Column presence, timestamp validity,
#' and travel.time positivity are expected to have already been checked by
#' \code{\link{mm_validate_data}}.
#'
#' @param data data.frame as returned by \code{\link{mm_validate_data}} for
#'   \code{\link{metab_bayes_2s}}: must contain \code{solar.time} and
#'   \code{travel.time}, sorted ascending by \code{solar.time}.
#' @keywords internal
mm_validate_data_2station <- function(data) {

  data_v <- v(data)

  # shared with prepdata_bayes_2s() and metab_bayes_2s() via mm_lag_2s(), so
  # the lag computed here can't drift from the one they apply
  lagged <- mm_lag_2s(data_v$solar.time, data_v$travel.time)

  if(!any(lagged$has_leadin)) {
    max_lag <- max(lagged$lag)
    stop(paste0(
      'insufficient lead-in data for upstream DO: no row has enough preceding upstream ',
      'observations to cover its own travel.time. The longest travel.time implies a lag of ',
      max_lag, ' timestep(s), but only ', nrow(data), ' row(s) were supplied; ',
      'supply lead-in upstream data before the first row to be modeled'))
  }

  invisible(NULL)
}
