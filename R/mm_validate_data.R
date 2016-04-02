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
        stop(paste0(data_type, " is NULL but required"))
      }
    }
    
    # check for missing or extra columns
    if('missing_cols' %in% data_tests) {
      missing.columns <- setdiff(names(expected.data), names(dat))
      missing.columns <- setdiff(missing.columns, optional.data) # optional cols don't count
      if(length(missing.columns) > 0) {
        stop(paste0(data_type, " is missing these columns: ", paste0(missing.columns, collapse=", ")))
      }
    }
    if('extra_cols' %in% data_tests) {
      extra.columns <- setdiff(names(dat), names(expected.data))
      if(length(extra.columns) > 0) {
        stop(paste0(data_type, " should omit these extra columns: ", paste0(extra.columns, collapse=", ")))
      }
    }
    
    # check for NA timestamps, better to run after missing_cols so the best 
    # error can be thrown if timecol is missing. check here, too, in case 
    # missing_cols was not among the data_tests or the metab_model data were 
    # specified without a timestamp column
    if('na_times' %in% data_tests) {
      timecol <- grep('date|time', names(dat), value=TRUE)
      if(length(timecol) != 1) stop("in ", data_type, " found ", length(timecol), " possible timestamp columns")
      na.times <- which(is.na(dat[,timecol]))
      if(length(na.times) > 0) {
        stop(paste0(data_type, " has NA date stamps in these rows: ", paste0(na.times, collapse=", ")))
      }
      if(timecol=='solar.time' && !lubridate::is.POSIXct(dat[,timecol])) stop("expecting 'solar.time' to be of class 'POSIXct'")
      if(timecol=='date' && !lubridate::is.Date(dat[,timecol])) stop("expecting 'date' to be of class 'Date'")
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
          collapse="; ")))
      }
    }
    
    # return the data, whose columns may be reordered/filtered
    dat
  })
  
  # return the data.frames, which may have had their columns reordered during validation and are packaged as a list
  return(dat_all)
}