#' Evaluate whether the data argument is properly formatted.
#' 
#' Will most often be called from within a metab_model constructor.
#' 
#' @param data a data.frame that has been (or will be) passed into a metab_model
#'   constructor
#' @param metab_class character the class name of the metab_model constructor
#' @import dplyr
#' @examples
#' \dontrun{
#' mm_validate_data(dplyr::select(mm_data(),-temp.water), "metab_mle")
#' }
#' @export
mm_validate_data <- function(data, data_daily, metab_class, tests=c('missing_cols','extra_cols','units')) {
  
  data_types <- setNames(c("data","data_daily"),c("data","data_daily"))
  dat_all <- lapply(data_types, function(data_type) {
    
    # pick out the data.frame for this loop
    dat <- get(data_type)
      
    # the data expectation is set by the default data argument to the specific metabolism class
    expected.data <- formals(metab_class)[[data_type]] %>% eval()
    optional.data <- attr(expected.data, 'optional')
    
    # quick return if dat is NULL
    if(is.null(v(dat))) {
      if('all' %in% optional.data) {
        return(dat)
      } else {
        stop(paste0(data_type, " is NULL but required"))
      }
    }
    
    # check for missing or extra columns
    if('missing_cols' %in% tests) {
      missing.columns <- setdiff(names(expected.data), names(dat))
      if(length(missing.columns) > 0) {
        stop(paste0(data_type, " is missing these columns: ", paste0(missing.columns, collapse=", ")))
      }
    }
    if('extra_cols' %in% tests) {
      extra.columns <- setdiff(names(dat), names(expected.data))
      if(length(extra.columns) > 0) {
        stop(paste0(data_type, " should omit these extra columns: ", paste0(extra.columns, collapse=", ")))
      }
    }
    
    # put the data columns in the same order as expected.data and eliminate any extra columns
    dat <- dat[names(expected.data)]
    
    # check for units mismatches. column names will already match exactly.
    if('units' %in% tests) {
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