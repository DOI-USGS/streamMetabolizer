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