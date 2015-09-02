#' A prototype for the fitting function of each metab_model type
#' 
#' This function does nothing but has the proper form for the constructor of an 
#' object inheriting from metab_model. Other constructor functions may call 
#' inheritParams to use the parameter definitions given here.
#' 
#' @param data data.frame of input data at the temporal resolution of raw 
#'   observations. Columns must have the same names, units, and format as the 
#'   default. See \code{\link{mm_data}} for a full data description.
#' @param data_daily data.frame containing inputs with a daily timestep. See 
#'   \code{\link{mm_data}} for a full data description.
#' @param info any information, in any format, that you would like to store 
#'   within the metab_model object
#' @param day_start start time of a day's data in number of hours from the 
#'   midnight that begins the modal date. For example, day_start=-1.5 indicates 
#'   that data describing 2006-06-26 begin at 2006-06-25 22:30, or at the first 
#'   observation time that occurs after that time if day_start doesn't fall 
#'   exactly on an observation time.
#' @param day_end end time of a day's data in number of hours from the midnight 
#'   that begins the modal date. For example, day_end=30 indicates that data 
#'   describing 2006-06-26 end at 2006-06-27 06:00, or at the last observation 
#'   time that occurs before that time if day_end doesn't fall exactly on an 
#'   observation time.
#' @param model_specs a list of model specifications and parameters for a model.
#'   Although this may be specified manually, it is easier to use a predefined 
#'   function from the \code{specs_xxx} family with a name beginning with 
#'   "specs_xxx", where "xxx" is the model type - e.g., "mle", "bayes", "night",
#'   etc. The help files for those functions list the necessary parameters,
#'   describe them in detail, and give default values.
metab_model_prototype <- function(data, data_daily, info, day_start, day_end, model_specs) {
  message("this function does nothing and should only be called for testing")
}