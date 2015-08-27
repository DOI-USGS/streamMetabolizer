#' A prototype for the model_fun argument to mm_model_by_ply
#' 
#' This function does nothing but has the proper form for a model_fun passed to 
#' mm_model_by_ply. Other functions to be used as model_fun may call 
#' inheritParams to use the parameter definitions given here.
#' 
#' @param data_ply a data.frame containing all relevant, validated modeling data
#'   for a single ply of data. (1 ply ~= 1 date, although the day length has
#'   been specified by day_start and day_end and may not be exactly 24 hours)
#' @param data_daily_ply NULL or a data.frame containing inputs with a daily 
#'   timestep. The local.date column must be present.
#' @param ... other args that were passed untouched from the function calling 
#'   mm_model_by_ply, through mm_model_by_ply, and finally to this function.
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
#' @param local_date the modal date of this ply of data and data_daily, and the 
#'   date by which this ply should be referred
#' @export
mm_model_by_ply_prototype <- function(data_ply, data_daily_ply, ..., day_start, day_end, local_date) {
  message("this function does nothing and should only be called for testing")
}