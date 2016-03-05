#' Fit a metabolism model to data
#' 
#' Runs the metabolism model specified by the \code{specs} argument. 
#' Returns a fitted model.
#' 
#' @author Alison Appling
#'   
#' @param specs a list of model specifications and parameters for a model.
#'   Although this may be specified manually, it is easier and safer to use
#'   \code{specs} to generate the list. The help file for that functions lists
#'   the necessary parameters, describes them in detail, and gives default
#'   values.
#' @param data data.frame of input data at the temporal resolution of raw 
#'   observations. Columns must have the same names, units, and format as the 
#'   default. See \code{\link{mm_data}} for a full data description.
#' @param data_daily data.frame containing inputs with a daily timestep. See 
#'   \code{\link{mm_data}} for a full data description.
#' @param info any information, in any format, that you would like to store 
#'   within the metab_model object
#' @export
#' @return An object inheriting from metab_model and containing the fitted 
#'   model. This object can be inspected with the functions in the 
#'   \code{\link{metab_model_interface}}.
#' @examples
#' 
#' @export
metab <- function(specs=specs(mm_name()), data=mm_data(NULL), data_daily=mm_data(NULL), info=NULL) {

  # determine which model function to call
  model_type <- mm_parse_name(specs$model_name)$type
  metab_fun <- switch(
    model_type,
    bayes  = metab_bayes,
    Kmodel = metab_Kmodel,
    mle    = metab_mle,
    night  = metab_night,
    sim    = metab_sim)

  # run the model
  metab_fun(specs=specs, data=data, data_daily=data_daily, info=info)
}
