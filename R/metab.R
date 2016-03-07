#' Fit a metabolism model to data
#' 
#' Runs the metabolism model specified by the \code{specs} argument. Returns a 
#' fitted model.
#' 
#' @author Alison Appling
#'   
#' @param specs a list of model specifications and parameters for a model. 
#'   Although this may be specified manually (it's just a list), it is easier
#'   and safer to use \code{\link{specs}} to generate the list, because the set
#'   of required parameters and their defaults depends on the model given in the
#'   \code{model_name} argument to \code{specs}. The help file for 
#'   \code{\link{specs}} lists the necessary parameters, describes them in 
#'   detail, and gives default values.
#' @param data data.frame of input data at the temporal resolution of raw 
#'   observations (unit-value). Columns must have the same names, units, and 
#'   format as the default. See the \strong{'Formatting \code{data}'} section 
#'   below for a full description.
#' @param data_daily data.frame containing inputs with a daily timestep. See the
#'   \strong{'Formatting \code{data_daily}'} section below for a full 
#'   description.
#' @param info any information, in any format, that you would like to store 
#'   within the metab_model object
#' @return An object inheriting from metab_model and containing the fitted 
#'   model. This object can be inspected with the functions in the 
#'   \code{\link{metab_model_interface}}.
#'   
#' @template metab_data
#'   
#' @examples
#' dat <- data_metab(num_days='3')
#' 
#' # fit a basic MLE model
#' mm <- metab(specs(mm_name('mle')), data=dat, info='my info')
#' predict_metab(mm)
#' get_info(mm)
#' get_fitting_time(mm)
#' 
#' # with chaining & customization
#' library(dplyr)
#' mm <- mm_name('mle', ode_method='Euler') %>%
#'   specs(GPP_init=40) %>%
#'   metab(data=dat)
#' predict_metab(mm)
#' \dontrun{
#' plot_DO_preds(predict_DO(mm))
#' plot_DO_preds(predict_DO(mm), y_var='pctsat', style='dygraphs')
#' }
#' @export
metab <- function(specs=specs(mm_name()), data=v(mm_data(NULL)), data_daily=v(mm_data(NULL)), info=NULL) {
  
  if(missing(specs)) {
    # if specs is left to the default, it gets confused about whether specs() is
    # the argument or the function. tell it which:
    specs <- streamMetabolizer::specs(mm_name())
  }
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
