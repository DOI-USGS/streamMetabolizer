#' Functions implemented by any \code{streamMetabolizer}-compatible metabolism 
#' model.
#' 
#' Metabolism models in the \code{streamMetabolizer} package all implement a 
#' common set of core functions. These functions are conceptually packaged as 
#' the \code{metab_model_interface} defined here.
#' 
#' @section Functions in the interface: \itemize{ \item 
#'   \code{\link{show}(metab_model) \{ print(metab_model) \}} \item 
#'   \code{\link{get_args}(metab_model) \{ return(args.list) \}} \item 
#'   \code{\link{get_fitting_data}(metab_model) \{ return(data.frame) \}} \item 
#'   \code{\link{predict_metab}(metab_model, ...) \{ return(data.frame) \}} }
#'   
#' @name metab_model_interface
#' @rdname metab_model_interface
#' @docType data
#' @format A collection of functions which any metabolism model in 
#'   \code{streamMetabolizer} should implement.
NULL

#### show ####
# show() is already a generic S4 function.

#' Extract the fitting arguments from a metabolism model.
#' 
#' A function in the metab_model_interface. Returns the arguments that were 
#' passed to a metabolism model.
#' 
#' @param metab_model A metabolism model, implementing the
#'   metab_model_interface, for which to return the arguments
#' @return A list of arguments
#' @export
#' @family metab_model_interface
#' @family get_args
get_args <- function(metab_model) {
  UseMethod("get_args")
}


#' Extract the fitting data from a metabolism model.
#' 
#' A function in the metab_model_interface. Returns the data that were passed to
#' a metabolism model.
#' 
#' @param metab_model A metabolism model, implementing the
#'   metab_model_interface, for which to return the data
#' @return A data.frame
#' @export
#' @family metab_model_interface
#' @family get_fitting_data
get_fitting_data <- function(metab_model) {
  UseMethod("get_fitting_data")
}

#' Predict metabolism from a fitted model.
#' 
#' A function in the metab_model_interface. Returns predictions of GPP, ER, and 
#' NEP.
#' 
#' @param metab_model A metabolism model, implementing the 
#'   metab_model_interface, to use in predicting metabolism
#' @return A list of arguments
#' @export
#' @family metab_model_interface
#' @family predict_metab
predict_metab <- function(metab_model) {
  UseMethod("predict_metab")
}