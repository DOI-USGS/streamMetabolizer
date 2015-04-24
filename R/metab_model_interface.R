#' Functions implemented by any \code{streamMetabolizer}-compatible metabolism 
#' model.
#' 
#' Metabolism models in the \code{streamMetabolizer} package all implement a 
#' common set of core functions. These functions are conceptually packaged as 
#' the \code{metab_model_interface} defined here.
#' 
#' @section Functions in the interface:
#'   
#'   \itemize{
#'   
#'   \item \code{\link{show}(metab_model) \{ display(metab_model) \}}
#'   
#'   \item \code{\link{get_fit}(metab_model) \{ return(fitted.model) \}}
#'   
#'   \item \code{\link{get_args}(metab_model) \{ return(args.list) \}}
#'   
#'   \item \code{\link{get_data}(metab_model) \{ return(data.frame) \}}
#'   
#'   \item \code{\link{get_version}(metab_model) \{ return(version.string) \}}
#'   
#'   \item \code{\link{predict_metab}(metab_model, ...) \{ return(data.frame)
#'   \}}
#'   
#'   \item \code{\link{predict_DO}(metab_model, ...) \{ return(data.frame) \}}
#'   
#'   }
#'   
#' @name metab_model_interface
#' @rdname metab_model_interface
#' @docType data
#' @format A collection of functions which any metabolism model in 
#'   \code{streamMetabolizer} should implement.
NULL

#### show ####
# show() is already a generic S4 function.


#### S3 generics ####

#' Extract the internal model from a metabolism model.
#' 
#' A function in the metab_model_interface. Returns the internal model 
#' representation as fitted to the supplied data and arguments.
#' 
#' @param metab_model A metabolism model, implementing the 
#'   metab_model_interface, for which to return the data
#' @return An internal model representation; may have any class
#' @export
#' @family metab_model_interface
#' @family get_fit
get_fit <- function(metab_model) {
  UseMethod("get_fit")
}


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
#' @family get_data
get_data <- function(metab_model) {
  UseMethod("get_data")
}


#' Extract the version of streamMetabolizer that was used to fit the model.
#' 
#' A function in the metab_model_interface. Returns the version of
#' streamMetabolizer that was used to fit the model.
#' 
#' @param metab_model A metabolism model, implementing the 
#'   metab_model_interface, for which to return the data
#' @return character representation of the package version
#' @export
#' @family metab_model_interface
#' @family get_version
get_version <- function(metab_model) {
  UseMethod("get_version")
}


#' Predict metabolism from a fitted model.
#' 
#' A function in the metab_model_interface. Returns predictions (estimates) of 
#' GPP, ER, and NEP.
#' 
#' @param metab_model A metabolism model, implementing the 
#'   metab_model_interface, to use in predicting metabolism
#' @return A data.frame of daily metabolism estimates. Columns: \describe{
#'   
#'   \item{GPP}{numeric estimate of Gross Primary Production, always positive,
#'   \eqn{mg O_2 L^{-1} d^{-1}}{mg O2 / L / d}}
#'   
#'   \item{R}{numeric estimate of Respiration, always negative, \eqn{mg O_2
#'   L^{-1} d^{-1}}{mg O2 / L / d}}
#'   
#'   \item{NEP}{numeric estimate of Net Ecosystem Production, positive for net
#'   autotrophy, \eqn{mg O_2 L^{-1} d^{-1}}{mg O2 / L / d}}
#'   
#'   }
#' @export
#' @family metab_model_interface
#' @family predict_metab
predict_metab <- function(metab_model) {
  UseMethod("predict_metab")
}


#' Predict DO from a fitted model.
#' 
#' A function in the metab_model_interface. Returns predictions of dissolved 
#' oxygen.
#' 
#' @param metab_model A metabolism model, implementing the 
#'   metab_model_interface, to use in predicting metabolism
#' @return A data.frame of dissolved oxygen predictions at the temporal
#'   resolution of the input data
#' @export
#' @family metab_model_interface
#' @family predict_DO
predict_DO <- function(metab_model) {
  UseMethod("predict_DO")
}