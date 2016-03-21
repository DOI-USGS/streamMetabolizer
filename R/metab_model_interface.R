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
#'   \item \code{\link{get_fitting_time}(metab_model) \{ return(proc_time) \}}
#'   
#'   \item \code{\link{get_info}(metab_model) \{ return(info) \}}
#'   
#'   \item \code{\link{get_specs}(metab_model) \{ return(specs.list) \}}
#'   
#'   \item \code{\link{get_data}(metab_model) \{ return(data.frame) \}}
#'   
#'   \item \code{\link{get_data_daily}(metab_model) \{ return(data.frame) \}}
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
#' @examples
#' methods(class="metab_model")
NULL

#### show ####
# show() is already a generic S4 function.


#### S3 generics ####

#' Extract the user-supplied metadata about a metabolism model.
#' 
#' A function in the metab_model_interface. Returns any user-supplied metadata.
#' 
#' @param metab_model A metabolism model, implementing the metab_model_interface, for which to return the metadata information.
#' @return The user-supplied metadata in the original format.
#' @export
#' @family metab_model_interface
#' @family get_info
get_info <- function(metab_model) {
  UseMethod("get_info")
}

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

#' Extract the amount of time that was required to fit the metabolism model.
#' 
#' A function in the metab_model_interface. Returns the time that was taken to
#' fit the model; see \code{\link{proc.time}} for details.
#' 
#' @param metab_model A metabolism model, implementing the 
#'   metab_model_interface, for which to return the time
#' @return An proc_time object
#' @export
#' @family metab_model_interface
#' @family get_fitting_time
get_fitting_time <- function(metab_model) {
  UseMethod("get_fitting_time")
}

#' Extract the fitting specifications from a metabolism model.
#' 
#' A function in the metab_model_interface. Returns the specifications that were
#' passed in when fitting the metabolism model.
#' 
#' @param metab_model A metabolism model, implementing the 
#'   metab_model_interface, for which to return the specifications
#' @return The list of specifications that was passed to \code{\link{metab}()}
#' @export
#' @family metab_model_interface
#' @family get_args
get_specs <- function(metab_model) {
  UseMethod("get_specs")
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

#' Extract the daily fitting data, if any, from a metabolism model.
#' 
#' A function in the metab_model_interface. Returns the daily data that were
#' passed to a metabolism model.
#' 
#' @param metab_model A metabolism model, implementing the 
#'   metab_model_interface, for which to return the data_daily
#' @return A data.frame
#' @export
#' @family metab_model_interface
#' @family get_data_daily
get_data_daily <- function(metab_model) {
  UseMethod("get_data_daily")
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
#' @param date_start Date or a class convertible with as.Date. The first date
#'   (inclusive) for which to report metabolism predictions. If NA, no filtering is 
#'   done.
#' @param date_end Date or a class convertible with as.Date. The last date 
#'   (inclusive) for which to report metabolism predictions. If NA, no filtering is 
#'   done.
#' @param ... Other arguments passed to class-specific implementations of
#'   \code{predict_metab}
#' @param use_saved logical. Is it OK to use predictions that were saved with
#'   the model?
#' @return A data.frame of daily metabolism estimates. Columns include:
#'   \describe{
#'   
#'   \item{GPP}{numeric estimate of Gross Primary Production, positive when 
#'   realistic, \eqn{mg O_2 L^{-1} d^{-1}}{mg O2 / L / d}}
#'   
#'   \item{ER}{numeric estimate of Ecosystem Respiration, negative when realistic, \eqn{mg
#'   O_2 L^{-1} d^{-1}}{mg O2 / L / d}}
#'   
#'   \item{K600}{numeric estimate of the reaeration rate \eqn{d^{-1}}{1 / d}}
#'   }
#' @examples 
#' dat <- data_metab('3', day_start=12, day_end=36)
#' mm <- metab_night(specs(mm_name('night')), data=dat)
#' predict_metab(mm)
#' predict_metab(mm, date_start=get_fit(mm)$date[2])
#' @export
#' @family metab_model_interface
#' @family predict_metab
predict_metab <- function(metab_model, date_start=NA, date_end=NA, ..., use_saved) {
  UseMethod("predict_metab")
}


#' Predict DO from a fitted model.
#' 
#' A function in the metab_model_interface. Returns predictions of dissolved 
#' oxygen.
#' 
#' @param metab_model A metabolism model, implementing the 
#'   metab_model_interface, to use in predicting metabolism
#' @param date_start Date or a class convertible with as.Date. The first date 
#'   (inclusive) for which to report DO predictions. If NA, no filtering is 
#'   done.
#' @param date_end Date or a class convertible with as.Date. The last date 
#'   (inclusive) for which to report DO predictions. If NA, no filtering is 
#'   done.
#' @param ... Other arguments passed to class-specific implementations of 
#'   \code{predict_DO}
#' @param use_saved logical. Is it OK to use predictions that were saved with 
#'   the model?
#' @return A data.frame of dissolved oxygen predictions at the temporal 
#'   resolution of the input data
#' @examples 
#' dat <- data_metab('3', day_start=12, day_end=36)
#' mm <- metab_night(specs(mm_name('night')), data=dat)
#' preds <- predict_DO(mm, date_start=get_fit(mm)$date[3])
#' head(preds)
#' @export
#' @family metab_model_interface
#' @family predict_DO
predict_DO <- function(metab_model, date_start=NA, date_end=NA, ..., use_saved) {
  UseMethod("predict_DO")
}