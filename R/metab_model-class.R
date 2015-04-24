#### class definition ####

#' A generic metabolism model class.
#' 
#' Class and function definitions for a metabolism model
#' 
#' @rdname metab_model-class
#' @name metab_model-class
#' @slot fit An internal representation of a fitted model.
#' @slot pkg_version A string indicating the package version used to create this metab_model object.
#' @slot args A list of arguments, excluding data, that were supplied to the fitting function.
#' @slot data The data that were used to fit the model.
#' @exportClass metab_model
#' @family metab.model.classes
setClass(
  "metab_model",
  slots=c(
    fit="ANY",
    pkg_version="character",
    args="list",
    data="data.frame"),
  
  prototype=c(
    fit=NULL,
    pkg_version="",
    args=NULL,
    data=NULL),
  
  # returns TRUE if valid, vector of error strings otherwise
  validity=function(object) {
    errorstrs <- character()
    
    #     if(is.null(fit)) {
    #       errorstrs <- c(errorstrs, paste("fit shouldn't be null")))
    #     }
    
    # Return
    if(length(errorstrs) == 0) {
      TRUE
    } else {
      errorstrs
    }
  }
)

#### initialize ####

#' Create a metab_model object.
#' 
#' Generates a new model of class metab_model (\code{\link{metab_model-class}}).
#' 
#' @param fit An internal representation of a fitted model.
#' @param args A list of arguments, excluding data, that were supplied to the fitting function.
#' @param data The data that were used to fit the model.
#' @param pkg_version A string indicating the package version used to create this metab_model object.
#' @return A metab_model object.
#'   
#' @examples
#' metab_model() 
#' metab_model(fit=1:5, args=list(length=5))
#' @export
metab_model <- function(
  fit="generic metab_model class; no actual fit",
  args=list(),
  data=data.frame(date.time=as.POSIXct(NA), DO.obs=as.numeric(NA)),
  pkg_version=as.character(packageVersion("streamMetabolizer"))) {
  
  # Create a dummy metab_model object
  new("metab_model", fit=fit, args=args, data=data, pkg_version=pkg_version)
}

#### loadModelInterface ####


#' Display the metab_model object
#' 
#' Print a metab_model object to the console.
#' 
#' @param object metab_model to be displayed.
setMethod(
  "show", "metab_model", 
  function(object) {
    cat("metab_model", if(class(object)[1] != "metab_model") paste0(" of type ", class(object)[1]), "\n")
    cat("  Fitted model (access with get_fit()):\n")
    cat("    class: ", paste0(class(object@fit), collapse=" & "),"\n")
    cat("  Fitting arguments (access with get_args()):\n")
    for(arg in names(object@args)) {
      cat(paste0("    ", arg, ": ", object@args[[arg]], "\n"))
    }
    cat("  Fitting data (truncated; access with get_data()):\n")
    print(head(object@data))
    cat("  Created with streamMetabolizer version", object@pkg_version, "\n")
  }
)


#' Retrieve the inner fitted model
#' 
#' @inheritParams get_fit
#' @export 
#' @family get_fit
get_fit.metab_model <- function(metab_model) {
  metab_model@fit
}


#' Retrieve the arguments that were used to fit the metab_model
#' 
#' @inheritParams get_args
#' @export 
#' @family get_args
get_args.metab_model <- function(metab_model) {
  metab_model@args
}


#' Retrieve the data that were used to fit the model
#' 
#' @inheritParams get_data
#' @export
#' @family get_data
get_data.metab_model <- function(metab_model) {
  metab_model@data
}


#' Retrieve the version of streamMetabolizer that was used to fit the model
#' 
#' @inheritParams get_version
#' @export
#' @family get_version
get_version.metab_model <- function(metab_model) {
  metab_model@pkg_version
}


#' Make metabolism predictions from a fitted metab_model.
#' 
#' Makes daily predictions of GPP, ER, and NEP.
#' 
#' @inheritParams predict_metab
#' @return A data.frame of predictions, as for the generic 
#'   \code{\link{predict_metab}}.
#' @export
#' @family predict_metab
predict_metab.metab_model <- function(metab_model) {
  
  # Generate dummy output
  set.seed(3000)
  date <- as.Date(c("2015-04-15","2015-04-16","2015-04-17"))
  GPP <- abs(rnorm(3, 4, 1))
  ER <- -1.6*GPP
  data.frame(date=date, GPP=GPP, ER=ER, NEP=GPP-ER)
}



#' Make dissolved oxygen predictions from a fitted metab_model.
#' 
#' Makes fine-scale predictions of dissolved oxygen using fitted coefficients,
#' etc. from the metabolism model.
#' 
#' @inheritParams predict_DO
#' @return A data.frame of predictions, as for the generic 
#'   \code{\link{predict_DO}}.
#' @export
#' @family predict_DO
predict_DO.metab_model <- function(metab_model) {
  
  # Generate dummy output
  date.time <- do.call(seq, c(as.list(as.POSIXct(strptime(c("2015-04-16", "2015-04-17"), format="%Y-%m-%d"))), list(by=as.difftime(15, units="mins"))))[-97]
  DO <- sin(as.numeric(date.time-min(date.time))/86500*2*pi)
  data.frame(date.time=date.time, DO.mod=DO)
}
