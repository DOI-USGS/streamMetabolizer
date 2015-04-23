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
#' @return An empty metab_model.
#'   
#' @examples
#' metab_model() 
#' @export
metab_model <- function(
  fit=lm(GPP~DO, data=data.frame(GPP=1:5, DO=1:5+0.1)), # trivial and wrong
  args=list(),
  data=data.frame(GPP=1:5, DO=1:5+0.1),
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
    cat("metab_model of type ",class(object),"\n")
  }
)


#' Retrieve the arguments that were used to fit the metab_model
#' 
#' @inheritParams get_args
#' @export 
#' @family get_metadata
get_args.metab_model <- function(metab_model) {
  metab_model@args
}


#' Retrieve the data that were used to fit the model
#' 
#' @inheritParams get_fitting_data
#' @export
#' @family get_fitting_data
get_fitting_data.metab_model <- function(metab_model) {
  metab_model@data
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
  data.frame(GPP=GPP, ER=ER, NEP=GPP-ER)
}
