#' @include specs-class.R
NULL

#### class definition ####

#' A metabolism model class.
#' 
#' @slot info Any metadata the user chooses to package with metabolism model.
#' @slot fit An internal representation of a fitted model.
#' @slot fitting_time Usually stored as a proc_time; the time required to fit
#'   the model.
#' @slot pkg_version A string indicating the package version used to create this
#'   metab_model object.
#' @slot specs A list of model specifications that were supplied to the fitting
#'   function.
#' @slot data The data that were used to fit the model.
#' @slot data_daily The daily data, if any, that were used to fit the model.
#' @exportClass metab_model
#' @family metab.model.classes
setClass(
  "metab_model",
  slots=c(
    info="ANY",
    metab_daily="ANY",
    fit="ANY",
    fitting_time="ANY",
    specs="specs",
    data="ANY",
    data_daily="ANY",
    pkg_version="character"),
  
  prototype=c(
    info=NULL,
    metab_daily=NULL,
    fit=NULL,
    fitting_time=NULL,
    specs=NULL,
    data=NULL,
    data_daily=NULL,
    pkg_version=""),
  
  # returns TRUE if valid, vector of error strings otherwise
  validity=function(object) {
    errorstrs <- character()
    
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
#' @param model_class character name of a class inheriting from metab_model - 
#'   the type of object to create
#' @param info User-supplied metadata of any form.
#' @param metab_daily a data.frame of daily metabolism estimates produced from
#'   the fit.
#' @param fit An internal representation of a fitted model.
#' @param fitting_time A proc_time object giving the time taken to fit the 
#'   model.
#' @param specs A list of model specifications that were supplied to the fitting
#'   function.
#' @param data The data that were used to fit the model.
#' @param data_daily The data_daily that were used to fit the model. May be 
#'   NULL.
#' @param pkg_version A string indicating the package version used to create 
#'   this metab_model object. The default should almost always be appropriate.
#' @param ... other arguments passed to new() for this particular model_class
#' @return A metab_model object.
#'   
#' @examples
#' metab_model() 
#' metab_model("metab_mle", fit=1:5, specs=list(length=5))
#' @importFrom utils packageVersion
#' @export
metab_model <- function(
  model_class="metab_model",
  info="user metadata goes here",
  metab_daily=NULL,
  fit="generic metab_model class; no actual fit",
  fitting_time=system.time({}),
  specs=list(),
  data=mm_data(solar.time, DO.obs, DO.sat, depth, temp.water, light),
  data_daily=mm_data(date, optional="all"),
  pkg_version=as.character(packageVersion("streamMetabolizer")),
  ...) {
  
  # print the fitting time if requested
  if(isTRUE(specs$verbose)) {
    fitting_mins <- floor(fitting_time[['elapsed']] / 60)
    fitting_secs <- round(fitting_time[['elapsed']] %% 60)
    message("model fit in ", if(fitting_mins > 0) paste0(fitting_mins, " min, "), fitting_secs, " sec")
  }
  
  # Create a metab_model object
  new(model_class, info=info, metab_daily=metab_daily,
      fit=fit, fitting_time=fitting_time, 
      specs=if(class(specs)[1] == 'specs') specs else add_specs_class(specs), 
      data=data, data_daily=data_daily, pkg_version=pkg_version, ...)
}


#### loadModelInterface ####

# Complex methods have their own files, e.g., metab_model.predict_DO.R

#' Retrieve the user-supplied model metadata
#' 
#' @inheritParams get_info
#' @export 
#' @family get_info
get_info.metab_model <- function(metab_model) {
  metab_model@info
}

#' Retrieve the inner fitted model
#' 
#' @inheritParams get_fit
#' @export 
#' @family get_fit
get_fit.metab_model <- function(metab_model) {
  metab_model@fit
}

#' Retrieve the time it took to fit the model
#' 
#' @inheritParams get_fitting_time
#' @export 
#' @family get_fitting_time
get_fitting_time.metab_model <- function(metab_model) {
  metab_model@fitting_time
}

#' Retrieve the specifications that were used to fit the metab_model
#' 
#' @inheritParams get_specs
#' @export 
#' @family get_specs
get_specs.metab_model <- function(metab_model) {
  metab_model@specs
}

#' Retrieve the data that were used to fit the model
#' 
#' @inheritParams get_data
#' @export
#' @family get_data
get_data.metab_model <- function(metab_model) {
  metab_model@data
}

#' Retrieve the data_daily that were used to fit the model
#' 
#' @inheritParams get_data_daily
#' @export
#' @family get_data_daily
get_data_daily.metab_model <- function(metab_model) {
  tryCatch(
    metab_model@data_daily,
    error=function(e) {
      warning('this metab_model is out of date and has no data_daily slot')
      as.data.frame(NULL)
    }
  )
}

#' Retrieve the version of streamMetabolizer that was used to fit the model
#' 
#' @inheritParams get_version
#' @export
#' @family get_version
get_version.metab_model <- function(metab_model) {
  metab_model@pkg_version
}
