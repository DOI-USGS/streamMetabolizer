#### class definition ####

#' A metabolism model class.
#' 
#' @slot info Any metadata the user chooses to package with metabolism model.
#' @slot fit An internal representation of a fitted model.
#' @slot pkg_version A string indicating the package version used to create this metab_model object.
#' @slot args A list of arguments, excluding data, that were supplied to the fitting function.
#' @slot data The data that were used to fit the model.
#' @slot data_daily The daily data, if any, that were used to fit the model.
#' @exportClass metab_model
#' @family metab.model.classes
setClass(
  "metab_model",
  slots=c(
    info="ANY",
    fit="ANY",
    args="list",
    data="data.frame",
    data_daily="data.frame",
    pkg_version="character"),
  
  prototype=c(
    info=NULL,
    fit=NULL,
    args=NULL,
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
#' @param model_class character name of a class inheriting from metab_model - the
#'   type of object to create
#' @param info User-supplied metadata of any form.
#' @param fit An internal representation of a fitted model.
#' @param args A list of arguments, excluding data, that were supplied to the 
#'   fitting function.
#' @param data The data that were used to fit the model.
#' @param data_daily The data_daily that were used to fit the model. May be NULL.
#' @param pkg_version A string indicating the package version used to create 
#'   this metab_model object. The default should almost always be appropriate.
#' @param ... other arguments passed to new() for this particular model_class
#' @return A metab_model object.
#'   
#' @examples
#' metab_model() 
#' metab_model(fit=1:5, args=list(length=5))
#' metab_model("metab_mle", fit=1:5, args=list(length=5))
#' @importFrom utils packageVersion
#' @export
metab_model <- function(
  model_class="metab_model",
  info="user metadata goes here",
  fit="generic metab_model class; no actual fit",
  args=list(day_start=-1.5, day_end=30, calc_DO_fun=calc_DO_mod),
  data=mm_data(local.time, DO.obs, DO.sat, depth, temp.water, light),
  data_daily=mm_data(local.date, K600, discharge.daily, velocity.daily, optional="all"),
  pkg_version=as.character(packageVersion("streamMetabolizer")),
  ...) {
  
  # Create a dummy metab_model object
  new(model_class, info=info, fit=fit, args=args, data=data, pkg_version=pkg_version, ...)
}

#### loadModelInterface ####


#' Display the metab_model object
#' 
#' Print a metab_model object to the console.
#' 
#' @param object metab_model to be displayed.
#' @importFrom utils head
setMethod(
  "show", "metab_model", 
  function(object) {
    cat("metab_model", if(class(object)[1] != "metab_model") paste0("of type ", class(object)[1]), "\n")
    cat("  User-supplied metadata (access with get_info()):\n")
    print(object@info)
    cat("  Fitted model (access with get_fit()):\n")
    cat("    class: ", paste0(class(object@fit), collapse=" & "),"\n")
    cat("  Fitting arguments (access with get_args()):\n")
    for(arg in names(object@args)) {
      arg_char <- tryCatch(
        if(is.null(object@args[[arg]])) {
          'NULL' 
        } else {
          paste0(as.character(object@args[[arg]]), collapse=", ")
        }, 
        error=function(e) paste0(paste0(class(object@args[[arg]]), collapse=","),". see get_args(...)[['",arg,"']]"))
      if(nchar(arg_char) > 100) arg_char <- paste0(substr(arg_char, 1, 100), "...")
      cat(paste0("    ", arg, ": ", arg_char, "\n"))
    }
    cat("  Fitting data (truncated; access with get_data()):\n")
    print(head(object@data))
    cat("  Created with streamMetabolizer version", object@pkg_version, "\n")
  }
)


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
#' @importFrom stats qnorm setNames
#' @export
#' @family predict_metab
predict_metab.metab_model <- function(metab_model, ci_level=0.95, ...) {
  
  fit <- get_fit(metab_model)
  vars <- c("GPP","ER","K600")
  if(all(vars %in% names(fit))) {
    crit <- qnorm((1 + ci_level)/2)
    c(list(fit['local.date']),
      lapply(vars, function(var) {
        est <- fit[[var]]
        sd <- fit[[paste0(var,".sd")]]
        data.frame(
          est,
          lower = est - crit * sd,
          upper = est + crit * sd) %>% 
          setNames(c(var, paste0(var, ".", c("lower","upper"))))
      })) %>%
      bind_cols() %>%
      as.data.frame()
  } else {
    warning("model does not contain all columns ", paste0(vars, collapse=", "))
    data.frame(
      local.date=as.Date(NA)[NULL], 
      GPP=numeric(), GPP.lower=numeric(), GPP.upper=numeric(),
      ER=numeric(), ER.lower=numeric(), ER.upper=numeric(),
      K600=numeric(), K600.lower=numeric(), K600.upper=numeric()
    )
  }
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
  
  # pull args from the model
  calc_DO_fun <- get_args(metab_model)$calc_DO_fun
  day_start <- get_args(metab_model)$day_start
  day_end <- get_args(metab_model)$day_end
  
  # get the metabolism (GPP, ER) data and estimates
  metab_ests <- predict_metab(metab_model)
  data <- get_data(metab_model)
  
  # re-process the input data with the metabolism estimates to predict DO
  mm_model_by_ply(
    mm_predict_1ply, data=data, data_daily=metab_ests, # for mm_model_by_ply
    day_start=day_start, day_end=day_end, # for mm_model_by_ply
    calc_DO_fun=calc_DO_fun) # for mm_predict_1ply
  
}
