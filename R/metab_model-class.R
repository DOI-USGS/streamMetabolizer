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
    fit="ANY",
    fitting_time="ANY",
    specs="specs",
    data="ANY",
    data_daily="ANY",
    pkg_version="character"),
  
  prototype=c(
    info=NULL,
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
  fit="generic metab_model class; no actual fit",
  fitting_time=system.time({}),
  specs=list(),
  data=mm_data(solar.time, DO.obs, DO.sat, depth, temp.water, light),
  data_daily=mm_data(date, K600, discharge, velocity, optional="all"),
  pkg_version=as.character(packageVersion("streamMetabolizer")),
  ...) {
  
  # Create a metab_model object
  new(model_class, info=info, fit=fit, fitting_time=fitting_time, 
      specs=if(class(specs)[1] == 'specs') specs else add_specs_class(specs), 
      data=data, data_daily=data_daily, pkg_version=pkg_version, ...)
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
    print(get_info(object))
    cat("  Fitted model (access with get_fit()):\n")
    cat("    class: ", paste0(class(get_fit(object)), collapse=" & "),"\n")
    print(get_specs(object), header="  Fitting specifications (access with get_specs()):\n", prefix="    ")
    cat("  Fitting time (access with get_fitting_time()):\n")
    print(get_fitting_time(object))
    cat("  Fitting data (truncated; access with get_data()):\n")
    print(head(get_data(object)))
    cat("  Fitting data_daily (truncated; access with get_data_daily()):\n")
    print(head(get_data_daily(object)))
    cat("  streamMetabolizer version", object@pkg_version, "\n")
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


#' Make metabolism predictions from a fitted metab_model.
#' 
#' Makes daily predictions of GPP, ER, and K600 with upper and lower bounds
#' reflecting a 95% CI.
#' 
#' @inheritParams predict_metab
#' @return A data.frame of predictions, as for the generic 
#'   \code{\link{predict_metab}}.
#' @importFrom stats qnorm setNames
#' @import dplyr
#' @export
#' @family predict_metab
predict_metab.metab_model <- function(metab_model, date_start=NA, date_end=NA, ...) {
  
  fit <- get_fit(metab_model) %>%
    mm_filter_dates(date_start=date_start, date_end=date_end)
  
  var_vec <- c('GPP','ER','K600')
  precalc_cis <- grep("^(GPP|ER|K600)_daily_(50pct|2\\.5pct|97\\.5pct)$", names(fit), value=TRUE)
  calcnow_cis <- grep("^(GPP|ER|K600)(\\.sd)*$", names(fit), value=TRUE)
  if(length(precalc_cis) == 9) {
    # the fit includes columns for GPP_daily_2.5pct, etc.; use these as CIs. first get the columns in order
    vars <- medlohi <- '.dplyr.var'
    grepstrs <- expand.grid(
      medlohi=c('50pct', '2\\.5pct', '97\\.5pct'),
      vars=var_vec) %>%
      with(paste0("^(",vars,")_daily_(",medlohi,")$"))
    precalc_cis <- sapply(grepstrs, grep, names(fit), value=TRUE, USE.NAMES=FALSE)
    fit[c('date', precalc_cis, 'warnings', 'errors')] %>% 
      setNames(c('date', paste0(rep(c('GPP','ER','K600'), each=3), rep(c('','.lower','.upper'), times=3)), 'warnings', 'errors'))

  } else if(length(calcnow_cis) == 6) {
    # the fit includes columns for GPP, GPP.sd, ER, ER.sd, K600, and K600.sd; calculate CIs
    ci_level <- 0.95
    crit <- qnorm((1 + ci_level)/2)
    c(list(fit['date']),
      lapply(var_vec, function(var) {
        est <- fit[[var]]
        sd <- fit[[paste0(var,".sd")]]
        if(is.null(sd)) sd <- NA
        data.frame(
          est,
          lower = est - crit * sd,
          upper = est + crit * sd) %>% 
          setNames(c(var, paste0(var, ".", c("lower","upper"))))
      }),
      list(fit[c('warnings','errors')])) %>%
      bind_cols() %>%
      as.data.frame()
    
  } else {
    warning("model is missing columns for estimates and/or CIs")
    data.frame(
      date=as.Date(NA)[NULL], 
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
#' @param calc_DO_fun The function to use in predicting DO from metabolism 
#'   estimates. The default is generally recommended, even if the model was fit 
#'   with assumptions of process error.
#' @param calc_DO_args Any arguments to pass to calc_DO_fun. If not supplied, 
#'   the original args used when fitting the model will be used; this is the 
#'   recommended option.
#' @return A data.frame of predictions, as for the generic 
#'   \code{\link{predict_DO}}.
#' @import dplyr
#' @importFrom lubridate tz
#' @export
#' @family predict_DO
predict_DO.metab_model <- function(metab_model, date_start=NA, date_end=NA, calc_DO_fun=calc_DO_mod, calc_DO_args, ..., use_saved=TRUE) {
  
  # pull args from the model
  if(missing(calc_DO_args)) calc_DO_args <- get_specs(metab_model)$calc_DO_args # OK to be NULL
  day_start <- get_specs(metab_model)$day_start
  day_end <- get_specs(metab_model)$day_end
  
  # get the input data; filter if requested
  data <- get_data(metab_model) %>% 
    v() %>% # units are lost here. someday we'll get them worked through the whole system.
    mm_filter_dates(date_start=date_start, date_end=date_end, day_start=day_start, day_end=day_end)
  
  # if allowed and available, use previously stored values for DO.mod rather than re-predicting them now
  if(isTRUE(use_saved)) {
    if(!is.null(data) && "DO.mod" %in% names(data)) {
      return(data)
    }
  }
  
  # get the metabolism estimates; filter as we did for data
  metab_ests <- predict_metab(metab_model, date_start=date_start, date_end=date_end)
    
  # re-process the input data with the metabolism estimates to predict DO
  mm_model_by_ply(
    mm_predict_1ply, data=data, data_daily=metab_ests, # for mm_model_by_ply
    day_start=day_start, day_end=day_end, day_tests=c(), # for mm_model_by_ply
    calc_DO_fun=calc_DO_fun, calc_DO_args=calc_DO_args) %>% # for mm_predict_1ply
    mm_filter_dates(date_start=date_start, date_end=date_end) # trim off the extra
}
