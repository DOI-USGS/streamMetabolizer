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
  
  # print the fitting time if requested
  if(isTRUE(specs$verbose)) {
    fitting_mins <- floor(fitting_time[['elapsed']] / 60)
    fitting_secs <- round(fitting_time[['elapsed']] %% 60)
    message("model fit in ", if(fitting_mins > 0) paste0(fitting_mins, " min, "), fitting_secs, " sec")
  }
  
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
    if(!is.null(get_info(object))) {
      cat("  User-supplied metadata:\n")
      print(get_info(object))
    }
    cat("streamMetabolizer version", object@pkg_version, "\n")
    print(get_specs(object), header="Specifications:\n", prefix="  ")
    cat("Fitting time: ", get_fitting_time(object)[['elapsed']], " secs elapsed\n", sep="")
    withCallingHandlers(
      tryCatch({
        metab_preds <- predict_metab(object)
        cat("Predictions (", nrow(metab_preds), " date", if(nrow(metab_preds)!=1) "s", "):\n", sep='')
        print(head(metab_preds, 10))
        if(nrow(metab_preds) > 10) cat("  ...")
      }, error=function(err) {
        cat("Prediction errors:\n")
        cat(paste0("  ", err$message))
      }), warning=function(war) {
        cat("Prediction warnings:", paste0('  ', war$message, collapse='\n'), sep='\n')
        invokeRestart("muffleWarning")
      })
    if(!is.data.frame(object@fit)) {
      if(('warnings' %in% names(object@fit)) && length(object@fit$warnings) > 0)
        cat("General fitting warnings:", paste0('  ', object@fit$warnings, collapse='\n'), sep='\n')
      if(('errors' %in% names(object@fit)) && length(object@fit$errors) > 0)
        cat("General fitting errors:", paste0('  ', object@fit$errors, collapse='\n'), sep='\n')
    }
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
#' @param attach.units logical. Should units be attached to the output?
#' @return A data.frame of predictions, as for the generic 
#'   \code{\link{predict_metab}}.
#' @importFrom stats qnorm setNames
#' @import dplyr
#' @importFrom unitted u v get_units
#' @export
#' @family predict_metab
predict_metab.metab_model <- function(metab_model, date_start=NA, date_end=NA, ..., attach.units=FALSE) {
  
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
    preds <- fit[c('date', precalc_cis, 'warnings', 'errors')] %>% 
      setNames(c('date', paste0(rep(var_vec, each=length(var_vec)), rep(c('','.lower','.upper'), times=3)), 'warnings', 'errors'))

  } else if(length(calcnow_cis) == 6) {
    # the fit includes columns for GPP, GPP.sd, ER, ER.sd, K600, and K600.sd; calculate CIs
    ci_level <- 0.95
    crit <- qnorm((1 + ci_level)/2)
    preds <- c(list(fit['date']),
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
      as.data.frame(stringsAsFactors=FALSE)
    
  } else {
    warning("model is missing columns for estimates and/or CIs")
    preds <- data.frame(
      date=as.Date(NA)[NULL], 
      GPP=numeric(), GPP.lower=numeric(), GPP.upper=numeric(),
      ER=numeric(), ER.lower=numeric(), ER.upper=numeric(),
      K600=numeric(), K600.lower=numeric(), K600.upper=numeric(),
      warnings=character(), errors=character(),
      stringsAsFactors=FALSE
    )
  }
  
  if(is.list(fit)) {
    if(!exists('valid_day', fit)) fit$valid_day <- TRUE
    attr(preds, 'warnings') <- format_stopwarn_msgs(fit$warnings[fit$valid_day])
    attr(preds, 'errors') <- format_stopwarn_msgs(fit$errors[fit$valid_day])
  }
  class(preds) <- c('preds_metab', class(preds))
  if(attach.units) {
    p <- u(preds, get_units(mm_data())[sapply(names(preds), function(x) strsplit(x, '\\.')[[1]][1], USE.NAMES=FALSE)])
  }
  preds
}

#' Summarize a vector of warning or error messages
#' 
#' Split ;-separated warning/error messages and condense into counts of each
#' unique message
#' 
#' @keywords internal
format_stopwarn_msgs <- function(msgs) {
  split_msgs <- unlist(strsplit(msgs, '; '))
  tbl_msgs <- sort(table(split_msgs))
  if(length(tbl_msgs) == 0) {
    c()
  } else {
    msgs_w_counts <- paste0(
      unname(tbl_msgs), " date", ifelse(unname(tbl_msgs)==1,"","s"), ": ",
      names(tbl_msgs))
    sort(msgs_w_counts)
  }
}

# Define unitted_preds_metab. Pretty limited functionality; hopefully enough for
# typical user needs
setOldClass('preds_metab')
setClass('unitted_preds_metab', contains=c('preds_metab','unitted_data.frame'))
setMethod("unitted", 'preds_metab', function(object, units) {
  if(!is.data.frame(object)) stop('expected object to be a data.frame')
  class(object) <- class(object)[class(object) != 'preds_metab']
  x <- u(object, units)
  as(x, 'unitted_preds_metab')
})
setMethod("show", 'unitted_preds_metab', function(object) {
  print.preds_metab(object)
})
setMethod("deunitted", 'unitted_preds_metab', function(object) {
  x <- v(as(x, 'unitted_data.frame'))
  class(x) <- c('preds_metab', class(x))
  x
})

#' Print metab preds, warnings, & errors
#' 
#' Print metab predictions with their warnings & errors if applicable
#' @inheritParams base::print
#' @export
print.preds_metab <- function(x, ...) {
  pristine_x <- x
  class(x) <- class(x)[class(x) != 'preds_metab']
  print(x)
  if(!is.null(attr(x, 'warnings')))
    cat("Summary of date-specific fitting warnings:", paste0('  ', attr(x, 'warnings'), collapse='\n'), sep='\n')
  if(!is.null(attr(x, 'errors')))
    cat("Summary of date-specific fitting errors:", paste0('  ', attr(x, 'errors'), collapse='\n'), sep='\n')
  invisible(pristine_x)
}


#' Make dissolved oxygen predictions from a fitted metab_model.
#' 
#' Makes fine-scale predictions of dissolved oxygen using fitted coefficients, 
#' etc. from the metabolism model.
#' 
#' @inheritParams predict_DO
#' @param attach.units logical. Should units be attached to the output?
#' @return A data.frame of predictions, as for the generic 
#'   \code{\link{predict_DO}}.
#' @import dplyr
#' @importFrom lubridate tz
#' @importFrom unitted u v get_units
#' @export
#' @family predict_DO
predict_DO.metab_model <- function(metab_model, date_start=NA, date_end=NA, ..., attach.units=FALSE, use_saved=TRUE) {
  
  # pull args from the model
  specs <- get_specs(metab_model)
  day_start <- specs$day_start
  day_end <- specs$day_end
  
  # get the input data; filter if requested
  data <- get_data(metab_model) %>% 
    v() %>% # units are discarded here. someday we'll get them worked through the whole system.
    mm_filter_dates(date_start=date_start, date_end=date_end, day_start=day_start, day_end=day_end)
  
  # if allowed and available, use previously stored values for DO.mod rather than re-predicting them now
  if(isTRUE(use_saved)) {
    if(!is.null(data) && "DO.mod" %in% names(data)) {
      return(data)
    }
  }
  
  # get the metabolism estimates; filter as we did for data
  warnings <- errors <- '.dplyr.var'
  metab_ests <- predict_metab(metab_model, date_start=date_start, date_end=date_end) %>%
    select(-warnings, -errors)
    
  # re-process the input data with the metabolism estimates to predict DO
  preds <- mm_model_by_ply(
    mm_predict_1ply, data=data, data_daily=metab_ests, # for mm_model_by_ply
    day_start=day_start, day_end=day_end, day_tests=c(), timestep_days=FALSE, # for mm_model_by_ply
    model_name=specs$model_name) %>% # for mm_predict_1ply
    mm_filter_dates(date_start=date_start, date_end=date_end) # trim off the extra
  
  if(attach.units) preds <- u(preds, get_units(mm_data())[gsub('DO.mod','DO.obs',names(preds))])
  preds
}
