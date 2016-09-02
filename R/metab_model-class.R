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


#' Display the metab_model object
#' 
#' Print a metab_model object to the console.
#' 
#' @param object metab_model to be displayed.
#' @importFrom utils head
#' @import dplyr
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
    
    # print fitting warnings if present
    params <- get_params(object, uncertainty='ci', fixed='stars', messages=TRUE)
    fixinfo <- if(any(grepl('\\*', params))) "(* = fixed value)" else ""
    cat("Parameters (", nrow(params), " date", if(nrow(params)!=1) "s", ")", fixinfo, ":\n", sep='')
    if(!exists('warnings', params)) params$warnings <- ''
    if(!exists('errors', params)) params$errors <- ''
    params %>%
      head(10) %>%
      mutate(messages = paste0(ifelse(errors != '', ifelse(warnings != '', 'e w', 'e  '), ifelse(warnings != '', '  w', '   ')))) %>%
      select(-warnings, -errors) %>%
      print()
    if(nrow(params) > 10) cat("  ...\n")
    fit <- get_fit(object)
    if(is.data.frame(fit) && all(c('warnings','errors') %in% names(fit))) {
      if(!exists('valid_day', fit)) fit$valid_day <- TRUE
      warn_msgs <- format_stopwarn_msgs(fit$warnings[fit$valid_day])
      stop_msgs <- format_stopwarn_msgs(fit$errors[fit$valid_day])
    } else if(is.list(fit)) {
      warn_msgs <- fit$warnings
      stop_msgs <- fit$errors
    } else {
      warn_msgs <- stop_msgs <- NULL
    }
    if(length(warn_msgs) > 0) cat("Fitting warnings:", paste0('  ', warn_msgs, collapse='\n'), sep='\n')
    if(length(stop_msgs) > 0) cat("Fitting errors:", paste0('  ', stop_msgs, collapse='\n'), sep='\n')
    
    # print metabolism predictions
    withCallingHandlers(
      tryCatch({
        metab_preds <- predict_metab(object)
        cat("Predictions (", nrow(metab_preds), " date", if(nrow(metab_preds)!=1) "s", "):\n", sep='')
        if(!exists('warnings', metab_preds)) metab_preds$warnings <- ''
        if(!exists('errors', metab_preds)) metab_preds$errors <- ''
        metab_preds %>%
          head(10) %>%
          mutate(messages = paste0(ifelse(errors != '', ifelse(warnings != '', 'e w', 'e  '), ifelse(warnings != '', '  w', '   ')))) %>%
          select(-warnings, -errors) %>%
          print()
        if(nrow(metab_preds) > 10) cat("  ...\n")
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

#' Collect the daily fitted parameters needed to predict GPP, ER, D, and DO
#' 
#' Returns a data.frame of parameters needed to predict GPP, ER, D, and DO
#' 
#' @inheritParams get_params
#' @return A data.frame of fitted parameters, as for the generic 
#'   \code{\link{get_params}}.
#' @export
#' @family get_params
get_params.metab_model <- function(
  metab_model, date_start=NA, date_end=NA, 
  uncertainty=c('sd','ci','none'), messages=TRUE, fixed=c('none','columns','stars'), 
  ..., attach.units=FALSE) {
  
  # process arguments
  uncertainty <- match.arg(uncertainty)
  fixed <- match.arg(fixed)
  
  # build the dDOdt function in order to pull out the metab.needs
  features <- mm_parse_name(get_specs(metab_model)$model_name)
  dDOdt <- create_calc_dDOdt(
    v(get_data(metab_model)[1,]), ode_method=features$ode_method, GPP_fun=features$GPP_fun,
    ER_fun=features$ER_fun, deficit_src=features$deficit_src)
  metab.needs <- environment(dDOdt)$metab.needs
  metab.optional <- c('DO.mod.1') # maybe should embed this in create_calc_DO?
  metab.all <- union(metab.needs, metab.optional)
  metab.search <- c(paste0(c('date','warnings','errors'),'$'), metab.all) %>%
    paste0('^', .) %>%
    paste0(collapse='|')
  
  # extract the daily parameters plus whatever else is daily (sds, gradients, etc.) 
  fit <- get_fit(metab_model)
  ddat <- get_data_daily(metab_model)
  
  # make sure we've got everything we need
  if(length(missing.metabs <- metab.needs[!metab.needs %in% union(names(fit), names(ddat))]) > 0) {
    stop(paste0("can't find metabolism parameter", if(length(missing.metabs)>1) "s", " ", paste0(missing.metabs, collapse=', ')))
  }
  
  # combine all daily values into one data.frame
  if(!is.null(fit) && !is.null(ddat) && nrow(ddat) > 0) {
    pars <- full_join(fit, ddat, by='date', copy=TRUE) 
  } else {
    if(!is.null(fit)) 
      pars <- fit
    else if(!is.null(ddat))
      pars <- ddat
    else
      return(NULL) # nothing available
  }
  pars <- pars %>%
    mm_filter_dates(date_start=date_start, date_end=date_end) %>% 
    { .[grep(metab.search, names(.), value=TRUE)] }
  
  # track provenance of each metab parameter. if any variables were available in
  # both x and y forms, combine them to minimize NAs
  metab.fit <- names(fit) %>% {.[. %in% metab.all]}
  metab.ddat <- names(ddat) %>% {.[. %in% metab.all]}
  metab.both <- intersect(metab.fit, metab.ddat)
  metab.either <- union(metab.fit, metab.ddat)
  for(a in metab.either) {
    if(a %in% metab.both) {
      a.x <- paste0(a,'.x')
      a.y <- paste0(a,'.y')
      pars[[a]] <- coalesce(pars[[a.x]], pars[[a.y]])
      pars[[paste0(a,'.fixed')]] <- coalesce(ifelse(is.na(pars[[a.x]]), NA, FALSE), ifelse(is.na(pars[[a.y]]), NA, TRUE))
    } else {
      pars[[paste0(a,'.fixed')]] <- a %in% metab.ddat
    }
  }

  # identify what we actually have, in the order we want it
  metab.out <- metab.all[metab.all %in% names(pars)]
  
  # add uncertainty columns if requested
  if(uncertainty != 'none') {
    metab.uncert <- paste0(rep(metab.out, each=1), rep(c('.sd'), times=length(metab.out)))
    metab.out <- c(rbind(metab.out, metab.uncert)) %>% { .[. %in% names(pars)]}
  }
  
  # add .fixed columns to the list of exported columns if requested
  if(fixed %in% c('columns','stars')) {
    for(a in metab.either) {
      add.after <- tail(grep(paste0('^', a), metab.out), 1)
      metab.out <- append(metab.out, paste0(a,'.fixed'), after=add.after)
    }
  }
  
  # select and order those columns of pars that match metab.needs, 
  # metab.optional, or other columns we've added. useful to order now because 
  # mm_sd_to_ci will swap columns in place
  params <- pars[c('date', metab.out)]
  
  # convert sds to CIs if requested
  if(uncertainty == 'ci') {
    params <- mm_sd_to_ci(params)
  }
  
  # convert .fixed columns to stars if requested (do this after mm_sd_to_ci b/c converts to character)
  if(fixed == 'stars') {
    params <- bind_cols(select(params, date), format.data.frame(select(params, -date)))
    for(a in metab.either) {
      params[[a]] <- paste0(params[[a]], ifelse(params[[paste0(a,'.fixed')]], '*', ' '))
      params[[paste0(a,'.fixed')]] <- NULL
    }
  }
  
  # attach warnings and errors if requested
  if(messages && exists('date', pars) && any(exists(c('warnings','errors'), pars))) {
    messages <- pars[c('date','warnings','errors') %>% { .[. %in% names(pars)] }]
    params <- left_join(params, messages, by='date', copy=TRUE)
  }
  
  # attach units if requested and available in mm_data
  if(attach.units) {
    param.units <- get_units(mm_data())[names(params)]
    params <- u(params, param.units)
  }
  
  # return
  params
}

#' Make dissolved oxygen predictions from a fitted metab_model.
#' 
#' Makes fine-scale predictions of dissolved oxygen using fitted coefficients, 
#' etc. from the metabolism model.
#' 
#' @inheritParams predict_DO
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
  if(isTRUE(use_saved) && !is.null(data) && "DO.mod" %in% names(data)) {
    return(data)
  }
  
  # get the metabolism estimates; filter as we did for data
  metab_ests <- get_params(metab_model, date_start=date_start, date_end=date_end, uncertainty='none', messages=FALSE)
  
  # re-process the input data with the metabolism estimates to predict DO
  preds <- mm_model_by_ply(
    mm_predict_DO_1ply, data=data, data_daily=metab_ests, # for mm_model_by_ply
    day_start=day_start, day_end=day_end, day_tests=c(), timestep_days=FALSE, # for mm_model_by_ply
    model_name=specs$model_name) %>% # for mm_predict_DO_1ply
    mm_filter_dates(date_start=date_start, date_end=date_end) # trim off the extra
  
  if(attach.units) preds <- u(preds, get_units(mm_data())[gsub('DO.mod','DO.obs',names(preds))])
  preds
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
predict_metab.metab_model <- function(metab_model, date_start=NA, date_end=NA, ..., attach.units=FALSE, use_saved=TRUE) {
  
  if(isTRUE(use_saved) && !is.null(metab_model@metab_daily)) {
    # if allowed and available, use previously stored values rather than
    # re-predicting them now
    preds <- metab_model@metab_daily
  
  } else {
    # otherwise predict them now
    
    # get the metabolism parameters; filter if requested
    metab_ests <- get_params(metab_model, date_start=date_start, date_end=date_end, uncertainty='sd', messages=FALSE)
    
    # pull args from the model
    specs <- get_specs(metab_model)
    day_start <- specs$day_start
    day_end <- specs$day_end
    
    # get the instantaneous data, including DO.mod; filter if requested
    data <- predict_DO(metab_model, date_start=date_start, date_end=date_end, attach.units=FALSE, use_saved=TRUE)
    
    # re-process the input data with the metabolism estimates to predict daily mean metabolism & reaeration (D)
    preds <- mm_model_by_ply(
      mm_predict_metab_1ply, data=data, data_daily=metab_ests, # for mm_model_by_ply
      day_start=day_start, day_end=day_end, day_tests=c(), timestep_days=FALSE, # for mm_model_by_ply
      model_name=specs$model_name) # for mm_predict_DO_1ply
    
    # attach warnings and errors if available
    fit <- get_fit(metab_model)
    if(!is.null(fit) && exists('date', fit) && any(exists(c('warnings','errors'), fit))) {
      messages <- fit[c('date','warnings','errors') %>% { .[. %in% names(fit)] }]
      preds <- full_join(preds, messages, by='date', copy=TRUE)
    }
  }
  
  if(attach.units) {
    pred.units <- get_units(mm_data())[sapply(names(preds), function(x) strsplit(x, '\\.')[[1]][1], USE.NAMES=FALSE)]
    preds <- u(preds, pred.units)
  }
  preds
}
