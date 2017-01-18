#' @include metab_model-class.R
NULL

#' Combine a time series of K estimates to predict consistent values
#' 
#' Takes daily estimates of K, usually from nighttime regression, and regresses 
#' against predictors such as discharge.daily. Returns a metab_Kmodel object
#' that only predicts daily K, nothing else.
#' 
#' Possible approaches:
#' 
#' \itemize{
#' 
#' \item{"mean"}{Predict K as the mean of all K values}
#' 
#' \item{"weighted mean"}{Predict K as the mean of all K values, weighted by the
#' inverse of the confidence intervals in the input K values}
#' 
#' \item{"KvQ"}{Regress K versus Q, tending toward overall mean in ranges of Q 
#' with sparse data}
#' 
#' \item{"weighted KvQ"}{Regress K versus Q, tending toward overall mean in 
#' ranges of Q with sparse data, weighting high-confidence K values more 
#' heavily}
#' 
#' \item{"T smoother"}{Predict K using a loess or spline smoother over time}
#' 
#' \item{"Q smoother"}{Predict K using a loess or spline smoother over 
#' discharge.daily}
#' 
#' \item{"TQ smoother"}{Predict K using a loess or spline smoother over both 
#' time and discharge.daily}
#' 
#' }
#' 
#' @author Alison Appling
#' 
#' @inheritParams metab
#' @return A metab_Kmodel object containing the fitted model. This object can be 
#'   inspected with the functions in the \code{\link{metab_model_interface}}.
#' @import dplyr
#' @importFrom magrittr %<>%
#' 
#' @examples
#' library(dplyr)
#' # create example data
#' set.seed(24842)
#' example_Ks <- data.frame(date=seq(as.Date("2012-08-15"),as.Date("2012-09-15"),
#'   as.difftime(1,units='days')), discharge.daily=exp(rnorm(32,2,1)), K600.daily=rnorm(32,30,4)) %>%
#'   mutate(K600.daily.lower=K600.daily-5, K600.daily.upper=K600.daily+6)
#' 
#' # mean
#' mm <- metab_Kmodel(
#'   specs(mm_name('Kmodel', engine='mean')), 
#'   data_daily=example_Ks) # two warnings expected for engine='mean'
#' get_params(mm)
#' \dontrun{
#' plot(get_params(mm)$date, get_params(mm)$K600.daily)
#' }
#' 
#' # linear model
#' mm <- metab_Kmodel(
#'   specs(mm_name('Kmodel', engine='lm'), predictors='discharge.daily'),
#'   data_daily=example_Ks)
#' get_params(mm)
#' \dontrun{
#' plot(get_data_daily(mm)$discharge.daily, get_params(mm)$K600.daily)
#' }
#' 
#' # loess
#' mm <- metab_Kmodel(    ### breaks ###
#'   specs(mm_name('Kmodel', engine='loess'), predictors='date', other_args=list(span=0.4)),
#'   data_daily=example_Ks)
#' get_params(mm)
#' \dontrun{
#' plot(get_params(mm)$date, get_params(mm)$K600.daily)
#' }
#' 
#' ## 3-phase workflow (sort of like complete pooling) for estimating K within 
#' ## days, then K across days, then GPP and ER within days
#' 
#' # 1. data and specifications for both of the MLE models
#' dat <- data_metab('10','15')
#' mle_specs <- specs(mm_name('mle'))
#' 
#' # fit a first-round MLE and extract the K estimates
#' mm1 <- metab_mle(mle_specs, data=dat)
#' K600_mm1 <- get_params(mm1, uncertainty='ci') %>% 
#'   select(date, K600.daily, K600.daily.lower, K600.daily.upper)
#' 
#' # smooth the K600s
#' mm2 <- metab_Kmodel(specs(mm_name('Kmodel', engine='mean'), 
#'   day_start=-1, day_end=23), data_daily=K600_mm1)
#' K600_mm2 <- get_params(mm2) %>% select(date, K600.daily)
#' 
#' # refit the MLE with fixed K
#' mm3 <- metab_mle(mle_specs, data=dat, data_daily=K600_mm2)
#' get_params(mm3, fixed='stars')
#' predict_metab(mm3)
#' \dontrun{
#' plot_metab_preds(mm3)
#' }
#' @export
#' @family metab_model
metab_Kmodel <- function(
  specs=specs(mm_name('Kmodel')),
  data=mm_data(solar.time, discharge, velocity, optional=c("all")), 
  data_daily=mm_data(date, K600.daily, K600.daily.lower, K600.daily.upper, discharge.daily, velocity.daily,
                     optional=c("K600.daily.lower", "K600.daily.upper", "discharge.daily", "velocity.daily")),
  info=NULL
) {
  
  if(missing(specs)) {
    # if specs is left to the default, it gets confused about whether specs() is
    # the argument or the function. tell it which:
    specs <- streamMetabolizer::specs(mm_name('Kmodel'))
  }
  # Fit the model
  fitting_time <- system.time({
    # Check and reformat arguments
    specs$engine <- match.arg(specs$engine, choices=c("mean", "lm", "loess"))
    if(length(specs$weights) > 0) specs$weights <-  match.arg(specs$weights, choices=c("1/CI","K600/CI"))
    if(length(specs$predictors) > 0) specs$predictors <- match.arg(specs$predictors, choices=c("date", "velocity.daily", "discharge.daily"), several.ok=TRUE)
    if('date' %in% specs$predictors) 
      specs$predictors[specs$predictors=='date'] <- 'as.numeric(date)'
    if(!(!('K600' %in% names(specs$transforms)) || 
         is.na(specs$transforms[['K600']]) || 
         isTRUE(specs$transforms[['K600']]=='log'))) 
      stop("specs$transforms['K600'] should be NA or 'log'")
    
    # Check data for correct column names & units
    dat_list <- mm_validate_data(if(missing(data)) NULL else data, if(missing(data_daily)) NULL else data_daily, "metab_Kmodel")
    data <- dat_list[['data']]
    data_daily <- dat_list[['data_daily']]
    
    # Prepare data_daily by aggregating any daily data, renaming K600.daily to 
    # K600.daily.obs, & setting data_daily$weight to reflect user weights & filters
    data_list <- prepdata_Kmodel(
      data=data, data_daily=data_daily, weights=specs$weights, filters=specs$filters, 
      day_start=specs$day_start, day_end=specs$day_end, day_tests=specs$day_tests)
    
    # Fit the model
    Kmodel_all <- do.call(Kmodel_allply, c(
      list(data_daily_all=data_list$filtered), 
      specs[c('engine','weights','predictors','transforms','other_args')]))
    
  })
  
  # Package the results for prediction
  mm <- metab_model(
    "metab_Kmodel", 
    info=info,
    fit=Kmodel_all,
    fitting_time=fitting_time,
    specs=specs,
    data=data,
    data_daily=data_list$unfiltered)
  
  # Add daily predictions of K
  mm@metab_daily <- get_params(mm)
  
  # Return
  mm
}


#### helpers ####

#' Prepare data_daily by aggregating any daily data, renaming K600.daily to K600.daily.obs, 
#' & setting data_daily$weight to reflect user weights & filters
#' 
#' @param data unit data to aggregate to daily_data. may be NULL.
#' @param data_daily daily data to prepare for K modeling
#' @param weights For Kmodel, character vector indicating the type of weighting 
#'   to use. Set to c() for no weights. One of c("1/CI", "K600/CI", c()).
#' @param filters For Kmodel, named numeric vector of limits to use in filtering
#'   data_daily. Elements may include
#'   c("CI.max","discharge.daily.max","velocity.daily.max"). If an element is
#'   given, the corresponding filter is applied: K600.daily.upper-K600.daily.lower <=
#'   CI.max, discharge.daily <= discharge.daily.max, velocity.daily <= 
#'   velocity.daily.max
#' @inheritParams metab
#' @inheritParams mm_model_by_ply
prepdata_Kmodel <- function(data, data_daily, weights, filters, day_start, day_end, day_tests) {
  # Aggregate unit data to daily timesteps if needed
  if(!is.null(data) && nrow(data) > 0) {
    columns <- c('discharge', 'velocity')[c('discharge', 'velocity') %in% names(data)]
    if(length(columns) == 0) stop("data arg is pointless without at least one of c('discharge', 'velocity')")
    aggs_daily <- mm_model_by_ply(
      Kmodel_aggregate_day, 
      data=data, data_daily=NULL, day_start=day_start, day_end=day_end, day_tests=day_tests, 
      columns=columns)
    names(aggs_daily)[match(columns, names(aggs_daily))] <- paste0(columns, '.daily')
    data_daily <- left_join(data_daily, aggs_daily, by="date")
  }
  
  # Avoid global variable warnings in R CMD check
  K600.daily.obs <- K600.daily.lower.obs <- K600.daily.upper.obs <- weight <- '.dplyr.var'
  
  # Rename the data so it's clear which K values are inputs and which are new model predictions
  if('K600.daily.lower' %in% names(data_daily)) {
    data_daily %<>% rename(K600.daily.obs=K600.daily, K600.daily.lower.obs=K600.daily.lower, K600.daily.upper.obs=K600.daily.upper)
  } else {
    data_daily %<>% rename(K600.daily.obs=K600.daily)
  }
  
  # Set weights
  if(length(weights) > 0) {
    if(!all(c('K600.daily.lower.obs','K600.daily.upper.obs') %in% names(data_daily))) {
      stop("need 'K600.daily.lower', and 'K600.daily.upper' in data_daily to set weights by CI or CI/K600")
    }
    data_daily %<>% 
      mutate(weight=switch(
        weights, 
        "1/CI"=1/(K600.daily.upper.obs - K600.daily.lower.obs), 
        "K600/CI"=pmax(K600.daily.obs, 0)/(K600.daily.upper.obs - K600.daily.lower.obs))) %>%
      mutate(weight = weight/sum(weight, na.rm=TRUE))
  } else {
    data_daily %<>% mutate(weight=1/length(which(!is.na(K600.daily.obs))))
  }
  
  out <- list(unfiltered = data_daily)
  
  # Filter out undesired days. Indicate filtering by weights
  if(('CI.max' %in% names(filters)) && !isTRUE(is.na(filters[['CI.max']])))
    data_daily %<>% mutate(weight=weight * if(K600.daily.upper.obs-K600.daily.lower.obs <= filters[['CI.max']]) 1 else 0)
  if(('discharge.daily.max' %in% names(filters)) && !isTRUE(is.na(filters[['discharge.daily.max']])))
    data_daily %<>% mutate(weight=weight * if(discharge.daily <= filters[['discharge.daily.max']]) 1 else 0)
  if(('velocity.daily.max' %in% names(filters)) && !isTRUE(is.na(filters[['velocity.daily.max']])))
    data_daily %<>% mutate(weight=weight * if(velocity.daily <= filters[['velocity.daily.max']]) 1 else 0)
  
  out <- c(out, list(filtered = data_daily))
  out
}

#' Aggregate unit values to daily values
#' 
#' For use in predicting K
#' 
#' @inheritParams mm_model_by_ply_prototype
#' @param columns character vector of names of columns to aggregate
#' @keywords internal
Kmodel_aggregate_day <- function(
  data_ply, ply_validity, ..., # inheritParams mm_model_by_ply_prototype
  columns=c('discharge', 'velocity') # new arg for this function
) {
  # Collect errors & warnings as character vectors
  stop_strs <- if(isTRUE(ply_validity)) character(0) else ply_validity
  
  # Calculate means of requested columns
  if(length(stop_strs) == 0) {
    cmeans_1day <- as.data.frame(t(colMeans(data_ply[columns])))
  } else {
    cmeans_1day <- as.data.frame(t(rep(as.numeric(NA),length(columns)))) %>% setNames(columns)
  }
  
  # Return, reporting any results, warnings, and errors
  data.frame(cmeans_1day,
             errors.agg=paste0(unique(stop_strs), collapse="; "),
             stringsAsFactors=FALSE)
}

#' Fit a K model
#' 
#' The model will predict daily K estimates from preliminary daily K estimates. 
#' Called from metab_Kmodel().
#' 
#' @param data_daily_all data to use as input, with columns including K600.daily.obs, 
#'   weight, and any predictors
#' @param predictors For Kmodel, character vector of variables (column names in 
#'   data or data_daily) to use in predicting K. Leave blank or set to c() for 
#'   no predictors. Otherwise, one or more of these may be included: c("date", 
#'   "velocity.daily", "discharge.daily").
#' @param transforms For Kmodel, a named character vector of names of functions 
#'   (probably 'log' or NA) to apply to K600.daily and/or the predictors. K600.daily should 
#'   probably be logged. The vector names must match the values of 
#'   \code{predictors}, although not all elements of \code{predictors} must be 
#'   included in \code{transforms}. Recommended transforms include 
#'   \code{c(K600.daily='log', date=NA, velocity.daily="log", discharge.daily="log")}
#' @param other_args Other arguments passed to the fitting function given by 
#'   \code{specs$engine}. \code{na.rm=TRUE} is already passed to 
#'   \code{mean} (which is actually implemented as \code{sum}, anyway).
#' @inheritParams metab_Kmodel
#' @inheritParams specs
#' @importFrom stats sd formula loess
#' @keywords internal
Kmodel_allply <- function(data_daily_all, engine, weights, predictors, transforms, other_args) {
  
  # Provide ability to skip a poorly-formatted dataset for calculating 
  # metabolism. Collect problems/errors as a list of strings and proceed. Also
  # collect warnings.
  stop_strs <- warn_strs <- character(0)
  
  # remove rows whose weights signal they should be filtered out
  weight <- '.dplyr.var'
  data_daily_all <- dplyr::filter(data_daily_all, weight > 0)
  # add transformations
  names(transforms) <- replace(names(transforms), names(transforms)=="K600", "K600.daily.obs")
  trans_preds <- sapply(c("K600.daily.obs", predictors), function(pred) {
    if(!is.na(transforms[pred])) {
      paste0(transforms[pred], "(", pred, ")")
    } else {
      pred
    }
  })
  # fit & return the model
  other_args[['possible_args']] <- NULL # might be left over from specs(), but don't want it now
  
  switch(
    engine,
    'mean' = {
      if(length(predictors) > 0) {
        warn_strs <- c(warn_strs, "predictors ignored for engine='mean'")
      }
      Kobs <- eval(parse(text=trans_preds[['K600.daily.obs']]), envir=data_daily_all)
      overall <- list(
        mean=sum(Kobs * data_daily_all$weight, na.rm=TRUE),
        se=(if(length(weights) > 0) {
          warn_strs <- c(warn_strs, "omitting sd for weighted mean")
          NA 
        } else { 
          sd(Kobs, na.rm=TRUE)/sqrt(nrow(data_daily_all)) # SE of the mean
        }))
    }, 
    'lm' = {
      if(length(predictors) == 0) trans_preds[2] <- "1"
      formul <- formula(paste0(trans_preds["K600.daily.obs"], " ~ ", paste0(trans_preds[-1], collapse=" + ")))
      wts <- if(length(weights) > 0) data_daily_all$weight else NULL
      overall <- do.call(lm, c(list(formula=formul, data=data_daily_all, weights=wts), other_args))
    }, 
    'loess' = {
      if(length(predictors) < 1) {
        stop("need at least one predictor for engine='loess'") # stop rather than stop_strs because it's a poorly formatted request
      }
      formul <- formula(paste0(trans_preds["K600.daily.obs"], " ~ ", paste0(trans_preds[-1], collapse=" + ")))
      if(length(weights) > 0) {
        overall <- do.call(loess, c(list(formula=formul, data=data_daily_all, weights=data_daily_all$weight), other_args))
      } else {
        overall <- do.call(loess, c(list(formula=formul, data=data_daily_all), other_args))
      }
    })
  
  list(overall=overall, warnings=warn_strs, errors=stop_strs)
}


#### helpers to the helper ####


#### metab_Kmodel class ####

#' Interpolation model of daily K for metabolism
#' 
#' \code{metab_Kmodel} models use initial daily estimates of K, along with 
#' predictors such as Q (discharge.daily) or U (velocity.daily) or T (time) to leverage all
#' available data to reach better, less variable daily estimates of K
#' 
#' @exportClass metab_Kmodel
#' @family metab.model.classes
setClass(
  "metab_Kmodel", 
  contains="metab_model"
)

#' @describeIn get_params Make daily re-predictions of K600.daily based on the 
#'   across-days model of K600.daily versus predictors. Only returns estimates
#'   for K600.daily, not any of the other daily parameters
#' @export
#' 
#' @inheritParams predict_metab
#' @import dplyr
#' @importFrom magrittr %<>%
#' @importFrom stats predict
get_params.metab_Kmodel <- function(
  metab_model, date_start=NA, date_end=NA, 
  uncertainty=c('sd','ci','none'), messages=TRUE, fixed=c('none','columns','stars'), 
  ..., attach.units=FALSE, use_saved=TRUE) {
  
  # re-predict K600.daily.mod if saved values are disallowed or unavailable; otherwise
  # use previously stored values for K600.daily.mod
  if(isTRUE(use_saved) && is.list(metab_model@fit) && exists('daily', metab_model@fit)) {
    params <- metab_model@fit$daily %>% 
      mm_filter_dates(date_start=date_start, date_end=date_end)
  } else {
    data_daily <- get_data_daily(metab_model) %>%
      mm_filter_dates(date_start=date_start, date_end=date_end)
    ktrans <- get_specs(metab_model)$transforms[['K600']]
    do_ktrans <- !is.na(ktrans) && ktrans=='log'
    engine <- get_specs(metab_model)$engine
    fit <- get_fit(metab_model)$overall
    . <- '.dplyr.var'
    params <- select(data_daily, date)
    switch(
      engine,
      mean = {
        params %<>% mutate(
          K600.daily = fit[['mean']],
          K600.daily.sd = fit[['se']])
      },
      lm = {
        preds <- predict(fit, newdata=data_daily, interval='confidence', level=0.95)
        params %<>% mutate(
          K600.daily = preds[,'fit'], # only the approx mean if ktrans=='log'
          K600.daily.sd = NA, # this shouldn't be necessary after resolving #238
          K600_daily_50pct = preds[,'fit'],
          K600_daily_2.5pct = preds[,'lwr'],
          K600_daily_97.5pct = preds[,'upr'])
      },
      loess = {
        preds <- predict(fit, newdata=data_daily, se=TRUE)
        params %<>% mutate(
          K600.daily = {preds$fit}, # {} avoids "Error: invalid subscript type 'list'"
          K600.daily.sd = {preds$se.fit})
      })
    if(do_ktrans) {
      for(col in setdiff(names(params), 'date')) {
        params[[col]] <- exp(params[[col]])
      }
    }
  }
  metab_model@fit <- c(list(daily=mutate(params, warnings='', errors='')), get_fit(metab_model)) # temporarily assign to @fit for compatibility with get_params.metab_model
  # code duplicated in get_params.metab_bayes:
  if(length(metab_model@fit$warnings) > 0) {
    omsg <- 'overall warnings'
    dmsg <- metab_model@fit$daily$warnings
    metab_model@fit$daily$warnings <- ifelse(dmsg == '', omsg, paste(omsg, dmsg, sep=';'))
  }
  if(length(metab_model@fit$errors) > 0) {
    omsg <- 'overall errors'
    dmsg <- metab_model@fit$daily$errors
    metab_model@fit$daily$errors <- ifelse(dmsg == '', omsg, paste(omsg, dmsg, sep=';'))
  }
  metab_model@fit <- metab_model@fit$daily # SUPER-TEMPORARY we're still converting fit$daily to fit until #247, #229
  NextMethod()
}

#' Override generic predict_metab for metab_Kmodel, which can't predict metab
#' 
#' metab_Kmodel predicts K (only) at daily timesteps and usually knows nothing
#' about GPP or ER. So it's not possible to predict metabolism from this model.
#' Try get_params() to retrieve the predicted values of K600.daily.
#' @inheritParams predict_metab
#' @export
predict_metab.metab_Kmodel <- function(
  metab_model, date_start=NA, date_end=NA, day_start=NA, day_end=NA, ..., attach.units=FALSE, use_saved=TRUE) {
  stop("can only predict K600.daily, not metabolism, from metab_Kmodel. try get_params() instead")
}

#' @describeIn predict_DO Throws an error because models of type 'Kmodel' can't 
#'   predict DO. \code{metab_Kmodel} predicts K at daily timesteps and usually
#'   knows nothing about GPP or ER. So it's not possible to predict DO from this
#'   model. Try passing the output to metab_mle and THEN predicting DO.
#' @export
predict_DO.metab_Kmodel <- function(metab_model, date_start=NA, date_end=NA, ..., use_saved=TRUE) {
  stop("can only predict K, not DO, from metab_Kmodel. try get_params() instead")
}
