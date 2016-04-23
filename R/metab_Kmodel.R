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
#'   as.difftime(1,units='days')), discharge.daily=exp(rnorm(32,2,1)), K600=rnorm(32,30,4)) %>%
#'   mutate(K600.lower=K600-5, K600.upper=K600+6)
#' 
#' # mean
#' mm <- metab_Kmodel(
#'   specs(mm_name('Kmodel', engine='mean')), 
#'   data_daily=example_Ks) # two warnings expected for engine='mean'
#' \dontrun{ plot_metab_preds(predict_metab(mm)) # flat b/c it's a mean }
#' 
#' # linear model
#' mm <- metab_Kmodel(
#'   specs(mm_name('Kmodel', engine='lm'), predictors='discharge.daily'),
#'   data_daily=example_Ks)
#' \dontrun{ plot_metab_preds(predict_metab(mm)) }
#' 
#' # loess
#' mm <- metab_Kmodel(
#'   specs(mm_name('Kmodel', engine='loess'), predictors='date', other_args=list(span=0.4)),
#'   data_daily=example_Ks)
#' \dontrun{ plot_metab_preds(predict_metab(mm)) }
#' @export
#' @family metab_model
metab_Kmodel <- function(
  specs=specs(mm_name('Kmodel')),
  data=mm_data(solar.time, discharge, velocity, optional=c("all")), 
  data_daily=mm_data(date, K600, K600.lower, K600.upper, discharge.daily, velocity.daily, optional=c("K600.lower", "K600.upper", "discharge.daily", "velocity.daily")),
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
    
    # Prepare data_daily by aggregating any daily data, renaming K600 to 
    # K600.obs, & setting data_daily$weight to reflect user weights & filters
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
  
  # Update data_daily with predictions
  preds <- predict_metab(mm)
  mm@data_daily %<>% left_join(select(preds, date, K600, K600.lower, K600.upper), by='date')
  
  # Return
  mm
}


#### helpers ####

#' Prepare data_daily by aggregating any daily data, renaming K600 to K600.obs, 
#' & setting data_daily$weight to reflect user weights & filters
#' 
#' @param data unit data to aggregate to daily_data. may be NULL.
#' @param data_daily daily data to prepare for K modeling
#' @param weights For Kmodel, character vector indicating the type of weighting 
#'   to use. Set to c() for no weights. One of c("1/CI", "K600/CI", c()).
#' @param filters For Kmodel, named numeric vector of limits to use in filtering
#'   data_daily. Elements may include
#'   c("CI.max","discharge.daily.max","velocity.daily.max"). If an element is
#'   given, the corresponding filter is applied: K600.upper-K600.lower <=
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
  K600.obs <- K600.lower.obs <- K600.upper.obs <- weight <- '.dplyr.var'
  
  # Rename the data so it's clear which K values are inputs and which are new model predictions
  if('K600.lower' %in% names(data_daily)) {
    data_daily %<>% rename(K600.obs=K600, K600.lower.obs=K600.lower, K600.upper.obs=K600.upper)
  } else {
    data_daily %<>% rename(K600.obs=K600)
  }
  
  # Set weights
  if(length(weights) > 0) {
    if(!all(c('K600.lower.obs','K600.upper.obs') %in% names(data_daily))) {
      stop("need 'K600.lower', and 'K600.upper' in data_daily to set weights by CI or CI/K600")
    }
    data_daily %<>% 
      mutate(weight=switch(
        weights, 
        "1/CI"=1/(K600.upper.obs - K600.lower.obs), 
        "K600/CI"=pmax(K600.obs, 0)/(K600.upper.obs - K600.lower.obs))) %>%
      mutate(weight = weight/sum(weight, na.rm=TRUE))
  } else {
    data_daily %<>% mutate(weight=1/length(which(!is.na(K600.obs))))
  }
  
  out <- list(unfiltered = data_daily)
  
  # Filter out undesired days. Indicate filtering by weights
  if(('CI.max' %in% names(filters)) && !isTRUE(is.na(filters[['CI.max']])))
    data_daily %<>% mutate(weight=weight * if(K600.upper.obs-K600.lower.obs <= filters[['CI.max']]) 1 else 0)
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
#' @param data_daily_all data to use as input, with columns including K600.obs, 
#'   weight, and any predictors
#' @param predictors For Kmodel, character vector of variables (column names in 
#'   data or data_daily) to use in predicting K. Leave blank or set to c() for 
#'   no predictors. Otherwise, one or more of these may be included: c("date", 
#'   "velocity.daily", "discharge.daily").
#' @param transforms For Kmodel, a named character vector of names of functions 
#'   (probably 'log' or NA) to apply to K600 and/or the predictors. K600 should 
#'   probably be logged. The vector names must match the values of 
#'   \code{predictors}, although not all elements of \code{predictors} must be 
#'   included in \code{transforms}. Recommended transforms include 
#'   \code{c(K600='log', date=NA, velocity.daily="log", discharge.daily="log")}
#' @param other_args Other arguments passed to the fitting function given by 
#'   \code{specs$engine}. \code{na.rm=TRUE} is already passed to 
#'   \code{mean} (which is actually implemented as \code{sum}, anyway).
#' @inheritParams metab_Kmodel
#' @inheritParams specs
#' @keywords internal
Kmodel_allply <- function(data_daily_all, engine, weights, predictors, transforms, other_args) {
  # remove rows whose weights signal they should be filtered out
  weight <- '.dplyr.var'
  data_daily_all <- dplyr::filter(data_daily_all, weight > 0)
  # add transformations
  names(transforms) <- replace(names(transforms), names(transforms)=="K600", "K600.obs")
  trans_preds <- sapply(c("K600.obs", predictors), function(pred) {
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
      if(length(predictors) > 0) warning("predictors ignored for engine='mean'")
      Kobs <- eval(parse(text=trans_preds[['K600.obs']]), envir=data_daily_all)
      list(
        mean=sum(Kobs * data_daily_all$weight, na.rm=TRUE),
        se=(if(length(weights) > 0) {
          warning("omitting sd for weighted mean")
          NA 
        } else { 
          sd(Kobs, na.rm=TRUE)/sqrt(nrow(data_daily_all)) # SE of the mean
        }))
    }, 
    'lm' = {
      if(length(predictors) == 0) trans_preds[2] <- "1"
      formul <- formula(paste0(trans_preds["K600.obs"], " ~ ", paste0(trans_preds[-1], collapse=" + ")))
      wts <- if(length(weights) > 0) data_daily_all$weight else NULL
      do.call(lm, c(list(formula=formul, data=data_daily_all, weights=wts), other_args))
    }, 
    'loess' = {
      if(length(predictors) < 1) stop("need at least one predictor for engine='loess'")
      formul <- formula(paste0(trans_preds["K600.obs"], " ~ ", paste0(trans_preds[-1], collapse=" + ")))
      if(length(weights) > 0) {
        do.call(loess, c(list(formula=formul, data=data_daily_all, weights=data_daily_all$weight), other_args))
      } else {
        do.call(loess, c(list(formula=formul, data=data_daily_all), other_args))
      }
    })
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

#' Make metabolism predictions from a fitted metab_model.
#' 
#' Makes daily re-predictions of K600.
#' 
#' @inheritParams predict_metab
#' @return A data.frame of predictions, as for the generic 
#'   \code{\link{predict_metab}}.
#' @export
#' @importFrom magrittr %<>%
#' @family predict_metab
predict_metab.metab_Kmodel <- function(metab_model, date_start=NA, date_end=NA, ..., use_saved=TRUE) {
  
  # re-predict K600.mod if saved values are disallowed or unavailable; otherwise
  # use previously stored values for K600.mod
  data_daily <- get_data_daily(metab_model) %>%
    mm_filter_dates(date_start=date_start, date_end=date_end)
  if(!isTRUE(use_saved) || is.null(data_daily) || !("K600" %in% names(data_daily))) {
    engine <- get_specs(metab_model)$engine
    ktrans <- get_specs(metab_model)$transforms['K600']
    fit <- get_fit(metab_model)
    . <- '.dplyr.var'
    switch(
      engine,
      mean = {
        data_daily %<>% mutate(
          GPP = NA,
          GPP.sd = NA,
          ER = NA,
          ER.sd = NA,
          K600 = (fit$mean %>% if(!is.na(ktrans) && ktrans=='log') exp(.) else .), # the parens are needed to avoid "Error: invalid subscript type 'closure'"
          K600.sd = (fit$se %>% if(!is.na(ktrans) && ktrans=='log') {warning('no SE available for mean(log(K600))'); NA} else .))
      },
      lm = {
        preds <- predict(fit, newdata=data_daily, interval='confidence', level=0.95)
        data_daily %<>% mutate(
          GPP_daily_50pct = NA,
          GPP_daily_2.5pct = NA,
          GPP_daily_97.5pct = NA,
          ER_daily_50pct = NA,
          ER_daily_2.5pct = NA,
          ER_daily_97.5pct = NA,
          K600_daily_50pct = (preds[,'fit'] %>% if(!is.na(ktrans) && ktrans=='log') exp(.) else .),
          K600_daily_2.5pct = (preds[,'lwr'] %>% if(!is.na(ktrans) && ktrans=='log') exp(.) else .),
          K600_daily_97.5pct = (preds[,'upr'] %>% if(!is.na(ktrans) && ktrans=='log') exp(.) else .))
      },
      loess = {
        preds <- predict(fit, newdata=data_daily, se=TRUE)
        data_daily %<>% mutate(
          GPP = NA,
          GPP.sd = NA,
          ER = NA,
          ER.sd = NA,
          K600 = (preds$fit %>% if(!is.na(ktrans) && ktrans=='log') exp(.) else .),
          K600.sd = (preds$se.fit %>% if(!is.na(ktrans) && ktrans=='log') exp(.) else .))
      }
    )
  } else { # we're not storing GPP or ER, so add them back in for prediction
    data_daily %<>% 
      rename(
        K600_daily_50pct = K600,
        K600_daily_2.5pct = K600.lower,
        K600_daily_97.5pct = K600.upper) %>% 
      mutate(
        GPP_daily_50pct = NA,
        GPP_daily_2.5pct = NA,
        GPP_daily_97.5pct = NA,
        ER_daily_50pct = NA,
        ER_daily_2.5pct = NA,
        ER_daily_97.5pct = NA)
  }
  metab_model@fit <- mutate(data_daily, warnings='', errors='') # temporary for converting lower/upper/sd to standard colnames
  NextMethod()
}

#' Override generic predict_DO for metab_Kmodel, which can't predict DO
#' 
#' metab_Kmodel predicts K at daily timesteps and usually knows nothing about 
#' GPP or ER. So it's not possible to predict DO from this model. Try passing 
#' the ouptut to metab_mle and THEN predicting DO.
#' @inheritParams predict_DO
#' @examples
#' \dontrun{
#' vfrench <- streamMetabolizer:::load_french_creek(attach.units=FALSE)
#' 
#' # set the specifications for the MLE models
#' mle_specs <- specs(mm_name('mle'), day_start=-1, day_end=23)
#' 
#' # fit a first-round MLE and extract the K estimates
#' mm1 <- metab_mle(mle_specs, data=vfrench)
#' K600_mm1 <- predict_metab(mm1) %>% select(date, K600, K600.lower, K600.upper)
#' 
#' # smooth the K600s
#' mm2 <- metab_Kmodel(specs(mm_name('Kmodel', engine='mean'), 
#'   day_start=-1, day_end=23), data_daily=K600_mm1)
#' K600_mm2 <- predict_metab(mm2) %>% select(date, K600)
#' 
#' # refit the MLE with fixed K
#' mm3 <- metab_mle(mle_specs, data=vfrench, data_daily=K600_mm2)
#' predict_metab(mm3)
#' }
#' @export
predict_DO.metab_Kmodel <- function(metab_model, date_start=NA, date_end=NA, ..., use_saved=TRUE) {
  stop("can only predict K, not DO, from metab_Kmodel")
}