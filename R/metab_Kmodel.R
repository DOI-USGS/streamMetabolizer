#' @include metab_model-class.R
NULL

#' Combine a time series of K estimates to predict unified values
#' 
#' Takes daily estimates of K, usually from nighttime regression, and regresses 
#' against predictors such as discharge.daily. Returns a metab_model that only 
#' predicts daily K, nothing else.
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
#' @param method the name of an interpolation or regression method relating K to
#'   the predictor[s] of choice.
#' @inheritParams metab_model_prototype
#' @inheritParams prepdata_Kmodel
#' @inheritParams Kmodel_allply
#' @inheritParams mm_is_valid_day
#' @param ... Other arguments passed to the fitting function given by 
#'   \code{method}. \code{na.rm=TRUE} is already passed to \code{mean}
#' @return A metab_Kmodel object containing the fitted model.
#' @examples
#' \dontrun{
#' library(dplyr)
#' fr <- streamMetabolizer:::load_french_creek() %>% unitted::v() %>% 
#'   filter(format(local.time, "%Y-%m-%d") == "2012-08-24") %>% 
#'   select(local.time)
#' fr <- fr %>% mutate(discharge.daily=exp(rnorm(nrow(fr))))
#' frk2 <- data.frame(local.date=seq(as.Date("2012-08-15"),as.Date("2012-09-15"),
#'   as.difftime(1,units='days')), discharge.daily=exp(rnorm(32,2,1)), K600=rnorm(32,30,4))
#'   
#' # 1-day tests
#' frk <- data.frame(local.date=as.Date("2012-08-24"), K600=20, stringsAsFactors=FALSE)
#' metab_Kmodel(data=fr, data_daily=frk, method='mean')
#' #metab_Kmodel(data=fr, data_daily=frk, method='lm') # NaN warning
#' #metab_Kmodel(data=fr, data_daily=frk, method='loess', predictors='local.date') # span error
#' 
#' # mean
#' mm <- metab_Kmodel(data_daily=frk2, method='mean')
#' 
#' # lm
#' mm <- metab_Kmodel(data_daily=frk2, method='lm') # very similar to mean
#' mm <- metab_Kmodel(data_daily=frk2, method='lm', predictors=c('discharge.daily'))
#' mm <- metab_Kmodel(data_daily=frk2, method='lm', predictors=c('local.date'))
#' mm <- metab_Kmodel(data_daily=frk2, method='lm', predictors=c('discharge.daily','local.date'))
#' 
#' # loess
#' mm <- metab_Kmodel(data_daily=frk2, method='loess', predictors='local.date')
#' mm <- metab_Kmodel(data_daily=frk2, method='loess', predictors='local.date', span=0.4)
#' mm <- metab_Kmodel(data_daily=frk2, method='loess', predictors=c('discharge.daily','local.date'))
#' mm <- metab_Kmodel(data_daily=frk2, method='loess', predictors='discharge.daily', span=2)
#' 
#' library(ggplot2)
#' ggplot(get_data_daily(mm), aes(x=local.date)) + geom_line(aes(y=K600)) + 
#'   geom_point(aes(y=K600.obs)) + geom_ribbon(aes(ymin=K600.lower, ymax=K600.upper), alpha=0.4)
#' }
#' @export
#' @import dplyr
#' @importFrom magrittr %<>%
#' @family metab_model
metab_Kmodel <- function(
  data=mm_data(local.time, discharge, velocity, optional=c("all")), 
  data_daily=mm_data(local.date, K600, K600.lower, K600.upper, discharge.daily, velocity.daily, optional=c("K600.lower", "K600.upper", "discharge.daily", "velocity.daily")), # inheritParams metab_model_prototype
  method=c("mean", "lm", "loess"), weights=c("CI", "CI/K600"), predictors=c("local.date", "velocity.daily", "discharge.daily"), 
  filters=c(CI.max=NA, discharge.daily.max=NA, velocity.daily.max=NA), transforms=c(K600='log', local.date=NA, velocity.daily="log", discharge.daily="log"),
  info=NULL, day_start=4, day_end=27.99, # inheritParams metab_model_prototype. day_start and day_end only relevant if !is.null(data)
  tests=c('full_day', 'even_timesteps', 'complete_data'), # args for mm_is_valid_day. only relevant if !is.null(data)
  ...
) {
  
  method <- match.arg(method)
  weights <- if(missing(weights)) c() else match.arg(weights)
  predictors <- if(missing(predictors)) c() else match.arg(predictors, several.ok=TRUE)
  if('local.date' %in% predictors) predictors[predictors=='local.date'] <- 'as.numeric(local.date)'
  if(!(!('K600' %in% names(transforms)) || is.na(transforms[['K600']]) || isTRUE(transforms[['K600']]=='log'))) 
    stop("transforms['K600'] should be NA or 'log'")
  
  fitting_time <- system.time({
    # Check data for correct column names & units
    dat_list <- mm_validate_data(if(missing(data)) NULL else data, data_daily, "metab_Kmodel")
    data <- dat_list[['data']]
    data_daily <- dat_list[['data_daily']]
    
    # Prepare data_daily by aggregating any daily data, renaming K600 to 
    # K600.obs, & setting data_daily$weight to reflect user weights & filters
    data_daily <- prepdata_Kmodel(data=data, data_daily=data_daily, weights=weights, filters=filters, day_start=day_start, day_end=day_end, tests=tests)
    
    # Fit the model
    Kmodel_all <- Kmodel_allply(data_daily_all=data_daily, method=method, weights=weights, predictors=predictors, transforms=transforms, ...)
    
  })
  
  # Package the results for prediction
  mm <- metab_model(
    "metab_Kmodel", 
    info=info,
    fit=Kmodel_all,
    fitting_time=fitting_time,
    args=list(
      method=method, weights=weights, predictors=predictors, filters=filters, transforms=transforms,
      day_start=day_start, day_end=day_end, tests=tests),
    data=data,
    data_daily=data_daily)
  
  # Update data_daily with predictions
  preds <- predict_metab(mm)
  mm@data_daily %<>% left_join(select(preds, local.date, K600, K600.lower, K600.upper), by='local.date')
  
  # Return
  mm
}


#### helpers ####

#' Prepare data_daily by aggregating any daily data, renaming K600 to K600.obs, 
#' & setting data_daily$weight to reflect user weights & filters
#' 
#' @param data unit data to aggregate to daily_data. may be NULL.
#' @param data_daily daily data to prepare for K modeling
#' @param weights character vector indicating the type of weighting to use. 
#'   Leave blank or set to c() for no weights.
#' @param filters named numeric vector of limits to use in filtering data_daily.
#'   If an element with a name in c("CI.max","discharge.daily.max","velocity.daily.max") is 
#'   given, the corresponding filter is applied: K600.upper-K600.lower >= 
#'   CI.max, discharge.daily <= discharge.daily.max, velocity.daily <= velocity.daily.max
#' @inheritParams metab_model_prototype
#' @inheritParams mm_is_valid_day
prepdata_Kmodel <- function(data, data_daily, weights, filters, day_start, day_end, tests) {
  # Aggregate unit data to daily timesteps if needed
  if(!is.null(data) && nrow(data) > 0) {
    columns <- c('discharge', 'velocity')[c('discharge', 'velocity') %in% names(data)]
    if(length(columns) == 0) stop("data arg is pointless without at least one of c('discharge', 'velocity')")
    aggs_daily <- mm_model_by_ply(Kmodel_aggregate_day, data=data, data_daily=NULL, day_start=day_start, day_end=day_end, tests=tests, columns=columns)
    names(aggs_daily)[match(columns, names(aggs_daily))] <- paste0(columns, '.daily')
    data_daily <- left_join(data_daily, aggs_daily, by="local.date")
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
    data_daily %<>% mutate(weight=switch(
      weights, 
      "CI"=1/(K600.upper.obs - K600.lower.obs), 
      "CI/K600"=1/(K600.upper.obs - K600.lower.obs)/K600.obs)) %>%
      mutate(weight = weight/sum(weight, na.rm=TRUE))
  } else {
    data_daily %<>% mutate(weight=1/length(which(!is.na(K600))))
  }
  
  # Filter out undesired days. Indicate filtering by weights
  if(!isTRUE(is.na(filters[['CI.max']])))
    data_daily %<>% mutate(weight=weight * if(K600.upper.obs-K600.lower.obs <= filters[['CI.max']]) 1 else 0)
  if(!isTRUE(is.na(filters[['discharge.daily.max']])))
    data_daily %<>% mutate(weight=weight * if(discharge.daily <= filters[['discharge.daily.max']]) 1 else 0)
  if(!isTRUE(is.na(filters[['velocity.daily.max']])))
    data_daily %<>% mutate(weight=weight * if(velocity.daily <= filters[['velocity.daily.max']]) 1 else 0)
  
  data_daily
}

#' Aggregate unit values to daily values
#' 
#' For use in predicting K
#' 
#' @inheritParams mm_model_by_ply_prototype
#' @inheritParams mm_is_valid_day
#' @param columns character vector of names of columns to aggregate
Kmodel_aggregate_day <- function(
  data_ply, data_daily_ply, day_start, day_end, local_date, # inheritParams mm_model_by_ply_prototype
  tests=c('full_day', 'even_timesteps', 'complete_data'), # inheritParams mm_is_valid_day
  columns=c('discharge', 'velocity')
) {
  # Provide ability to skip a poorly-formatted day for calculating 
  # metabolism, without breaking the whole loop. Just collect 
  # problems/errors as a list of strings and proceed. Also collect warnings.
  validity <- mm_is_valid_day(data_ply, day_start=day_start, day_end=day_end, tests=tests)
  stop_strs <- if(isTRUE(validity)) character(0) else validity
  
  # Calculate means of requested columns
  if(length(stop_strs) == 0) {
    cmeans_1day <- as.data.frame(t(colMeans(data_ply[columns])))
  } else {
    cmeans_1day <- as.data.frame(t(rep(as.numeric(NA),length(columns)))) %>% setNames(columns)
  }
  
  # Return, reporting any results, warnings, and errors
  data.frame(cmeans_1day,
             errors.agg=paste0(stop_strs, collapse="; "),
             stringsAsFactors=FALSE)
}

#' Fit a K model
#' 
#' The model will predict daily K estimates from preliminary daily K estimates.
#' Called from metab_Kmodel().
#' 
#' @param data_daily_all data to use as input, with columns including K600.obs,
#'   weight, and any predictors
#' @param predictors character vector of variables (column names in data or
#'   data_daily) to use in predicting K. Leave blank or set to c() for no
#'   predictors.
#' @inheritParams metab_Kmodel
#' @keywords internal
Kmodel_allply <- function(data_daily_all, method, weights, predictors, transforms, ...) {
  # remove rows whose weights signal they should be filtered out
  weight <- '.dplyr.var'
  data_daily_all <- dplyr::filter(data_daily_all, weight > 0)
  # add transformations
  trans_preds <- sapply(c("K600.obs", predictors), function(pred) {
    if(!is.na(transforms[pred])) {
      paste0(transforms[pred], "(", pred, ")")
    } else {
      pred
    }
  })
  # fit & return the model
  switch(
    method,
    mean = {
      if(length(predictors) > 0) warning("predictors ignored for method='mean'")
      list(
        mean=sum(data_daily_all$K600.obs * data_daily_all$weight, na.rm=TRUE, ...),
        se=(if(length(weights) > 0) {
          warning("omitting sd for weighted mean")
          NA 
        } else { 
          sd(data_daily_all$K600.obs, na.rm=TRUE)/sqrt(nrow(data_daily_all)) # SE of the mean
        }))
    }, 
    lm = {
      if(length(predictors) == 0) predictors <- "1"
      formul <- formula(paste0(trans_preds["K600.obs"], " ~ ", paste0(trans_preds[-1], collapse=" + ")))
      wts <- (if(length(weights) > 0) data_daily_all$weight else NULL)
      lm(formul, data=data_daily_all, weights=wts, ...)
    }, 
    loess = {
      if(length(predictors) < 1) stop("need at least one predictor for method='loess'")
      formul <- formula(paste0(trans_preds["K600.obs"], " ~ ", paste0(trans_preds[-1], collapse=" + ")))
      if(length(weights) > 0) {
        loess(formul, data=data_daily_all, weights=data_daily_all$weight, ...)
      } else {
        loess(formul, data=data_daily_all, ...)
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
  if(!isTRUE(use_saved) || is.null(data_daily) || !("K600.mod" %in% names(data_daily))) {
    method <- get_args(metab_model)$method
    ktrans <- get_args(metab_model)$transforms['K600']
    fit <- get_fit(metab_model)
    switch(
      method,
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
  }
  metab_model@fit <- data_daily # temporary for converting lower/upper/sd to standard colnames
  NextMethod()
}

#' Override generic predict_DO for metab_Kmodel, which can't
#' 
#' metab_Kmodel predicts K at daily timesteps and usually knows nothing about 
#' GPP or ER. So it's not possible to predict DO from this model. Try passing 
#' the ouptut to metab_mle and THEN predicting DO.
#' @inheritParams predict_DO
#' @examples
#' \dontrun{
#' vfrench <- streamMetabolizer:::load_french_creek(attach.units=FALSE)
#' 
#' # fit a first-round MLE and extract the K estimates
#' mm1 <- metab_mle(data=vfrench, day_start=-1, day_end=23)
#' K600_mm1 <- predict_metab(mm1) %>% select(local.date, K600, K600.lower, K600.upper)
#' 
#' # smooth the K600s
#' mm2 <- metab_Kmodel(data_daily=K600_mm1, method='mean', day_start=-1, day_end=23)
#' K600_mm2 <- predict_metab(mm2) %>% select(local.date, K600)
#' 
#' # refit the MLE with fixed K
#' mm3 <- metab_mle(data=vfrench, data_daily=K600_mm2, day_start=-1, day_end=23)
#' predict_metab(mm3)
#' }
predict_DO.metab_Kmodel <- function(metab_model, date_start=NA, date_end=NA, ..., use_saved=TRUE) {
  stop("can only predict K, not DO, from metab_Kmodel")
}