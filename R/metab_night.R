#' @include metab_model-class.R
NULL

#' Nighttime regression for K estimation
#' 
#' Fits a model to estimate K from nighttime input data on DO, temperature, 
#' light, etc. The default day start & end are 12 noon on the preceding to 
#' present day; the algorithm then filters the data to just those time points 
#' for which light is very low.
#' 
#' @author Alison Appling, Maite Arroita, Bob Hall
#'   
#' @inheritParams metab
#' @return A metab_night object containing the fitted model. This object can be 
#'   inspected with the functions in the \code{\link{metab_model_interface}}.
#'   
#' @examples
#' dat <- data_metab('3', day_start=12, day_end=35)
#' mm <- metab_night(data=dat)
#' predict_metab(mm)
#' \dontrun{
#' plot_DO_preds(predict_DO(mm))
#' }
#' @export
#' @family metab_model
metab_night <- function(
  specs=specs(mm_name('night')),
  data=mm_data(solar.time, DO.obs, DO.sat, depth, temp.water, light), 
  data_daily=mm_data(NULL),
  info=NULL
) {
  
  if(missing(specs)) {
    # if specs is left to the default, it gets confused about whether specs() is
    # the argument or the function. tell it which:
    specs <- streamMetabolizer::specs(mm_name('night'))
  }
  fitting_time <- system.time({
    # Check data for correct column names & units
    dat_list <- mm_validate_data(if(missing(data)) NULL else data, if(missing(data_daily)) NULL else data_daily, "metab_night")
    data <- v(dat_list[['data']])
    data_daily <- v(dat_list[['data_daily']])
    
    # model the data, splitting into potentially overlapping 'plys' for each date
    night_all <- mm_model_by_ply(
      nightreg_1ply, 
      data=data, data_daily=data_daily, day_start=specs$day_start, day_end=specs$day_end, # for mm_model_by_ply
      day_tests=c(), timestep_days=FALSE, # for mm_model_by_ply - bypass the regular timestep calc & validity check on the full day ply
      night_tests=specs$day_tests, specs=specs) # for nightreg_1ply
  })
  
  # Package and return results
  mm <- metab_model(
    model_class="metab_night", 
    info=info,
    fit=night_all,
    fitting_time=fitting_time,
    specs=specs,
    data=dat_list[['data']], # keep the units if given
    data_daily=dat_list[['data_daily']])
  
  # Update data with DO predictions
  mm@data <- predict_DO(mm)
  
  # Return
  mm
}


#### helpers ####

#' Make daily reaeration estimates from input parameters
#' 
#' Called from metab_night().
#' 
#' @inheritParams mm_model_by_ply_prototype
#' @param night_tests character vector of validity tests to conduct on the data 
#'   after subsetting to just nighttime
#' @inheritParams metab
#' @return data.frame of estimates and \code{\link[stats]{lm}} model diagnostics
#' @keywords internal
#' @references Hornberger, George M., and Mahlon G. Kelly. Atmospheric 
#'   Reaeration in a River Using Productivity Analysis. Journal of the 
#'   Environmental Engineering Division 101, no. 5 (October 1975): 729-39.
#'   
#'   Raymond, Peter A., Christopher J. Zappa, David Butman, Thomas L. Bott, Jody
#'   Potter, Patrick Mulholland, Andrew E. Laursen, William H. McDowell, and 
#'   Denis Newbold. Scaling the gas transfer velocity and hydraulic geometry in 
#'   streams and small rivers. Limnology & Oceanography: Fluids & Environments 2
#'   (2012); 41:53.
#' @importFrom utils head tail
#' @importFrom stats lm coef setNames
nightreg_1ply <- function(
  data_ply, data_daily_ply, day_start, day_end, ply_date, ..., # inheritParams mm_model_by_ply_prototype
  night_tests=TRUE, # new arg here
  specs=specs('n_np_pi_eu_.lm') # inheritParams metab
) {
  
  # Try to run the model. Collect warnings/errors as a list of strings and
  # proceed to the next data_ply if there are stop-worthy issues.
  stop_strs <- warn_strs <- character(0)
  night.1d <- withCallingHandlers(
    tryCatch({
      # subset to times of darkness. require that these are all consecutive -
      # it'd be meaningless to look at the diff from 1 night to the next 2 nights
      which_night <- which(v(data_ply$light) < v(u(0.1, "umol m^-2 s^-1")))
      if(length(which_night) == 0) stop("no nighttime rows in data_ply")
      if(any(diff(which_night) > 1)) 
        stop_strs <- c(stop_strs, "need exactly one night per data_ply")
      night_dat <- data_ply[which_night,]
      
      # check for data validity here, after subsetting to the time window we 
      # really care about
      which_twilight <- c(head(which_night,1) - 1, tail(which_night,1) + 1)
      has_sunset <- which_twilight[1] >= 1
      has_sunrise <- which_twilight[2] <= nrow(data_ply)
      if('include_sunset' %in% night_tests) {
        if(!has_sunset)
          stop_strs <- c(stop_strs, "data don't include day-night transition")
        night_tests <- night_tests[night_tests != 'include_sunset']
      }
      # for the full_day test, let the twilight hours define a narrower required
      # window than day_start, day_end did (if the twilight hours can be
      # determined from the data_ply)
      date_start <- as.POSIXct(paste0(as.character(ply_date), " 00:00:00"), tz=lubridate::tz(v(data_ply$solar.time)))
      night_hours <- as.numeric(data_ply[which_twilight+c(1,-1),'solar.time'] - date_start, units='hours')
      night_start <- if(has_sunset) max(day_start, night_hours[1]) else day_start
      night_end <- if(has_sunrise) min(day_end, night_hours[2]) else day_end
      # run the validity tests (except include_sunset, handled above)
      validity <- mm_is_valid_day(
        night_dat, # data split by mm_model_by_ply and subsetted here
        day_start=night_start, day_end=night_end, day_tests=night_tests,
        ply_date=ply_date) # args passed from metab_night (after modifying the tests)
      if(!isTRUE(validity)) stop_strs <- c(stop_strs, validity)
      
      # actually stop if anything has broken so far; the tryCatch will catch it,
      # and our stop_strs will be retained for later reporting
      if(length(stop_strs) > 0) stop("")
      
      # smooth DO data
      night_dat$DO.obs.smooth <- u(c(stats::filter(night_dat$DO.obs, rep(1/3, 3), sides=2)), get_units(night_dat$DO.obs))
      
      # calculate dDO/dt
      dDO <- u(diff(v(night_dat$DO.obs.smooth)), get_units(night_dat$DO.obs.smooth))
      dt <- u(as.numeric(diff(v(night_dat$solar.time)), units='days'), 'd')
      night_dat$dDO.dt <- (dDO / dt)[c(NA, 1:length(dDO))]
      
      # calculate saturation deficit
      night_dat$DO.sat.def <- night_dat$DO.sat - night_dat$DO.obs.smooth
      
      # fit model & extract the important stuff (see Chapra & DiToro 1991)
      lm_dDOdt <- lm(dDO.dt ~ DO.sat.def, data=v(night_dat))
      KO2_units <-get_units(night_dat$dDO.dt[1] / night_dat$DO.sat.def[1])
      out <- list(
        row.first = which_night[1],
        row.last = which_night[length(which_night)],
        KO2 = u(coef(lm_dDOdt)[["DO.sat.def"]], KO2_units),
        KO2.sd = u(summary(lm_dDOdt)$coefficients[["DO.sat.def","Std. Error"]], KO2_units),
        ER.vol = u(coef(lm_dDOdt)[["(Intercept)"]], get_units(night_dat$dDO.dt)),
        ER.sd.vol = u(summary(lm_dDOdt)$coefficients[["(Intercept)","Std. Error"]], get_units(night_dat$dDO.dt)),
        rsq = summary(lm_dDOdt)$r.squared,
        p = summary(lm_dDOdt)$coefficients[["DO.sat.def","Pr(>|t|)"]]
      )
      # convert ER from volumetric to area-based
      mean_depth <- mean(night_dat$depth)
      out <- c(out, list(
        ER = out$ER.vol * mean_depth,
        ER.sd = out$ER.sd.vol * mean_depth
      ))
      
      # convert from KO2 to K600. the coefficients used here differ for 
      # Maite's code and that in LakeMetabolizer. Maite and Bob use K600 = 
      # ((600/(1800.6-(120.1*temp)+(3.7818*temp^2)-(0.047608*temp^3)))^-0.5)*KO2.
      # LakeMetabolizer & streamMetabolizer use K600 <- ((600/(1568 - 
      # 86.04*temp + 2.142*temp^2 - 0.0216*temp^3))^-0.5)*KO2. The resulting 
      # numbers are believable by either method but do differ. I'll stick with
      # LakeMetabolizer since those are the coefficients from Raymond et al. 
      # 2012, which is probably the newer source. the sd conversion works
      # because the conversion is linear in KO2.
      out <- c(out, list(
        K600 = convert_kGAS_to_k600(kGAS=out$KO2, temperature=u(mean(night_dat$temp.water), 'degC'), gas="O2"),
        K600.sd = convert_kGAS_to_k600(kGAS=out$KO2.sd, temperature=u(mean(night_dat$temp.water), 'degC'), gas="O2")
      ))
      
      # return everything. remove units since we can't fully support them
      # elsewhere yet
      lapply(out, v)
      
    }, error=function(err) {
      # on error: give up, remembering error
      if(nchar(err$message) > 0) stop_strs <<- c(stop_strs, err$message)
      NA
    }), warning=function(war) {
      # on warning: record the warning and run again
      warn_strs <<- c(warn_strs, war$message)
      invokeRestart("muffleWarning")
    })
  
  # stop_strs may have accumulated during nlm() call. If failed, use dummy data 
  # to fill in the model output with NAs.
  if(length(stop_strs) > 0) {
    night.1d <- setNames(as.list(rep(NA, 10)), c('KO2','KO2.sd','ER','ER.sd','rsq','p','K600','K600.sd','row.first','row.last'))
  }
  
  # Return, reporting any results, warnings, and errors
  data.frame(GPP=NA, GPP.sd=NA, 
             ER=night.1d$ER, ER.sd=night.1d$ER.sd, 
             K600=night.1d$K600, K600.sd=night.1d$K600.sd,
             r.squared=night.1d$rsq, p.value=night.1d$p,
             row.first=night.1d$row.first, row.last=night.1d$row.last,
             warnings=paste0(unique(warn_strs), collapse="; "), 
             errors=paste0(unique(stop_strs), collapse="; "),
             stringsAsFactors=FALSE)
}



#### metab_night class ####

#' Reaeration model fitted by nighttime regression
#' 
#' \code{metab_night} models use nighttime regression to fit values of K for a
#' given DO time series.
#' 
#' @exportClass metab_night
#' @family metab.model.classes
setClass(
  "metab_night", 
  contains="metab_model"
)


#' Nighttime dissolved oxygen predictions from a metab_night.
#' 
#' metab_night only fits ER and K, and only for the darkness hours. We will
#' therefore make predictions just for those hours.
#' 
#' @inheritParams predict_DO
#' @return A data.frame of predictions, as for the generic 
#'   \code{\link{predict_DO}}.
#' @export
#' @import dplyr
#' @family predict_DO
predict_DO.metab_night <- function(metab_model, date_start=NA, date_end=NA, ..., use_saved=TRUE) {
  
  # pull args from the model
  day_start <- get_specs(metab_model)$day_start
  day_end <- get_specs(metab_model)$day_end
  
  # get the DO, temperature, etc. data; filter if requested
  data <- get_data(metab_model) %>%
    mm_filter_dates(date_start=date_start, date_end=date_end, day_start=day_start, day_end=day_end)
  
  # if allowed and available, use previously stored values for DO.mod rather than re-predicting them now
  if(isTRUE(use_saved)) {
    if(!is.null(data) && "DO.mod" %in% names(data)) {
      return(data)
    }
  }

  # get the metabolism (GPP, ER) data and estimates; filter if requested
  date <- ER <- K600 <- row.first <- row.last <- ".dplyr.var"
  metab_ests <- get_fit(metab_model) %>% 
    dplyr::select(date, ER, K600, row.first, row.last) %>%
    mm_filter_dates(date_start=date_start, date_end=date_end) %>%
    mutate(GPP = 0)
  
  # re-process the input data with the metabolism estimates to predict DO, using
  # our special nighttime regression prediction function
  mm_model_by_ply(
    model_fun=metab_night_predict_1ply, data=data, data_daily=metab_ests, # for mm_model_by_ply
    day_start=day_start, day_end=day_end, day_tests=c()) %>% # for mm_model_by_ply
    mm_filter_dates(date_start=date_start, date_end=date_end)
}

#' Helper to predict_DO.metab_model
#' 
#' Usually assigned to model_fun within mm_model_by_ply, called from there
#' 
#' @inheritParams mm_model_by_ply_prototype
#' @return a data.frame of predictions
#' @importFrom stats complete.cases
metab_night_predict_1ply <- function(
  data_ply, data_daily_ply, day_start, day_end, ply_date, timestep_days, ... # inheritParams mm_model_by_ply_prototype
) {
  
  # subset to times of darkness, just as we did in nightreg_1ply
  if(nrow(data_daily_ply)>0 && !is.na(data_daily_ply$row.first) && !is.na(data_daily_ply$row.last)) {
    which_night <- seq(min(nrow(data_ply), data_daily_ply$row.first), 
                       min(nrow(data_ply), data_daily_ply$row.last))
    night_dat <- data_ply[which_night,]
  } else {
    # for metab_night, sometimes data_plys (especially night_dat) are entirely empty. if that's today,
    # return right quick now.
    night_dat <- data_ply[c(),]
    night_dat[1,'DO.mod'] <- as.numeric(NA)
    return(night_dat)
  }
  
  # apply the regular prediction function
  mm_predict_1ply(
    data_ply=night_dat, data_daily_ply=data_daily_ply, 
    day_start, day_end, ply_date, timestep_days=timestep_days, 
    calc_DO_fun=calc_DO_mod) # use default ODE_method
  
}