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
#' # set the date in several formats
#' start.chron <- chron::chron(dates="08/23/12", times="22:00:00")
#' end.chron <- chron::chron(dates="08/25/12", times="06:00:00")
#' start.posix <- as.POSIXct(format(start.chron, "%Y-%m-%d %H:%M:%S"), tz="UTC")
#' end.posix <- as.POSIXct(format(end.chron, "%Y-%m-%d %H:%M:%S"), tz="UTC")
#' mid.date <- as.Date(start.posix + (end.posix - start.posix)/2, tz=lubridate::tz(start.posix))
#' 
#' # get, format, & subset data
#' vfrench <- streamMetabolizer:::load_french_creek(attach.units=FALSE)
#' vfrenchshort <- vfrench[vfrench$solar.time >= start.posix & vfrench$solar.time <= end.posix, ]
#' 
#' # dates & subsetting specific to nighttime regression
#' first.dark <- 100 + which(vfrenchshort$light[101:nrow(vfrenchshort)] < 0.1)[1]
#' stop.dark <- 100 + which(
#'   format(vfrenchshort$solar.time[101:nrow(vfrenchshort)], "%H:%M") == "23:00")[1]
#' vfrenchnight <- vfrenchshort[first.dark:stop.dark,]
#' night.start <- eval(parse(text=format(vfrenchnight$solar.time[1], "%H + %M/60")))
#' night.end <- eval(parse(text=format(vfrenchnight$solar.time[nrow(vfrenchnight)], "%H + %M/60")))
#' 
#' # fit
#' mm <- metab_night(data=vfrenchnight, 
#'   specs=specs('n_np_pi_eu_.lm'), 
#'   day_start=night.start, day_end=night.end)
#'   
#' # give estimates
#' predict_metab(mm)
#' streamMetabolizer:::load_french_creek_std_mle(vfrenchnight, estimate='K')
#' 
#' # predict DO
#' plot_DO_preds(predict_DO(mm))
#' 
#' @export
#' @family metab_model
metab_night <- function(
  specs=specs(mm_name('night')),
  data=mm_data(solar.time, DO.obs, DO.sat, depth, temp.water, light), 
  data_daily=mm_data(NULL),
  info=NULL, 
  day_start=12, day_end=36, tests=c('full_day', 'even_timesteps', 'complete_data')
) {
  
  fitting_time <- system.time({
    # Check data for correct column names & units
    dat_list <- mm_validate_data(data, if(missing(data_daily)) NULL else data_daily, "metab_night")
    data <- dat_list[['data']]
    data_daily <- dat_list[['data_daily']]
    
    # model the data, splitting into potentially overlapping 'plys' for each date
    night_all <- mm_model_by_ply(
      nightreg_1ply, data=data, data_daily=data_daily, # for mm_model_by_ply
      day_start=day_start, day_end=day_end, # for mm_model_by_ply and mm_is_valid_day
      tests=tests, # for mm_is_valid_day
      specs=specs) # for nightreg_1ply
  })
  
  # Package and return results
  mm <- metab_model(
    model_class="metab_night", 
    info=info,
    fit=night_all,
    fitting_time=fitting_time,
    args=list(specs=specs, day_start=day_start, day_end=day_end, tests=tests),
    data=data,
    data_daily=data_daily)
  
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
#' @inheritParams mm_is_valid_day
#' @inheritParams metab
#' @return data.frame of estimates and \code{\link[stats]{nlm}} model 
#'   diagnostics
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
  data_ply, data_daily_ply, day_start, day_end, ply_date, # inheritParams mm_model_by_ply_prototype
  tests=c('full_day', 'even_timesteps', 'complete_data'), # inheritParams mm_is_valid_day
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
        stop_strs <<- c(stop_strs, "need exactly one night per data_ply")
      night_dat <- data_ply[which_night,]
      # do the full_day test here, rather than in mm_is_valid_day, because
      # it's not about specific times but about whether night begins and ends
      # within the ply bounds
      if(!is.na(full_day_test <- match('full_day',tests))) {
        which_twilight <- c(head(which_night,1) - 1, tail(which_night,1) + 1)
        if(which_twilight[1] < 1)
          stop_strs <<- c(stop_strs, "data don't start before day-night transition")
        if(which_twilight[2] > nrow(data_ply))
          stop_strs <<- c(stop_strs, "data don't end after night-day transition")
        tests <- tests[-full_day_test]
      }
      
      # check for data validity here, after subsetting to the time window we
      # really care about
      timestep.days <- suppressWarnings(mean(as.numeric(diff(v(night_dat$solar.time)), units="days"), na.rm=TRUE))
      validity <- mm_is_valid_day(
        night_dat, # data split by mm_model_by_ply and subsetted here
        tests=tests, day_start=day_start, day_end=day_end, # args passed from metab_night (after modifying tests)
        timestep_days=timestep.days) # arg supplied here to avoid calculating twice
      if(!isTRUE(validity)) stop_strs <<- c(stop_strs, validity)
      
      # actually stop if anything has broken so far
      if(length(stop_strs) > 0) stop("invalid data ply")
      
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
        K600 = convert_kGAS_to_k600(kGAS=out$KO2, temperature=mean(v(night_dat$temp.water)), gas="O2"),
        K600.sd = convert_kGAS_to_k600(kGAS=out$KO2.sd, temperature=mean(v(night_dat$temp.water)), gas="O2")
      ))
      
      # return everything. remove units since we can't fully support them
      # elsewhere yet
      lapply(out, v)
      
    }, error=function(err) {
      # on error: give up, remembering error
      stop_strs <<- c(stop_strs, err$message)
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
             warnings=paste0(warn_strs, collapse="; "), 
             errors=paste0(stop_strs, collapse="; "),
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
  day_start <- get_args(metab_model)$day_start
  day_end <- get_args(metab_model)$day_end
  
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
    day_start=day_start, day_end=day_end) %>% # for mm_model_by_ply
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
  data_ply, data_daily_ply, day_start, day_end, ply_date # inheritParams mm_model_by_ply_prototype
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
  mm_predict_1ply(data_ply=night_dat, data_daily_ply=data_daily_ply, 
                  day_start, day_end, ply_date, calc_DO_fun=calc_DO_mod) # use default ODE_method
  
}