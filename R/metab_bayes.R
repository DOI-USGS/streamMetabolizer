#' @include metab_model-class.R
NULL

#' Basic Bayesian metabolism model fitting function
#' 
#' Fits a Bayesian model to estimate GPP and ER from input data on DO, 
#' temperature, light, etc. See \code{\link{mm_name}} to choose a Bayesian model
#' and \code{\link{specs}} for relevant options for the \code{specs} argument.
#' 
#' As of summer and fall 2016, a new compilation of any Stan model gives
#' deprecation warnings including \code{typedef 'size_type' locally defined but
#' not used [-Wunused-local-typedefs]}, \code{typedef 'index_range' locally
#' defined but not used [-Wunused-local-typedefs]}, \code{typedef 'index'
#' locally defined but not used [-Wunused-local-typedefs]}, and \code{'void
#' stan::math::set_zero_all_adjoints()' defined but not used
#' [-Wunused-function]}. THESE ARE OKAY. Subsequent runs of the compiled Stan
#' model will be quieter, and the model will work.
#' 
#' @author Alison Appling, Bob Hall
#'   
#' @inheritParams metab
#' @return A metab_bayes object containing the fitted model. This object can be 
#'   inspected with the functions in the \code{\link{metab_model_interface}} and
#'   also \code{\link{get_mcmc}}.
#'   
#' @examples
#' \dontrun{
#' dat <- data_metab('3', res='30')
#' # fast-ish model version, but still too slow to auto-run in examples
#' mm <- metab_bayes(data=dat,
#'   specs(mm_name('bayes', err_proc_iid=FALSE), 
#'     n_cores=3, n_chains=3, burnin_steps=300, saved_steps=100))
#' mm
#' get_fitting_time(mm)
#' predict_metab(mm)
#' plot_DO_preds(predict_DO(mm))
#' 
#' # error and warning messages are printed with the mm object if present
#' dat <- data_metab('3', res='30', flaws=c('missing middle'))
#' mm <- metab(specs(mm_name('bayes', err_proc_iid=FALSE), 
#'   n_cores=3, n_chains=3, burnin_steps=300, saved_steps=100, verbose=FALSE),
#'   data=dat)
#' predict_metab(mm)
#' }
#' @export
#' @family metab_model
metab_bayes <- function(
  specs=specs(mm_name('bayes')), 
  data=mm_data(solar.time, DO.obs, DO.sat, depth, temp.water, light, discharge, optional='discharge'), 
  data_daily=mm_data(date, discharge.daily, optional='all'),
  info=NULL
) {
  
  if(missing(specs)) {
    # if specs is left to the default, it gets confused about whether specs() is
    # the argument or the function. tell it which:
    specs <- streamMetabolizer::specs(mm_name('bayes'))
  }
  fitting_time <- system.time({
    # Check data for correct column names & units
    dat_list <- mm_validate_data(if(missing(data)) NULL else data, if(missing(data_daily)) NULL else data_daily, "metab_bayes")
    num_discharge_cols <- length(grep('discharge', c(names(dat_list$data), names(dat_list$data_daily))))
    pool_K600 <- mm_parse_name(specs$model_name)$pool_K600
    if(xor(num_discharge_cols > 0, pool_K600 %in% c('linear','binned'))) 
      stop('discharge data should be included if & only if pool_K600 indicates hierarchy')
    if(num_discharge_cols > 1)
      stop('either discharge or discharge.daily may be specified, but not both')
    
    # Handle discharge. If K600 is a hierarchical function of discharge and 
    # data$discharge was given, compute daily discharge and store in data.daily 
    # where it'll be accessible for the user to inspect it after model fitting
    if((pool_K600 %in% c('linear','binned')) && ('discharge' %in% names(dat_list$data))) {
      # calculate daily discharge
      dailymean <- function(data_ply, data_daily_ply, day_start, day_end, ply_date, ply_validity, timestep_days, ...) {
        data.frame(discharge.daily = if(isTRUE(ply_validity[1])) mean(data_ply$discharge) else NA)
      }
      dischdaily <- mm_model_by_ply(
        model_fun=dailymean, data=v(dat_list$data), day_start=specs$day_start, day_end=specs$day_end)
      
      # add units if either of the input dfs were unitted
      if(is.unitted(dat_list$data) || is.unitted(dat_list$data_daily)) {
        dischdaily_units <- get_units(mm_data(date, discharge.daily))
        if(is.unitted(dat_list$data)) {
          if(dischdaily_units['discharge.daily'] != get_units(dat_list$data$discharge))
            stop('mismatch between discharge units in data (',get_units(dat_list$data$discharge),
                 ') and requirement for data_daily (', dischdaily_units['discharge.daily'] ,')')
        }
        dischdaily <- u(dischdaily, dischdaily_units)
      }

      # merge with any existing dat_list$data_daily      
      if(is.null(v(dat_list$data_daily))) {
        dat_list$data_daily <- dischdaily
      } else {
        # need both or neither dfs to be unitted for the full_join. if 
        # data_daily was unitted then dischdaily will already also be unitted 
        # (see add units chunk above), so just check/fix the case where
        # data_daily lacks units but data & therefore dischdaily has them
        if(is.unitted(dischdaily) && !is.unitted(dat_list$data_daily)) {
          dat_list$data_daily <- u(dat_list$data_daily, get_units(mm_data()[names(dat_list$data_daily)]))
        }
        dat_list$data_daily <- full_join(dat_list$data_daily, dischdaily, by='date')
      }
    }
    # If we have discharge.daily, then we need logged discharge.daily. compute
    # and store it now
    if('discharge.daily' %in% names(dat_list$data_daily)) {
      dat_list$data_daily$ln.discharge.daily <- log(v(dat_list$data_daily$discharge.daily))
    }
    # If we need discharge bins, compute & store those now, as well
    if(pool_K600 %in% c('binned')) {
      if(is.character(specs$K600_daily_beta_cuts)) {
        if(!requireNamespace('ggplot2', quietly=TRUE)) {
          stop("need ggplot2 for K600_pool='binned' when is.character(K600_daily_beta_cuts). ",
               "either install the ggplot2 package or switch to a numeric vector for K600_daily_beta_cuts")
        }
        cut_fun <- switch(
          specs$K600_daily_beta_cuts,
          interval = ggplot2::cut_interval,
          number = ggplot2::cut_number)
        # run once with high dig.lab to parse the breaks from levels(cuts) as numeric
        cuts <- cut_fun(v(dat_list$data_daily$ln.discharge.daily), n=specs$K600_daily_beta_num, dig.lab=20)
        specs$K600_daily_beta_breaks <- levels(cuts) %>%
          strsplit('\\[|\\(|\\]|,') %>%
          lapply(function(lev) as.numeric(lev[2:3])) %>%
          unlist() %>%
          unique()
        # run again with default dig.lab for prettier bin labels
        cuts <- cut_fun(v(dat_list$data_daily$ln.discharge.daily), n=specs$K600_daily_beta_num)
      } else {
        # we already know breaks precisely because we're supplying them
        specs$K600_daily_beta_breaks <- specs$K600_daily_beta_cuts
        cuts <- cut(dat_list$data_daily$ln.discharge.daily, breaks=specs$K600_daily_beta_cuts, ordered_result=TRUE)
      }
      dat_list$data_daily$discharge.bin.daily <- as.numeric(cuts)
      specs$K600_daily_beta_bins <- levels(cuts)
      # you should be able to retrieve the bin names with 
      # specs$K600_daily_beta_bins[dat_list[['data_daily']]$discharge.bin.daily]
      # or, later, 
      # get_specs(fit)$K600_daily_beta_bins[get_data_daily(fit)$discharge.bin.daily].
      # specs$K600_daily_beta_breaks and specs$K600_daily_beta_bins are created
      # here purely for manual inspection after the model has been run
    }
    
    # Use de-unitted version until we pack up the model to return
    data <- v(dat_list$data)
    data_daily <- v(dat_list$data_daily)
        
    # Check and parse model file path
    specs$model_path <- mm_locate_filename(specs$model_name)
    
    # check the format of keep_mcmcs (more checks, below, are split_dates-specific)
    if(is.logical(specs$keep_mcmcs)) {
      if(length(specs$keep_mcmcs) != 1) {
        stop("if keep_mcmcs is logical, it must have length 1")
      }
    } else if(specs$split_dates == FALSE) {
      stop("if split_dates==FALSE, keep_mcmcs must be a single logical value")
    }
    if(is.logical(specs$keep_mcmc_data)) {
      if(length(specs$keep_mcmc_data) != 1) {
        stop("if keep_mcmc_data is logical, it must have length 1")
      }
    } else if(specs$split_dates == FALSE) {
      stop("if split_dates==FALSE, keep_mcmc_data must be a single logical value")
    }
    
    # model the data. create outputs bayes_all (a data.frame) and bayes_mcmc (an MCMC object from tan)
    if(specs$split_dates == TRUE) {
      if(!is.logical(specs$keep_mcmcs)) specs$keep_mcmcs <- as.Date(specs$keep_mcmcs)
      if(!is.logical(specs$keep_mcmc_data)) specs$keep_mcmc_data <- as.Date(specs$keep_mcmc_data)
      # one day at a time, splitting into overlapping ~24-hr 'plys' for each date
      bayes_daily <- mm_model_by_ply(
        bayes_1ply, data=data, data_daily=data_daily, # for mm_model_by_ply
        day_start=specs$day_start, day_end=specs$day_end, day_tests=specs$day_tests, # for mm_model_by_ply
        specs=specs) # for bayes_1ply
      # if we saved the modeling object[s] in the df, pull them out now
      extract_object_list <- function(col) {
        # has side effects on bayes_daily!
        if(col %in% names(bayes_daily)) {
          val <- bayes_daily[[col]]
          names(val) <- bayes_daily$date
          bayes_daily[[col]] <<- NULL
          val
        } else NULL
      }
      . <- '.dplyr.var'
      compile_log <- extract_object_list('compile_log') %>% 
      {.[!sapply(., is.null)]} %>% 
      {if(!is.null(.)) setNames(., 'Compilation') else .}
      log <- extract_object_list('log') %>% { setNames(., paste0('MCMC_', names(.))) }
      bayes_log <- c(compile_log, log)
      bayes_mcmc <- extract_object_list('mcmcfit')
      bayes_mcmc_data <- extract_object_list('mcmc_data')
      bayes_all <- list(daily=bayes_daily)
      if(nrow(bayes_all$daily) == 0 || length(which(bayes_daily$valid_day)) == 0) 
        bayes_all$errors <- c(bayes_all$errors, "no valid days of data")
      
    } else if(specs$split_dates == FALSE) {
      # all days at a time, after first filtering out bad days
      filtered <- mm_filter_valid_days(
        data, data_daily, day_start=specs$day_start, day_end=specs$day_end, day_tests=specs$day_tests)
      if(length(unique(filtered$data$date)) > 1 && (specs$day_end - specs$day_start) > 24) 
        warning("multi-day models should probably have day_end - day_start <= 24 hours")
      bayes_all_list <- bayes_allply(
        data_all=filtered$data, data_daily_all=filtered$data_daily, removed=filtered$removed,
        specs=specs)
      # if we saved the modeling object[s] in the list, pull them out now
      . <- '.dplyr.var'
      bayes_log <- bayes_all_list[c('compile_log', 'log')] %>% 
        setNames(c('Compilation','MCMC_All_Days')) %>% { .[!sapply(., is.null)] }
      bayes_mcmc <- bayes_all_list$mcmcfit
      bayes_mcmc_data <- bayes_all_list$mcmc_data
      # now a list of dfs, log, warnings, and errors
      bayes_all <- bayes_all_list[!(names(bayes_all_list) %in% c('compile_log','log','mcmcfit','mcmc_data'))]
    }
  })
  
  # Package and return results
  mm <- metab_model(
    model_class="metab_bayes", 
    info=info,
    fit=bayes_all,
    log=bayes_log,
    mcmc=bayes_mcmc,
    mcmc_data=bayes_mcmc_data,
    fitting_time=fitting_time,
    specs=specs,
    data=dat_list$data, # keep the units if given
    data_daily=dat_list$data_daily)
  
  # Update data with DO predictions
  core_cols <- grepl("^(date|GPP|ER|K600)", names(bayes_all$daily))
  success <- any(complete.cases(bayes_all$daily[core_cols])) && 
    length(bayes_all$errors) == 0 &&
    length(bayes_all$fit$errors[bayes_all$fit$valid_day] != '') == 0
  if(success) {
    mm@data <- predict_DO(mm)
  } else {
    warning(paste0('Modeling failed: ', paste0(bayes_all$errors, collapse='\n')))
  }
  
  # Return
  mm
}


#### helpers ####

#' Make daily metabolism estimates from input parameters
#' 
#' Called from metab_bayes().
#' 
#' @inheritParams mm_model_by_ply_prototype
#' @inheritParams metab
#' @return data.frame of estimates and MCMC model diagnostics
#' @importFrom stats setNames
#' @keywords internal
bayes_1ply <- function(
  data_ply, data_daily_ply, ply_date, ply_validity, ..., # inheritParams mm_model_by_ply_prototype
  specs # inheritParams metab
) {
  
  # Provide ability to skip a poorly-formatted day for calculating 
  # metabolism, without breaking the whole loop. Just collect 
  # problems/errors as a list of strings and proceed. Also collect warnings.
  stop_strs <- if(isTRUE(ply_validity)) character(0) else ply_validity
  warn_strs <- character(0)
  
  specs$keep_mcmc <- if(is.logical(specs$keep_mcmcs)) {
    isTRUE(specs$keep_mcmcs)
  } else {
    isTRUE(ply_date %in% specs$keep_mcmcs)
  }
  
  # Calculate metabolism by Bayesian MCMC
  data_list <- NULL # (in case it doesn't get assigned in the tryCatch)
  if(length(stop_strs) == 0) {
    bayes_1day <- withCallingHandlers(
      tryCatch({
        # first: try to run the bayes fitting function
        data_list <- prepdata_bayes(
          data=data_ply, data_daily=data_daily_ply, ply_date=ply_date,
          specs=specs, engine=specs$engine, model_name=specs$model_name)
        do.call(mcmc_bayes, c(list(data_list=data_list), specs))
      }, error=function(err) {
        # on error: give up, remembering error. dummy values provided below
        stop_strs <<- c(stop_strs, err$message)
        NA
      }), warning=function(war) {
        # on warning: record the warning and run again
        warn_strs <<- c(warn_strs, war$message)
        invokeRestart("muffleWarning")
      })
  } 

  # stop_strs may have accumulated during prepdata_bayes() or mcmc_bayes()
  # calls. If failed, use dummy data to fill in the model output with NAs.
  if(length(stop_strs) > 0) {
    bayes_1day <- data.frame(GPP_daily_2.5pct=NA, GPP_daily_50pct=NA, GPP_daily_97.5pct=NA,
                             ER_daily_2.5pct=NA, ER_daily_50pct=NA, ER_daily_97.5pct=NA, 
                             K600_daily_2.5pct=NA, K600_daily_50pct=NA, K600_daily_97.5pct=NA)
  }
  
  # package the results, data, warnings, and errors
  outdf <- data.frame(
    bayes_1day[!(names(bayes_1day) %in% c('mcmcfit','log','compile_log'))],
    valid_day=isTRUE(ply_validity),
    warnings=paste0(trimws(unique(warn_strs)), collapse="; "), 
    errors=paste0(trimws(unique(stop_strs)), collapse="; "),
    stringsAsFactors=FALSE) %>%
    mutate(log = list(bayes_1day$log))
  
  # attach the compile_log, mcmcfit, & mcmcdata if requested/available
  if(exists('compile_log', bayes_1day)) outdf$compile_log <- list(bayes_1day$compile_log)
  if(specs$keep_mcmc) outdf$mcmcfit <- bayes_1day$mcmcfit
  keep_mcmc_dat <- 
    if(is.logical(specs$keep_mcmc_data)) {
      isTRUE(specs$keep_mcmc_data)
    } else {
      isTRUE(ply_date %in% specs$keep_mcmc_data)
    }
  if(keep_mcmc_dat) outdf$mcmc_data <- list(data_list)
  
  # return
  outdf
}


#' Make daily metabolism estimates from input parameters using a hierarchical 
#' approach.
#' 
#' Called from metab_bayes().
#' 
#' @param data_all data.frame of the form \code{mm_data(solar.time, DO.obs, 
#'   DO.sat, depth, temp.water, light)} and containing data for just one 
#'   estimation-day (this may be >24 hours but only yields estimates for one 
#'   24-hour period)
#' @param data_daily_all data.frame of daily priors, if appropriate to the given
#'   model_path
#' @param removed data.frame of dates that were removed and why
#' @inheritParams metab
#' @return data.frame of estimates and MCMC model diagnostics
#' @keywords internal
bayes_allply <- function(
  data_all, data_daily_all, removed,
  specs
) {
  # Provide ability to skip a poorly-formatted dataset for calculating 
  # metabolism. Collect problems/errors as a list of strings and proceed. Also
  # collect warnings.
  stop_strs <- warn_strs <- character(0)
  
  # Calculate metabolism by Bayesian MCMC
  data_list <- NULL # (in case it doesn't get assigned in the tryCatch)
  bayes_allday <- withCallingHandlers(
    tryCatch({
      if(is.null(data_all) || nrow(data_all) == 0) stop("no valid days of data")
      # first: try to run the bayes fitting function
      data_list <- prepdata_bayes(
        data=data_all, data_daily=data_daily_all, ply_date=NA,
        specs=specs, engine=specs$engine, model_name=specs$model_name)
      specs$keep_mcmc <- specs$keep_mcmcs
      do.call(mcmc_bayes, c(list(data_list=data_list), specs))
    }, error=function(err) {
      # on error: give up, remembering error. dummy values provided below
      stop_strs <<- c(stop_strs, err$message)
      NA
    }), warning=function(war) {
      # on warning: record the warning and run again
      warn_strs <<- c(warn_strs, war$message)
      invokeRestart("muffleWarning")
    })

  # match dates back to daily estimates
  date_vec <- unique(data_all$date)
  
  # stop_strs may have accumulated during prepdata_bayes() or mcmc_bayes()
  # calls. If failed, use dummy data to fill in the model output with NAs.
  if(length(stop_strs) > 0 || any(grepl("^Stan model .* does not contain samples", warn_strs))) {
    na_vec <- rep(as.numeric(NA), length(date_vec))
    bayes_allday <- c(
      list(daily=data.frame(date=date_vec, GPP_daily_2.5pct=na_vec, GPP_daily_50pct=na_vec, GPP_daily_97.5pct=na_vec,
                            ER_daily_2.5pct=na_vec, ER_daily_50pct=na_vec, ER_daily_97.5pct=na_vec, 
                            K600_daily_2.5pct=na_vec, K600_daily_50pct=na_vec, K600_daily_97.5pct=na_vec)),
      list(log=if(exists('bayes_allday') && is.list(bayes_allday)) {
        bayes_allday[c('compile_log', 'log')]
      } else NULL ))
  } else {
    # check this now in case we're not saving the data_list
    if(length(date_vec) != data_list$d || length(date_vec) != nrow(bayes_allday$daily))
      stop_strs <- c(stop_strs, "couldn't match dates to date indices")
    index <- '.dplyr.var'
    bayes_allday$daily <- bayes_allday$daily %>%
      rename(date=index) %>% 
      mutate(date=date_vec)
  }

  # add columns for compatibility with other dfs (e.g., removed) & methods
  # (e.g., predict_metab)
  bayes_allday$daily <- bayes_allday$daily %>%
    mutate(valid_day=TRUE, warnings='', errors='')
  
  # add back the dates that were removed during date filtering
  if(nrow(removed) > 0) {
    if(nrow(bayes_allday$daily) > 0) {
      bayes_allday$daily <- bayes_allday$daily %>%
        full_join(mutate(removed, valid_day=FALSE, warnings=''), by=c('date', 'valid_day', 'warnings', 'errors')) %>%
        arrange(date)
    } else {
      GPP_daily_2.5pct <- GPP_daily_50pct <- GPP_daily_97.5pct <- ER_daily_2.5pct <- ER_daily_50pct <- ER_daily_97.5pct <- 
        K600_daily_2.5pct <- K600_daily_50pct <- K600_daily_97.5pct <- valid_day <- warnings <- errors <- '.dplyr.var'
      bayes_allday$daily  <- 
        mutate(removed, 
               GPP_daily_2.5pct=NA, GPP_daily_50pct=NA, GPP_daily_97.5pct=NA,
               ER_daily_2.5pct=NA, ER_daily_50pct=NA, ER_daily_97.5pct=NA, 
               K600_daily_2.5pct=NA, K600_daily_50pct=NA, K600_daily_97.5pct=NA,
               valid_day=FALSE, warnings='') %>% 
        select(date, GPP_daily_2.5pct, GPP_daily_50pct, GPP_daily_97.5pct,
               ER_daily_2.5pct, ER_daily_50pct, ER_daily_97.5pct, 
               K600_daily_2.5pct, K600_daily_50pct, K600_daily_97.5pct,
               valid_day, warnings, errors)
    }
  }
  
  # Return, reporting any results, warnings, and errors
  c(bayes_allday,
    list(mcmc_data=if(specs$keep_mcmc_data) data_list else NULL,
         warnings=trimws(unique(warn_strs)),
         errors=trimws(unique(stop_strs))))
}


#### helpers to the helper ####

#' Prepare data for passing to Stan
#' 
#' This function accepts exactly one day's worth of data, (one ply, which might 
#' be 24 hrs or 31.5 or so), which should already be validated. It prepares the 
#' data needed to run a Bayesian MCMC method to estimate GPP, ER, and K600.
#' 
#' @inheritParams mm_model_by_ply_prototype
#' @inheritParams metab
#' @inheritParams specs
#' @return list of data for input to runstan_bayes
#' @importFrom unitted v
#' @keywords internal
prepdata_bayes <- function(
  data, data_daily, ply_date=NA, # inheritParams mm_model_by_ply_prototype
  specs, # inheritParams metab (for hierarchical priors)
  engine, model_name #inheritParams specs
) {
  
  # remove units if present
  data <- v(data)
  data_daily <- v(data_daily)
  if(length(ply_date) != 1) stop("ply_date must have length 1")
  if(!is.na(ply_date)) data$date <- as.Date(ply_date)
  
  # define a function to package 1+ days of obs of a variable into a time x date matrix
  date_table <- table(data$date)
  num_dates <- length(date_table)
  num_daily_obs <- unique(unname(date_table))
  if(length(num_daily_obs) > 1) {
    warning(paste0(
      sapply(num_daily_obs, function(ndo) {
        tslabel <- paste(ndo, 'rows per day')
        tsdates <- names(date_table)[date_table == ndo]
        paste0(tslabel, ': ', paste0(tsdates, collapse=', '))
      }),
      collapse='\n')
    )
    stop("dates have differing numbers of rows; observations cannot be combined in matrix")
  }
  time_by_date_matrix <- function(vec) {
    matrix(data=vec, nrow=num_daily_obs, ncol=num_dates, byrow=FALSE)
  }
  
  # double-check that our dates are going to line up with the input dates. this 
  # should be redundant w/ above date_table checks, so just being extra careful
  obs_dates <- time_by_date_matrix(as.character(data$date, "%Y-%m-%d"))
  unique_dates <- apply(obs_dates, MARGIN=2, FUN=function(timevec) unique(timevec))
  if(!all.equal(unique_dates, names(date_table))) stop("couldn't fit given dates into matrix")
  
  # confirm that every day has the same modal timestep and put a value on that timestep
  obs_times <- time_by_date_matrix(as.numeric(data$solar.time - data$solar.time[1], units='days'))
  unique_timesteps <- unique(apply(obs_times, MARGIN=2, FUN=function(timevec) unique(round(diff(timevec), digits=12)))) # 10 digits is 8/1000000 of a second. 14 digits exceeds machine precision for datetimes
  if(length(unique_timesteps) != 1) stop("could not determine a single timestep for all observations")
  timestep_days <- mean(apply(obs_times, MARGIN=2, FUN=function(timevec) mean(diff(timevec))))
  
  # parse model name into features for deciding what data to include
  features <- mm_parse_name(model_name)
  
  # Format the data for Stan. Stan disallows period-separated names, so
  # change all the input data to underscore-separated. parameters given in
  # specs are already underscore-separated for this reason
  data_list = c(
    list(
      
      # Overall
      d = num_dates,
      
      # Daily
      n = num_daily_obs # one value applicable to every day
    ),
      
    switch(
      features$pool_K600,
      linear = list(ln_discharge_daily = data_daily$ln.discharge.daily),
      binned = list(
        b = specs$K600_daily_beta_num,
        discharge_bin_daily = data_daily$discharge.bin.daily)
    ),
    
    list(
      DO_obs_1 = array(time_by_date_matrix(data$DO.obs)[1,], dim=num_dates), # duplication of effort below should be small compared to MCMC time
      
      # Every timestep
      frac_GPP = {
        mat_light <- time_by_date_matrix(data$light)
        if(isTRUE(mm_parse_name(model_name)$GPP_fun == 'linlight')) {
          # normalize light by the sum of light in the first 24 hours of the time window
          in_solar_day <- apply(obs_times, MARGIN=2, FUN=function(timevec) {timevec - timevec[1] <= 1} )
          sweep(mat_light, MARGIN=2, STATS=colSums(mat_light*in_solar_day), FUN=`/`)
        } else {
          mat_light
        }
      },
      frac_ER  = time_by_date_matrix(timestep_days),
      frac_D   = time_by_date_matrix(timestep_days), # the yackulic shortcut models rely on this being constant over time
      KO2_conv = time_by_date_matrix(convert_k600_to_kGAS(k600=1, temperature=data$temp.water, gas="O2")),
      depth    = time_by_date_matrix(data$depth),
      DO_sat   = time_by_date_matrix(data$DO.sat),
      DO_obs   = time_by_date_matrix(data$DO.obs)
    ),
    
    specs[specs$params_in]
  )
  
  data_list
}

#' Run an MCMC simulation on a formatted data ply
#' 
#' @param data_list a formatted list of inputs to the Stan model
#' @param engine character string indicating which software to use
#' @param model_path the Stan model file to use, as a full file path
#' @param params_out a character vector of parameters whose values in the MCMC 
#'   runs should be recorded and summarized
#' @param keep_mcmc logical. If TRUE, the Stan output object will be saved. Be 
#'   careful; these can be big, and a run with many models might overwhelm R's 
#'   memory.
#' @param n_chains the number of chains to run
#' @param n_cores the number of cores to apply to this run
#' @param burnin_steps the number of steps per chain to run and ignore before 
#'   starting to collect MCMC 'data'
#' @param saved_steps the number of MCMC steps per chain to save
#' @param thin_steps the number of steps to move before saving another step. 1 
#'   means save all steps.
#' @param verbose logical. give status messages?
#' @param ... ignored arguments
#' @return a data.frame of outputs
#' @import parallel
#' @keywords internal
mcmc_bayes <- function(data_list, engine='stan', model_path, params_out, split_dates, keep_mcmc=FALSE, n_chains=4, n_cores=4, burnin_steps=4000, saved_steps=40000, thin_steps=1, verbose=FALSE, ...) {
  engine <- match.arg(engine)
  bayes_function <- switch(engine, stan = runstan_bayes)
  
  tot_cores <- detectCores()
  if (!is.finite(tot_cores)) { tot_cores <- 1 } 
  n_cores <- min(tot_cores, n_cores)
  if(verbose) message(paste0("MCMC (",engine,"): requesting ",n_chains," chains on ",n_cores," of ",tot_cores," available cores"))
  
  bayes_function(
    data_list=data_list, model_path=model_path, params_out=params_out, split_dates=split_dates, keep_mcmc=keep_mcmc, n_chains=n_chains, n_cores=n_cores, 
    burnin_steps=burnin_steps, saved_steps=saved_steps, thin_steps=thin_steps, verbose=verbose)
}

#' Run Stan on a formatted data ply
#' 
#' @inheritParams mcmc_bayes
#' @param ... args passed to other runxx_bayes functions but ignored here
#' @import parallel
#' @import dplyr
#' @import tibble
#' @importFrom tidyr gather spread
#' @keywords internal
runstan_bayes <- function(data_list, model_path, params_out, split_dates, keep_mcmc=FALSE, n_chains=4, n_cores=4, burnin_steps=1000, saved_steps=1000, thin_steps=1, verbose=FALSE, ...) {
  
  # stan() can't find its own function cpp_object_initializer() unless the 
  # namespace is loaded. requireNamespace is somehow not doing this. Thoughts
  # (not solution):
  # https://stat.ethz.ch/pipermail/r-devel/2014-September/069803.html
  if(!suppressPackageStartupMessages(require(rstan))) {
    stop("the rstan package is required for Stan MCMC models")
  }
  
  # use auto_write=TRUE to recompile if needed, or load from existing .rda file
  # without recompiling if possible
  mobj_path <- gsub('.stan$', '.stanrda', model_path)
  if(!file.exists(mobj_path) || file.info(mobj_path)$mtime < file.info(model_path)$mtime) {
    if(verbose) message("compiling Stan model")
    compile_log <- capture.output({
      stan_mobj <- rstan::stan_model(file=model_path, auto_write=TRUE)
    }, type=c('output'), split=verbose)
    rm(stan_mobj)
    gc() # this humble line saves us from many horrible R crashes
    autowrite_path <- gsub('.stan$', '.rda', model_path)
    if(!file.exists(autowrite_path)) autowrite_path <- file.path(tempdir(), basename(autowrite_path))
    if(!file.exists(autowrite_path)) {
      warning('could not find saved rda model file')
    } else {
      tryCatch({
        file.copy(autowrite_path, mobj_path, overwrite=TRUE)
        file.remove(autowrite_path)
      }, error=function(e) {
        warning('could not copy Stan rda to .stanrda file: ', e$message)
        mobj_path <- autowrite_path
      })
    }
  } else {
    if(verbose) message("loading pre-compiled Stan model")
  }
  stan_mobj <- readRDS(mobj_path)
  
  # make note of existing log files so we don't read them later
  oldlogfiles <- normalizePath(file.path(tempdir(), grep("_StanProgress.txt", dir(tempdir()), value=TRUE)))
  
  # run Stan
  if(verbose) message("sampling Stan model")
  consolelog <- capture.output(
    runstan_out <- rstan::sampling(
      object=stan_mobj,
      data=data_list,
      pars=params_out,
      include=TRUE,
      chains=n_chains,
      warmup=burnin_steps,
      iter=saved_steps+burnin_steps,
      thin=thin_steps,
      init="random",
      verbose=verbose,
      open_progress=FALSE,
      cores=n_cores),
    split=verbose)

  # this is a good place for a breakpoint when running small numbers of models
  # manually (or keep_mcmc also helps with inspection)
  #   show(runstan_out)
  #   rstan::plot(runstan_out)
  #   pairs(runstan_out)
  #   traceplot(runstan_out)
  
  # format output (but first detect and handle a failed model run)
  if(runstan_out@mode == 2L) {
    # for failed model runs, we still want to keep the mcmc
    stan_out <- NULL
    warning(capture.output(print(runstan_out)))
  } else if(split_dates) {
    # for one-day models, use a 1-row data.frame. see ls('package:rstan')
    stan_mat <- rstan::summary(runstan_out)$summary
    names_params <- rep(gsub("\\[1\\]", "", rownames(stan_mat)), each=ncol(stan_mat)) # the GPP, ER, etc. part of the name
    names_stats <- rep(gsub("%", "pct", colnames(stan_mat)), times=nrow(stan_mat)) # the mean, sd, etc. part of the name
    stan_out <- format_mcmc_mat_split(stan_mat, names_params, names_stats, keep_mcmc, runstan_out)
  } else {
    # for multi-day or unsplit models, format output into a list of data.frames,
    # one per unique number of nodes sharing a variable name
    stan_mat <- rstan::summary(runstan_out)$summary
    stan_out <- format_mcmc_mat_nosplit(stan_mat, data_list$d, keep_mcmc, runstan_out)
  } 
  
  # attach the contents of the most recent logfile in tempdir(), which should be for this model
  newlogfiles <- normalizePath(file.path(tempdir(), grep("_StanProgress.txt", dir(tempdir()), value=TRUE)))
  logfile <- setdiff(newlogfiles, oldlogfiles)
  log <- if(length(logfile) > 0) readLines(logfile) else consolelog
  stan_out <- c(stan_out, c(list(log=log), if(exists('compile_log')) list(compile_log=compile_log)))
  
  return(stan_out)
}

#' Format MCMC output into a one-row data.frame
#' 
#' For split_dates models. Formats output into a one-row data.frame for 
#' row-binding with other such data.frames
#' 
#' @param mcmc_mat matrix as extracted from Stan
#' @param names_params character vector of the names of the parameters
#' @param names_stats character vector of the names of the statistics
#' @import dplyr
#' @keywords internal
format_mcmc_mat_split <- function(mcmc_mat, names_params, names_stats, keep_mcmc, runmcmc_out) {
  # format the matrix into a 1-row df
  mcmc_out <- mcmc_mat %>% t %>% c %>% # get a 1D vector of GPP_daily_mean, GPP_sd, ..., ER_daily_mean, ER_daily_sd, ... etc
    t %>% as.data.frame() %>% # convert from 1D vector to 1-row data.frame
    setNames(paste0(names_params, "_", names_stats))
  
  # add the model object as a df column if requested
  if(keep_mcmc == TRUE) {
    mcmc_out <- mutate(mcmc_out, mcmcfit=list(runmcmc_out))
  }
  
  mcmc_out
}

#' Format MCMC output into a list of data.frames
#' 
#' For multi-day or unsplit models. Formats output into a list of data.frames, 
#' one per unique number of nodes sharing a variable name
#' 
#' @param mcmc_mat matrix as extracted from Stan
#' @import dplyr
#' @keywords internal
format_mcmc_mat_nosplit <- function(mcmc_mat, data_list_d, keep_mcmc, runmcmc_out) {
  stat <- val <- . <- rowname <- variable <- index <- varstat <- '.dplyr_var'
  
  colnames(mcmc_mat) <- gsub("%", "pct", colnames(mcmc_mat))
  
  # determine how many unique nrows, & therefore data.frames, there should be
  var_table <- table(gsub("\\[[[:digit:]]+\\]", "", rownames(mcmc_mat)))
  all_dims <- unique(var_table) %>% setNames(., .)
  names(all_dims)[all_dims == 1] <- "overall"
  names(all_dims)[all_dims == data_list_d] <- "daily" # overrides 'overall' if d==1. that's OK
  
  # for each unique nrows, create the data.frame with vars in columns and indices in rows
  mcmc_out <- lapply(all_dims, function(odim) {
    dim_params <- names(var_table[which(var_table == odim)])
    dim_rows <- sort(do.call(c, lapply(dim_params, function(dp) grep(paste0("^", dp, "(\\[|$)"), rownames(mcmc_mat)))))
    row_order <- names(sort(sapply(dim_params, function(dp) grep(paste0("^", dp, "(\\[|$)"), rownames(mcmc_mat))[1])))
    varstat_order <- paste0(rep(row_order, each=ncol(mcmc_mat)), '_', rep(colnames(mcmc_mat), times=length(row_order)))
    
    as.data.frame(mcmc_mat[dim_rows,,drop=FALSE]) %>% # stan version didn't have drop=FALSE
      rownames_to_column() %>%
      gather(stat, value=val, 2:ncol(.)) %>%
      mutate(variable=gsub("\\[[[:digit:]]+\\]", "", rowname),
             index=if(odim == 1) 1 else as.numeric(sapply(strsplit(rowname, "\\[|\\]"), `[[`, 2)),
             varstat=ordered(paste0(variable, '_', stat), varstat_order)) %>%
      select(index, varstat, val) %>%
      spread(varstat, val)
  })
  
  # add the model object as a list item if requested
  if(keep_mcmc == TRUE) {
    mcmc_out <- c(mcmc_out, list(mcmcfit=runmcmc_out))
  }
  
  mcmc_out
}


#### metab_bayes class ####

#' Metabolism model fitted by Bayesian MCMC
#' 
#' \code{metab_bayes} models use Bayesian MCMC methods to fit values of GPP, ER,
#' and K for a given DO curve.
#' 
#' @exportClass metab_bayes
#' @family metab.model.classes
setClass(
  "metab_bayes", 
  contains="metab_model",
  slots=c(log="ANY", mcmc="ANY", mcmc_data="ANY")
)

#' Extract any MCMC model objects that were stored with the model
#' 
#' A function specific to metab_bayes models. Returns an MCMC object or, for 
#' nopool models, a list of MCMC objects. These objects are saved by default
#' because they should usually be inspected manually; see \code{keep_mcmcs}
#' argument to \code{\link{specs}} for options for saving space.
#' 
#' @param metab_model A Bayesian metabolism model (metab_bayes) from which to 
#'   return the MCMC model object[s]
#' @return The MCMC model object[s]
#' @export
get_mcmc <- function(metab_model) {
  UseMethod("get_mcmc")
}

#' Retrieve any MCMC model object[s] that were saved with a metab_bayes model
#' 
#' @inheritParams get_mcmc
#' @export 
#' @family get_mcmc
get_mcmc.metab_bayes <- function(metab_model) {
  metab_model@mcmc
}

#' Extract any MCMC data list[s] that were stored with the model
#' 
#' A function specific to metab_bayes models. Returns data as formatted to run
#' through the MCMC process or, for nopool models, a list of data lists. These
#' lists are not saved by default; see \code{keep_mcmc_data} argument to
#' \code{\link{specs}} for options.
#' 
#' @param metab_model A Bayesian metabolism model (metab_bayes) from which to 
#'   return the data list that was passed to the MCMC
#' @return The MCMC data list
#' @export
get_mcmc_data <- function(metab_model) {
  UseMethod("get_mcmc_data")
}

#' Retrieve any MCMC data list[s] that were saved with a metab_bayes model
#' 
#' @inheritParams get_mcmc_data
#' @export 
#' @family get_mcmc_data
get_mcmc_data.metab_bayes <- function(metab_model) {
  metab_model@mcmc_data
}

#' Return the log file[s] from a model run
#' 
#' If a log file was created during a model run, this function can retrieve it.
#' 
#' @param metab_model A Bayesian metabolism model (metab_bayes) from which to 
#'   return the log file, if available
#' @return The MCMC log file[s] lines
#' @export
get_log <- function(metab_model) {
  UseMethod("get_log")
}

#' Return the log file[s] from a Bayesian MCMC run
#' 
#' If a log file was created during the MCMC run, metab_bayes() attempted to 
#' capture it. Retrieve what was captured with this function.
#' 
#' @inheritParams get_log
#' @export
#' @family get_log
get_log.metab_bayes <- function(metab_model) {
  out <- metab_model@log
  if(!is.null(out) && length(out) > 0) {
    if(is.list(out)) out <- lapply(out, function(o) { 
      if(!is.null(o)) class(o) <- c('logs_metab', class(o)); o 
    } )
    class(out) <- c('logs_metab', class(out))
  } else {
    message('no log file[s] found')
  }
  out
}

#' Print metab logs
#' 
#' Print metab model compilation and/or fitting logs
#' @inheritParams base::print
#' @export
print.logs_metab <- function(x, ...) {
  if(is.list(x)) {
    outcat <- do.call(c, lapply(names(x), function(llname) { 
      c(paste0("### ", llname, " ###"), '', x[[llname]], '')
    }))
  } else {
    outcat <- x
  }
  cat(paste(outcat, collapse='\n'))
  invisible(x)
}

#' Make metabolism predictions from a fitted metab_model.
#' 
#' Makes daily predictions of GPP, ER, and K600.
#' 
#' @inheritParams predict_metab
#' @return A data.frame of predictions, as for the generic 
#'   \code{\link{predict_metab}}.
#' @import dplyr
#' @importFrom unitted get_units u
#' @export
#' @family predict_metab
predict_metab.metab_bayes <- function(metab_model, date_start=NA, date_end=NA, ..., attach.units=FALSE) {
  # with Bayesian models, the daily mean metabolism values of GPP, ER, and D
  # should have been produced during the model fitting
  
  # decide on the column names to pull and their new values. fit.names and metab.names should be parallel
  Var1 <- Var2 <- '.dplyr.var'
  fit.names <- expand.grid(c('50pct','2.5pct','97.5pct'), c('GPP','ER'), stringsAsFactors=FALSE) %>% #,'D'
    select(Var2, Var1) %>% # variables were in their expand.grid order; now reshuffle them into their paste order
    apply(MARGIN = 1, FUN=function(row) do.call(paste, c(as.list(row), list(sep='_daily_'))))
  metab.names <- expand.grid(c('','.lower','.upper'), c('GPP','ER'), stringsAsFactors=FALSE) %>% #,'D'
    select(Var2, Var1) %>% # variables were in their expand.grid order; now reshuffle them into their paste order
    apply(MARGIN = 1, FUN=function(row) do.call(paste0, as.list(row)))
  
  # pull and retrieve the columns
  fit <- metab_model@fit$daily %>%
    mm_filter_dates(date_start=date_start, date_end=date_end)
  preds <- fit[c('date', fit.names)] %>% 
    setNames(c('date', metab.names)) # these errors & warnings will mostly be date validity notes, unless split_dates==T
  
  # add date-specific fitting warnings and errors as msgs.fit. though these 
  # could also be prediction messages if split_dates==T, we're planning to force
  # split_dates to always be F in the near future. and whenever split_dates==F, 
  # date-specific messages are all just date validity notes and belong in 
  # fitting alone. general messages apply mostly to fitting so are noted here. 
  # get_params also handles general messages, but because we don't call
  # get_params from this predict_metab function, we need to add those messages
  # separately here
  warnings <- errors <- '.dplyr.var'
  if(!is.null(fit) && all(exists(c('date','warnings','errors'), fit))) {
    messages <- fit %>%
      select(date, warnings, errors) %>%
      compress_msgs('msgs.fit', warnings.overall=metab_model@fit$warnings, errors.overall=metab_model@fit$errors)
    preds <- full_join(preds, messages, by='date', copy=TRUE)
  } else {
    preds <- mutate(preds, msgs.fit=NA)
  }
  
  # add general fitting warnings and errors. almost always, general errors 
  # during fitting prohibit prediction and general warnings don't affect 
  # prediction; treat them here as if this is always the case (because 
  # prediction-specific errors or warnings would probably be due to a poorly 
  # written model, which I hope we'll have few of, and I don't know how I'd
  # distinguish since both types of messages come out of Stan)
  preds <- mutate(
    preds,
    warnings=if(length(metab_model@fit$errors) > 0) NA else '',
    errors=if(length(metab_model@fit$errors) > 0) NA else '')
  
  # attach.units if requested
  if(attach.units) {
    pred.units <- get_units(mm_data())[sapply(names(preds), function(x) strsplit(x, '\\.')[[1]][1], USE.NAMES=FALSE)]
    preds <- u(preds, pred.units)
  }
  preds
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
get_params.metab_bayes <- function(metab_model, date_start=NA, date_end=NA, uncertainty='ci', messages=TRUE, ..., attach.units=FALSE) {
  # Stan prohibits '.' in variable names, so we have to convert back from '_' to
  # '.' here to become consistent with the non-Bayesian models
  parnames <- setNames(gsub('_', '\\.', metab_model@specs$params_out), metab_model@specs$params_out)
  for(i in seq_along(parnames)) {
    names(metab_model@fit$daily) <- gsub(names(parnames[i]), parnames[[i]], names(metab_model@fit$daily))
  }
  names(metab_model@fit$daily) <- gsub('_mean$', '', names(metab_model@fit$daily))
  names(metab_model@fit$daily) <- gsub('_sd$', '.sd', names(metab_model@fit$daily))
  # code duplicated in get_params.metab_Kmodel:
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
