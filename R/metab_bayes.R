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
#'
#' # view the Stan model file as stored on your system
#' file.edit(get_specs(mm)$model_path)
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
    pool_K600_type <- mm_parse_name(specs$model_name, expand = TRUE)$pool_K600_type
    if(xor(num_discharge_cols > 0, pool_K600_type %in% c('linear','binned')))
      stop('discharge data should be included if & only if pool_K600_type indicates hierarchy')
    if(num_discharge_cols > 1)
      stop('either discharge or discharge.daily may be specified, but not both')

    # Handle discharge. If K600 is a hierarchical function of discharge and
    # data$discharge was given, compute daily discharge and store in data.daily
    # where it'll be accessible for the user to inspect it after model fitting
    if((pool_K600_type %in% c('linear','binned')) && ('discharge' %in% names(dat_list$data))) {
      # calculate daily discharge
      dailymean <- function(data_ply, data_daily_ply, day_start, day_end, ply_date, ply_validity, timestep_days, ...) {
        data.frame(discharge.daily = if(isTRUE(ply_validity[1])) mean(data_ply$discharge) else NA)
      }
      dischdaily <- mm_model_by_ply(
        model_fun=dailymean, data=v(dat_list$data), day_start=specs$day_start, day_end=specs$day_end,
        day_tests=specs$day_tests, required_timestep=specs$required_timestep)

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
      dat_list$data_daily$lnQ.daily <- log(v(dat_list$data_daily$discharge.daily))
    }
    # If we need discharge bins, compute & store those now, as well
    if(pool_K600_type %in% c('binned')) {
      # linear interpolation from node to node, horizontal at the edges
      bounds <- c(-Inf, specs$K600_lnQ_nodes_centers, Inf)
      cuts <- cut(dat_list$data_daily$lnQ.daily, breaks=bounds, ordered_result=TRUE)
      widths <- diff(bounds)[cuts]
      bins <- rbind(pmax(1, as.numeric(cuts) - 1), pmin(length(specs$K600_lnQ_nodes_centers), as.numeric(cuts)))
      weights <- ifelse(is.infinite(widths), 1, (bounds[as.numeric(cuts)+1] - dat_list$data_daily$lnQ.daily)/widths)

      # package info so it gets passed to specs
      dat_list$data_daily <- mutate(
        dat_list$data_daily,
        lnQ.bin1 = bins[1,],
        lnQ.bin2 = bins[2,],
        lnQ.bin1.weight = weights,
        lnQ.bin2.weight = 1-weights)
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

    # model the data. create outputs bayes_all (a data.frame) and bayes_mcmc (an
    # MCMC object from Stan)
    if(specs$split_dates == TRUE) {
      if(!is.logical(specs$keep_mcmcs)) specs$keep_mcmcs <- as.Date(specs$keep_mcmcs)
      if(!is.logical(specs$keep_mcmc_data)) specs$keep_mcmc_data <- as.Date(specs$keep_mcmc_data)
      # one day at a time, splitting into overlapping ~24-hr 'plys' for each date
      bayes_daily <- mm_model_by_ply(
        bayes_1ply, data=data, data_daily=data_daily, # for mm_model_by_ply
        day_start=specs$day_start, day_end=specs$day_end, day_tests=specs$day_tests, required_timestep=specs$required_timestep, # for mm_model_by_ply
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
      bayes_compile_time <- bayes_all_list[['compile_time']]
      bayes_mcmc <- extract_object_list('mcmcfit')
      bayes_mcmc_data <- extract_object_list('mcmc_data')
      bayes_all <- list(daily=bayes_daily)
      if(nrow(bayes_all$daily) == 0 || length(which(bayes_daily$valid_day)) == 0)
        bayes_all$errors <- c(bayes_all$errors, "no valid days of data")

    } else if(specs$split_dates == FALSE) {
      # all days at a time, after first filtering out bad days
      filtered <- mm_filter_valid_days(
        data, data_daily, day_start=specs$day_start, day_end=specs$day_end,
        day_tests=specs$day_tests, required_timestep=specs$required_timestep)
      if(length(unique(filtered$data$date)) > 1 && (specs$day_end - specs$day_start) > 24)
        warning("multi-day models should probably have day_end - day_start <= 24 hours")
      bayes_all_list <- bayes_allply(
        data_all=filtered$data, data_daily_all=filtered$data_daily, removed=filtered$removed,
        specs=specs)
      # if we saved the modeling object[s] in the list, pull them out now
      . <- '.dplyr.var'
      bayes_log <- bayes_all_list[c('compile_log', 'log')] %>%
        setNames(c('Compilation','MCMC_All_Days')) %>% { .[!sapply(., is.null)] }
      bayes_compile_time <- bayes_all_list[['compile_time']]
      bayes_mcmc <- bayes_all_list$mcmcfit
      bayes_mcmc_data <- bayes_all_list$mcmc_data
      # now a list of dfs, log, warnings, and errors
      bayes_all <- bayes_all_list[!(names(bayes_all_list) %in% c('compile_log','compile_time','log','mcmcfit','mcmc_data'))]
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
    fitting_time=fitting_time - bayes_compile_time,
    compile_time=bayes_compile_time,
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
      warntxt <- paste0(
        'Modeling failed\n',
        if(length(bayes_all$warnings) > 0) paste0('  Warnings:\n', paste0('    ', bayes_all$warnings, collapse='\n')),
        if(length(bayes_all$errors) > 0) paste0('  Errors:\n', paste0('    ', bayes_all$errors, collapse='\n')))
      warning(warntxt)
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
          data=data_ply, data_daily=data_daily_ply, ply_date=ply_date, specs=specs)
        do.call(runstan_bayes, c(list(data_list=data_list), specs))
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

  # stop_strs may have accumulated during prepdata_bayes() or runstan_bayes()
  # calls. If failed, use dummy data to fill in the model output with NAs.
  if(length(stop_strs) > 0) {
    bayes_1day <- data.frame(
      GPP_daily_2.5pct=NA, GPP_daily_50pct=NA, GPP_daily_97.5pct=NA,
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
  data_list <- list(d=1, n=1) # (in case it doesn't get assigned in the tryCatch)
  bayes_allday <- withCallingHandlers(
    tryCatch({
      if(is.null(data_all) || nrow(data_all) == 0) stop("no valid days of data")
      # first: try to run the bayes fitting function
      data_list <- prepdata_bayes(
        data=data_all, data_daily=data_daily_all, ply_date=NA, specs=specs)
      specs$keep_mcmc <- specs$keep_mcmcs
      do.call(runstan_bayes, c(list(data_list=data_list), specs))
    }, error=function(err) {
      # on error: give up, remembering error. dummy values provided below
      stop_strs <<- c(stop_strs, err$message)
      NA
    }), warning=function(war) {
      # on warning: record the warning and run again
      warn_strs <<- c(warn_strs, war$message)
      invokeRestart("muffleWarning")
    })

  # match date and time info to indices
  date_df <- tibble::tibble(
    date=as.Date(unique(data_all$date)),
    date_index=seq_len(data_list$d))
  datetime_df <- tibble::tibble(
    solar.time=data_all$solar.time,
    date_index=rep(seq_len(data_list$d), each=data_list$n),
    time_index=rep(seq_len(data_list$n), times=data_list$d)) %>%
    left_join(date_df, by='date_index')

  # stop_strs may have accumulated during prepdata_bayes() or runstan_bayes()
  # calls. If failed, use dummy data to fill in the model output with NAs.
  if(length(stop_strs) > 0 || any(grepl("^Stan model .* does not contain samples", warn_strs))) {
    na_vec <- rep(as.numeric(NA), nrow(date_df))
    bayes_allday <- c(
      list(daily=data.frame(
        date=date_df$date, GPP_daily_2.5pct=na_vec, GPP_daily_50pct=na_vec, GPP_daily_97.5pct=na_vec,
        ER_daily_2.5pct=na_vec, ER_daily_50pct=na_vec, ER_daily_97.5pct=na_vec,
        K600_daily_2.5pct=na_vec, K600_daily_50pct=na_vec, K600_daily_97.5pct=na_vec)),
      list(log=if(exists('bayes_allday') && is.list(bayes_allday)) {
        bayes_allday[c('compile_log', 'log')]
      } else NULL ))
  } else {
    # match dates back to daily estimates, datetimes back to inst
    date_index <- time_index <- index <- '.dplyr.var'
    bayes_allday$daily <- bayes_allday$daily %>%
      left_join(date_df, by='date_index') %>%
      select(-date_index, -time_index, -index) %>%
      select(date, everything())
    if(!is.null(bayes_allday$inst)) {
      bayes_allday$inst <- bayes_allday$inst %>%
        left_join(datetime_df, by=c('date_index','time_index')) %>%
        select(-date_index, -time_index, -index) %>%
        select(date, solar.time, everything())
    }
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
#' This function accepts pre-validated data (though more problems may be
#' discovered here). It prepares the data needed to run a Bayesian MCMC method
#' to estimate GPP, ER, and K600.
#'
#' @inheritParams mm_model_by_ply_prototype
#' @inheritParams metab
#' @return list of data for input to runstan_bayes
#' @importFrom unitted v
#' @keywords internal
prepdata_bayes <- function(
  data, data_daily, ply_date=NA, # inheritParams mm_model_by_ply_prototype
  specs # inheritParams metab (for hierarchical priors, model_name)
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
  obs_dates <- time_by_date_matrix(format(data$date, format="%Y-%m-%d"))
  unique_dates <- apply(obs_dates, MARGIN=2, FUN=function(timevec) unique(timevec))
  if(!all.equal(unique_dates, names(date_table))) stop("couldn't fit given dates into matrix")

  # confirm that every day has the same modal timestep and put a value on that
  # timestep. the tolerance for uniqueness within each day is set by the default
  # for mm_get_timestep. the tolerance for uniqueness across days is 10 digits
  # is 8/1000000 of a second. 14 digits exceeds machine precision for datetimes
  obs_times <- time_by_date_matrix(as.numeric(data$solar.time - data$solar.time[1], units='days'))
  timestep_eachday <- apply(obs_times, MARGIN=2, FUN=mm_get_timestep, format='mean', require_unique=TRUE)
  if(length(unique(round(timestep_eachday, digits=10))) != 1) stop("could not determine a single timestep for all observations")
  timestep_days <- mean(timestep_eachday)
  n24 <- round(1/timestep_days)

  # give message if day length is too short
  if(n24 > num_daily_obs) stop("day_end - day_start < 24 hours; aborting because daily metabolism could be wrong")

  # parse model name into features for deciding what data to include
  features <- mm_parse_name(specs$model_name, expand=TRUE)

  # Format the data for Stan. Stan disallows period-separated names, so
  # change all the input data to underscore-separated. parameters given in
  # specs are already underscore-separated for this reason
  data_list = c(
    list(

      # Overall
      d = num_dates,
      timestep = timestep_days, # length of each timestep in days
      n24 = n24, # number of observations in first 24 hours, for computing GPP & ER

      # Daily
      n = num_daily_obs # one value applicable to every day
    ),

    switch(
      features$pool_K600_type,
      linear = list(lnQ_daily = array(time_by_date_matrix(data_daily$lnQ.daily), dim=num_dates)),
      binned = list(
          b = length(specs$K600_lnQ_nodes_centers),
          lnQ_bins = rbind(data_daily$lnQ.bin1, data_daily$lnQ.bin2),
          lnQ_bin_weights = rbind(data_daily$lnQ.bin1.weight, data_daily$lnQ.bin2.weight))
    ),

    list(
      DO_obs_1 = array(time_by_date_matrix(data$DO.obs)[1,], dim=num_dates)), # duplication of effort below should be small compared to MCMC time

    # Every timestep
    switch(
      features$GPP_fun,
      linlight = list(
        # X_mult_Y syntax: X = process reflected by multiplier, Y = quantity
        # modified by multiplier
        light_mult_GPP = {
          mat_light <- time_by_date_matrix(data$light)
          # normalize light by the sum of light in the first 24 hours of the time window
          in_solar_day <- apply(obs_times, MARGIN=2, FUN=function(timevec) {timevec - timevec[1] <= 1} )
          daily_totals <- colSums(mat_light*in_solar_day)
          if(any(daily_totals <= 0)) {
            stop('daily light total is <= 0 on ', paste(names(date_table)[which(daily_totals <= 0)], collapse=', '))
          }
          sweep(mat_light, MARGIN=2, STATS=daily_totals, FUN=`/`) / timestep_days
        }),
      satlight = list(
        light = time_by_date_matrix(data$light)
      )
    ),

    list(
      # X_mult_Y syntax: X = process reflected by multiplier, Y = quantity
      # modified by multiplier
      const_mult_ER  = time_by_date_matrix(1),
      KO2_conv = {
        KO2_conv_vec <- suppressWarnings(convert_k600_to_kGAS(k600=1, temperature=data$temp.water, gas="O2"))
        if(any(is.nan(KO2_conv_vec))) {
          bad_rows <- which(is.nan(KO2_conv_vec))
          show_rows <- bad_rows[seq_len(min(length(bad_rows), 3))]
          more_rows <- if(length(bad_rows) > 3) length(bad_rows) - 3 else NA
          bad_times <- data$solar.time[show_rows]
          bad_temps <- data$temp.water[show_rows]
          stop(sprintf(
            'NaNs in KO2-K600 conversion at %s%s',
            paste0(sprintf('%s (temp.water=%0.2f)', format(bad_times, '%Y-%m-%d %H:%M:%S'), bad_temps), collapse=', '),
            if(!is.na(more_rows)) sprintf(', and %d more rows', more_rows) else ''
          ))
        }
        time_by_date_matrix(KO2_conv_vec)
      },
      depth    = time_by_date_matrix(data$depth),
      DO_sat   = time_by_date_matrix(data$DO.sat),
      DO_obs   = time_by_date_matrix(data$DO.obs)
    ),

    specs[specs$params_in]
  )
  if(features$pool_K600_type == 'binned') {
    data_list$K600_lnQ_nodes_meanlog <- array(data_list$K600_lnQ_nodes_meanlog, dim=data_list$b)
    data_list$K600_lnQ_nodes_sdlog <- array(data_list$K600_lnQ_nodes_sdlog, dim=data_list$b)
  }

  # check that the params_out are unique (non-unique messes up our parsing of
  # the stanfit output)
  if(length(specs$params_out) != length(unique(specs$params_out))) {
    stop('params_out must all be unique')
  }

  data_list
}

#' Run Stan on a formatted data ply
#'
#' @param data_list a formatted list of inputs to the Stan model
#' @param model_path the Stan model file to use, as a full file path
#' @param model_name the coded model name, as from mm_name, giving the model
#'   structure
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
#' @import parallel
#' @import dplyr
#' @import tibble
#' @importFrom tidyr gather spread
#' @keywords internal
runstan_bayes <- function(
  data_list, model_path, model_name, params_out, split_dates, keep_mcmc=FALSE,
  n_chains=4, n_cores=4, burnin_steps=1000, saved_steps=1000, thin_steps=1,
  verbose=FALSE, ...) {

  # determine how many cores to use
  tot_cores <- detectCores()
  if (!is.finite(tot_cores)) { tot_cores <- 1 }
  n_cores <- min(tot_cores, n_cores)
  if(verbose) message(paste0("MCMC (","Stan","): requesting ",n_chains," chains on ",n_cores," of ",tot_cores," available cores"))

  # stan() can't find its own function cpp_object_initializer() unless the
  # namespace is loaded. requireNamespace is somehow not doing this. Thoughts
  # (not solution):
  # https://stat.ethz.ch/pipermail/r-devel/2014-September/069803.html
  if(!suppressPackageStartupMessages(require(rstan))) {
    stop("the rstan package is required for Stan MCMC models")
  }

  # use auto_write=TRUE to recompile if needed, or load from existing .rds file
  # without recompiling if possible
  compile_time <- system.time({})
  mobj_path <- gsub('.stan$', '.stanrds', model_path)
  if(!file.exists(mobj_path) || file.info(mobj_path)$mtime < file.info(model_path)$mtime) {
    if(verbose) message("compiling Stan model")
    compile_time <- system.time({
      compile_log <- capture.output({
        stan_mobj <- rstan::stan_model(file=model_path, auto_write=TRUE)
      }, type=c('output'), split=verbose)
    })
    rm(stan_mobj)
    gc() # this humble line saves us from many horrible R crashes
    autowrite_path <- gsub('.stan$', '.rds', model_path)
    if(!file.exists(autowrite_path)) autowrite_path <- gsub('.stan$', '.rda', model_path) # for backwards compatibility with rstan < 2.13
    if(!file.exists(autowrite_path)) autowrite_path <- file.path(tempdir(), basename(autowrite_path))
    if(!file.exists(autowrite_path)) {
      warning('could not find saved rds model file')
    } else {
      tryCatch({
        file.copy(autowrite_path, mobj_path, overwrite=TRUE)
        file.remove(autowrite_path)
      }, error=function(e) {
        warning('could not copy Stan rds to .stanrds file: ', e$message)
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
    stan_out <- format_mcmc_mat_nosplit(stan_mat, data_list$d, data_list$n, model_name, keep_mcmc, runstan_out)
  }

  # attach the contents of the most recent logfile in tempdir(), which should be for this model
  newlogfiles <- normalizePath(file.path(tempdir(), grep("_StanProgress.txt", dir(tempdir()), value=TRUE)))
  logfile <- setdiff(newlogfiles, oldlogfiles)
  log <- if(length(logfile) > 0) readLines(logfile) else consolelog
  stan_out <- c(stan_out, c(
    list(log=log),
    if(exists('compile_log')) list(compile_log=compile_log),
    list(compile_time=compile_time)))

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
format_mcmc_mat_nosplit <- function(mcmc_mat, data_list_d, data_list_n, model_name, keep_mcmc, runmcmc_out) {

  # assign parameters to appropriately sized data.frames. list every anticipated
  # parameter here, but also catch other parameters below (for custom models, or
  # unusual parameters like err_proc_iid_sigma_scaled). pull the list (manually)
  # from the parameters block of mm_generate_mcmc_file. it's only important to
  # distinguish among model types when the resulting dimensions differ between
  # types and would imply conflicting ideas for nrows of the named data.frame
  features <- mm_parse_name(model_name, expand=TRUE)
  par_homes <- list(
    overall = c(
      'err_obs_iid_sigma', 'err_obs_iid_sigma_scaled',
      'err_proc_iid_sigma', 'err_proc_iid_sigma_scaled',
      'err_proc_acor_phi', 'err_proc_acor_sigma_scaled',
      'lp__'),
    KQ_overall = c( # n=1
      switch(
        features$pool_K600_type,
        none=c(),
        normal=c('K600_daily_predlog'),
        linear=c('lnK600_lnQ_intercept', 'lnK600_lnQ_slope'),
        binned=c()),
      'K600_daily_sdlog', 'K600_daily_sdlog_scaled', 'K600_daily_sigma', 'K600_daily_sigma_scaled'),
    KQ_binned = c( # n=few to many
      'lnK600_lnQ_nodes'),
    daily = c(
      'GPP', 'ER',
      'GPP_daily', 'Pmax', 'alpha', 'ER_daily', 'K600_daily', 'DO_R2',
      if(features$pool_K600_type %in% c('linear','binned')) 'K600_daily_predlog',
      if(features$err_proc_GPP) 'GPP_pseudo_R2'),
    inst = c(
      'DO_mod', 'DO_mod_partial', # d*n
      'DO_mod_partial_sigma', # d*n
      'GPP_inst', 'ER_inst', 'KO2_inst', # d*n
      'GPP_inst_partial', 'err_proc_GPP', # d*n for GPP process error
      'err_proc_acor_inc', 'err_proc_acor', # can be d*n (trapezoid) or d*(n-1) (euler), same timestamp indexing as GPP_inst
      'err_obs_iid', 'err_proc_iid' # d*(n-1), timestamp[i+1] relates to var[i:i+1]
    )
  )

  # declare dplyr variables
  stat <- val <- . <- rowname <- variable <- varstat <-
    indexstr <- date_index <- time_index <- index <- '.dplyr_var'

  # determine which data.frames to create and which params to include in each
  var_table <- table(gsub("\\[[[:digit:]|,]+\\]", "", rownames(mcmc_mat)))
  par_dfs <- sapply(names(var_table), function(parname) {
    home <- names(par_homes)[which(sapply(par_homes, function(vd) parname %in% vd))]
    if(length(home) == 0) home <- var_table[[parname]]
    return(home)
  })
  all_dims <- lapply(setNames(nm=unique(par_dfs)), function(upd) names(par_dfs)[which(par_dfs == upd)])

  # for each unique nrows, create the data.frame with vars in columns and indices in rows
  colnames(mcmc_mat) <- gsub("%", "pct", colnames(mcmc_mat))
  mcmc_out <- lapply(setNames(nm=names(all_dims)), function(dfname) {
    df_params <- all_dims[[dfname]]
    dim_rows <- sort(do.call(c, lapply(df_params, function(dp) grep(paste0("^", dp, "(\\[|$)"), rownames(mcmc_mat)))))
    row_order <- names(sort(sapply(df_params, function(dp) grep(paste0("^", dp, "(\\[|$)"), rownames(mcmc_mat))[1]))) # use the same order as mcmc_mat
    varstat_order <- paste0(rep(row_order, each=ncol(mcmc_mat)), '_', rep(colnames(mcmc_mat), times=length(row_order)))
    par_dims <- sapply(df_params, function(dp) length(grep(paste0("^", dp, "(\\[|$)"), rownames(mcmc_mat))))

    tibble::as_tibble(mcmc_mat[dim_rows,,drop=FALSE]) %>%
      mutate(rowname=rownames(mcmc_mat[dim_rows,,drop=FALSE])) %>%
      select(rowname, everything()) %>%
      gather(stat, value=val, 2:ncol(.)) %>%
      mutate(variable=gsub("\\[[[:digit:]|,]+\\]", "", rowname),
             indexstr=if(1 %in% par_dims) '1' else sapply(strsplit(rowname, "\\[|\\]"), `[[`, 2),
             # parse/factorify the index for ordering
             index=
               if(dfname %in% c('daily','inst') || any(grepl(',', indexstr))) {
                 indexstr
               } else {
                 as.numeric(indexstr)
               },
             date_index=
               if(dfname=='daily') {
                 as.numeric(indexstr)
               } else if(dfname=='inst' || any(grepl(',', indexstr))) {
                 sapply(strsplit(indexstr, ','), function(ind) as.numeric(ind[2]))
               } else {
                 NA
               },
             time_index=
               if(dfname=='inst' || any(grepl(',', indexstr))) {
                 sapply(strsplit(indexstr, ','), function(ind) as.numeric(ind[1]))
               } else {
                 NA
               },
             # determine the order of statistics for each variable (mean, se_mean, etc.)
             varstat=ordered(paste0(variable, '_', stat), varstat_order)) %>%
      select(date_index, time_index, index, varstat, val) %>%
      arrange(date_index, time_index, index) %>%
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
  slots=c(log="ANY", mcmc="ANY", mcmc_data="ANY", compile_time="ANY")
)

#' Extract any MCMC model objects that were stored with the model
#'
#' A function specific to metab_bayes models. Returns an MCMC object of class
#' `stanfit` ([rstan::stanfit-class]), which is saved in the metab_model by
#' default because you should almost always inspect it; see `keep_mcmcs`
#' argument to [specs()] for options for saving space. The \code{rstan} methods
#' for [rstan::stanfit-class] objects include `summary()`, `get_stancode()`,
#' `stan_dens()`, `stan_diag()`, and many more. See
#' `?'rstan-plotting-functions'`, [rstan::stanfit-class] and the
#' \href{https://cran.r-project.org/web/packages/rstan/rstan.pdf}{rstan manual}.
#'
#' @md
#' @param metab_model A Bayesian metabolism model (metab_bayes) from which to
#'   return the MCMC model object[s]
#' @return The MCMC model object[s]
#' @export
get_mcmc <- function(metab_model) {
  UseMethod("get_mcmc")
}

#' @describeIn get_mcmc Get the Bayesian MCMC model object
#' @export
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

#' @describeIn get_mcmc_data Retrieve any MCMC data list[s] that were saved with
#'   a metab_bayes model
#' @export
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

#' @describeIn get_log If a log file was created during the Bayesian MCMC run,
#'   metab_bayes() attempted to capture it. Retrieve what was captured with this
#'   function.
#' @export
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
#' @param x an object to print
#' @param ... ignored; included only for compatibility with `base::print`
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

#' @describeIn predict_metab Pulls daily metabolism estimates out of the Stan
#'   model results; looks for \code{GPP} or \code{GPP_daily} and for \code{ER}
#'   or \code{ER_daily} among the \code{params_out} (see \code{\link{specs}}),
#'   which means you can save just one (or both) of those sets of daily
#'   parameters when running the Stan model. Saving fewer parameters can help
#'   models run faster and use less RAM.
#' @export
#' @import dplyr
#' @importFrom unitted get_units u
#' @importFrom lifecycle deprecated is_present
predict_metab.metab_bayes <- function(metab_model, date_start=NA, date_end=NA, ..., attach.units=deprecated()) {
  # with Bayesian models, the daily mean metabolism values of GPP, ER, and D
  # should have been produced during the model fitting

  # check units-related arguments
  if (lifecycle::is_present(attach.units)) {
    unitted_deprecate_warn("predict_metab(attach.units)")
  } else {
    attach.units <- FALSE
  }

  # decide on the column names to pull and their new values. fit.names and metab.names should be parallel
  Var1 <- Var2 <- '.dplyr.var'
  fit.names.metab <- expand.grid(c('50pct','2.5pct','97.5pct'), c('GPP','ER'), stringsAsFactors=FALSE) %>% #,'D'
    select(Var2, Var1) %>% # variables were in their expand.grid order; now reshuffle them into their paste order
    apply(MARGIN = 1, FUN=function(row) do.call(paste, c(as.list(row), list(sep='_'))))
  fit.names.param <- expand.grid(c('50pct','2.5pct','97.5pct'), c('GPP_daily','ER_daily'), stringsAsFactors=FALSE) %>% #,'D'
    select(Var2, Var1) %>% # variables were in their expand.grid order; now reshuffle them into their paste order
    apply(MARGIN = 1, FUN=function(row) do.call(paste, c(as.list(row), list(sep='_'))))
  metab.names <- expand.grid(c('','.lower','.upper'), c('GPP','ER'), stringsAsFactors=FALSE) %>% #,'D'
    select(Var2, Var1) %>% # variables were in their expand.grid order; now reshuffle them into their paste order
    apply(MARGIN = 1, FUN=function(row) do.call(paste0, as.list(row)))

  # pull and retrieve the columns
  fit <- metab_model@fit$daily %>%
    mm_filter_dates(date_start=date_start, date_end=date_end)
  fit.names <- if(all(fit.names.metab %in% names(fit))) {
    fit.names.metab
  } else if(all(fit.names.param %in% names(fit))) {
    fit.names.param
  } else {
    stop('could find neither GPP & ER nor GPP_daily & ER_daily in the model fit')
  }
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
  if(!is.null(fit) && all(c('date','warnings','errors') %in% names(fit))) {
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

#' @describeIn get_params Does a little formatting to convert from Stan output
#'   to streamMetabolizer parameter names; otherwise the same as
#'   \code{get_params.metab_model}
#' @importFrom lifecycle deprecated is_present
#' @export
get_params.metab_bayes <- function(
  metab_model, date_start=NA, date_end=NA, uncertainty='ci', messages=TRUE,
  ..., attach.units=deprecated()) {

  # check units-related arguments
  if (lifecycle::is_present(attach.units)) {
    unitted_deprecate_warn("get_params(attach.units)")
  } else {
    attach.units <- FALSE
  }

  # Stan prohibits '.' in variable names, so we have to convert back from '_' to
  # '.' here to become consistent with the non-Bayesian models
  parnames <- setNames(gsub('_', '\\.', metab_model@specs$params_out), metab_model@specs$params_out)
  parnames <- parnames[order(nchar(parnames), decreasing=TRUE)]
  for(i in seq_along(parnames)) {
    names(metab_model@fit$daily) <- gsub(names(parnames[i]), parnames[[i]], names(metab_model@fit$daily))
  }
  names(metab_model@fit$daily) <- gsub('_mean$', '', names(metab_model@fit$daily))
  names(metab_model@fit$daily) <- gsub('_sd$', '.sd', names(metab_model@fit$daily))
  names(metab_model@fit$daily) <- gsub('_50pct$', '.median', names(metab_model@fit$daily))
  names(metab_model@fit$daily) <- gsub('_2.5pct$', '.lower', names(metab_model@fit$daily))
  names(metab_model@fit$daily) <- gsub('_97.5pct$', '.upper', names(metab_model@fit$daily))
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
  metab_model@fit <- metab_model@fit$daily # TEMPORARY we're still converting fit$daily to fit until #247, #229
  NextMethod()
}
