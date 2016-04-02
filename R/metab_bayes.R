#' @include metab_model-class.R
NULL

#' Basic Bayesian metabolism model fitting function
#' 
#' Fits a Bayesian model to estimate GPP and ER from input data on DO, 
#' temperature, light, etc. See \code{\link{mm_name}} to choose a Bayesian model
#' and \code{\link{specs}} for relevant options for the \code{specs} 
#' argument.
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
#' dat <- data_metab() # 1 day of example data
#' # fast-ish model version, but still too slow to auto-run in examples
#' mm <- metab_bayes(
#'   specs(mm_name('bayes', err_proc_iid=FALSE, engine='jags'), 
#'         n_chains=1, burnin_steps=300, saved_steps=100),
#'   data=dat)
#' predict_metab(mm)
#' get_fitting_time(mm)
#' plot_DO_preds(predict_DO(mm))
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
        data.frame(discharge.daily = if(ply_validity) mean(data_ply$discharge) else NA)
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
          stop("need ggplot2 for K600_pool='binned' and character value for K600_daily_beta_cuts. ",
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
      # specs$K600_daily_beta_bins[dat_list[['data_daily']]$discharge.bin.daily] or, later,
      # get_specs(fit)$K600_daily_beta_bins[get_data_daily(fit)$discharge.bin.daily]
    }
    
    # Use de-unitted version until we pack up the model to return
    data <- v(dat_list$data)
    data_daily <- v(dat_list$data_daily)
        
    # Check and parse model file path. First try the streamMetabolizer models
    # dir, then try a regular path, then complain / continue depending on
    # whether we found a file. Add the complete path to specs. This is
    # best done here, in the metab_bayes call, so that models can be defined,
    # passed to another computer, and still run successfully.
    specs$model_path <- system.file(paste0("models/", specs$model_name), package="streamMetabolizer")
    if(!file.exists(specs$model_path)) 
      specs$model_path <- specs$model_name
    if(!file.exists(specs$model_path)) 
      stop("could not locate the model file at ", specs$model_path)
    
    # check the format of keep_mcmcs (more checks, below, are split_dates-specific)
    if(is.logical(specs$keep_mcmcs)) {
      if(length(specs$keep_mcmcs) != 1) {
        stop("if keep_mcmcs is logical, it must have length 1")
      }
    } else if(specs$split_dates == FALSE) {
      stop("if split_dates==FALSE, keep_mcmcs must be a single logical value")
    }
    
    # model the data. create outputs bayes_all (a data.frame) and bayes_mcmc (an MCMC object from JAGS or Stan)
    if(specs$split_dates == TRUE) {
      if(!is.logical(specs$keep_mcmcs)) {
        specs$keep_mcmcs <- as.Date(specs$keep_mcmcs)
      }
      # one day at a time, splitting into overlapping 31.5-hr 'plys' for each date
      bayes_daily <- mm_model_by_ply(
        bayes_1ply, data=data, data_daily=data_daily, # for mm_model_by_ply
        day_start=specs$day_start, day_end=specs$day_end, day_tests=specs$day_tests, # for mm_model_by_ply
        specs=specs) # for bayes_1ply
      # if we saved the modeling object[s] in the df, pull them out now
      if('mcmcfit' %in% names(bayes_daily)) {
        bayes_mcmc <- bayes_daily$mcmcfit
        names(bayes_mcmc) <- bayes_daily$date
        bayes_daily$mcmcfit <- NULL
      } else {
        bayes_mcmc <- NULL
      }
      bayes_all <- list(daily=bayes_daily)
      
    } else if(specs$split_dates == FALSE) {
      # all days at a time, after first filtering out bad days
      filtered <- mm_filter_valid_days(
        data, data_daily, day_start=specs$day_start, day_end=specs$day_end, day_tests=specs$day_tests)
      if(length(unique(filtered$data$date)) > 1 && (specs$day_end - specs$day_start) > 24) 
        warning("multi-day models should probably have day_end - day_start <= 24 hours")
      bayes_all_list <- bayes_allply(
        data_all=filtered$data, data_daily_all=filtered$data_daily, removed=filtered$removed,
        specs=specs)
      # if we saved the modeling object, pull it out now
      bayes_mcmc <- bayes_all_list$mcmcfit
      bayes_all_list$mcmcfit <- NULL
      bayes_all <- bayes_all_list # list of dfs
    }
  
  })
  
  # Package and return results
  mm <- metab_model(
    model_class="metab_bayes", 
    info=info,
    fit=bayes_all,
    mcmc=bayes_mcmc,
    fitting_time=fitting_time,
    specs=specs,
    data=dat_list$data, # keep the units if given
    data_daily=dat_list$data_daily)
  
  # Update data with DO predictions
  mm@data <- predict_DO(mm)
  
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
  
  # Calculate metabolism by Bayesian MCMC
  if(length(stop_strs) == 0) {
    bayes_1day <- withCallingHandlers(
      tryCatch({
        # first: try to run the bayes fitting function
        data_list <- prepdata_bayes(
          data=data_ply, data_daily=data_daily_ply, ply_date=ply_date,
          specs=specs, engine=specs$engine, model_name=specs$model_name, priors=specs$priors)
        specs$keep_mcmc <- if(is.logical(specs$keep_mcmcs)) {
          isTRUE(specs$keep_mcmcs)
        } else {
          isTRUE(ply_date %in% specs$keep_mcmcs)
        }
        all_mcmc_args <- c('engine','model_path','params_out','split_dates','keep_mcmc','n_chains','n_cores','adapt_steps','burnin_steps','saved_steps','thin_steps','verbose')
        do.call(mcmc_bayes, c(
          list(data_list=data_list),
          specs[all_mcmc_args[all_mcmc_args %in% names(specs)]]))
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
    bayes_1day <- data.frame(GPP_daily_mean=NA)
  }
  
  # Return, reporting any results, warnings, and errors
  data.frame(bayes_1day,
             warnings=paste0(unique(warn_strs), collapse="; "), 
             errors=paste0(unique(stop_strs), collapse="; "),
             stringsAsFactors=FALSE)
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
  bayes_allday <- withCallingHandlers(
    tryCatch({
      # first: try to run the bayes fitting function
      data_list <- prepdata_bayes(
        data=data_all, data_daily=data_daily_all, ply_date=NA,
        specs=specs, engine=specs$engine, model_name=specs$model_name, priors=specs$priors)
      all_mcmc_args <- c('engine','model_path','params_out','split_dates','keep_mcmc','n_chains','n_cores','adapt_steps','burnin_steps','saved_steps','thin_steps','verbose')
      do.call(mcmc_bayes, c(
        list(data_list=data_list),
        specs[all_mcmc_args[all_mcmc_args %in% names(specs)]]))
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
  if(length(stop_strs) > 0) {
    bayes_allday <- list(daily=data.frame(date=date_vec, GPP_daily_mean=NA))
  } else {
    # check this now so 
    if(length(date_vec) != data_list$d || length(date_vec) != nrow(bayes_allday$daily))
      stop_strs <- c(stop_strs, "couldn't match dates to date indices")
    index <- '.dplyr.var'
    bayes_allday$daily <- bayes_allday$daily %>%
      rename(date=index) %>% 
      mutate(date=date_vec)
  }
  
  # Return, reporting any results, warnings, and errors
  c(bayes_allday,
    list(warnings=unique(warn_strs),
         errors=unique(stop_strs)))
}


#### helpers to the helper ####

#' Prepare data for passing to JAGS or Stan
#' 
#' This function accepts exactly one day's worth of data, (one ply, which might 
#' be 24 hrs or 31.5 or so), which should already be validated. It prepares the 
#' data needed to run a Bayesian MCMC method to estimate GPP, ER, and K600.
#' 
#' @inheritParams mm_model_by_ply_prototype
#' @inheritParams metab
#' @inheritParams specs
#' @param priors logical. Should the data list be modified such that JAGS will 
#'   return priors rather than posteriors?
#' @return list of data for input to runjags_bayes or runstan_bayes
#' @importFrom unitted v
#' @keywords internal
prepdata_bayes <- function(
  data, data_daily, ply_date=NA, # inheritParams mm_model_by_ply_prototype
  specs, # inheritParams metab (for hierarchical priors)
  engine, model_name, #inheritParams specs
  priors=FALSE # inherited by specs
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
  if(length(num_daily_obs) > 1) stop("dates have differing numbers of rows; observations cannot be combined in matrix")
  time_by_date_matrix <- function(vec) {
    switch(
      engine,
      jags=matrix(data=vec, ncol=num_daily_obs, nrow=num_dates, byrow=TRUE),
      stan=matrix(data=vec, nrow=num_daily_obs, ncol=num_dates, byrow=FALSE)
    )
  }
  date_margin <- switch(engine, jags=1, stan=2)
  
  # double-check that our dates are going to line up with the input dates. this 
  # should be redundant w/ above date_table checks, so just being extra careful
  obs_dates <- time_by_date_matrix(as.character(data$date, "%Y-%m-%d"))
  unique_dates <- apply(obs_dates, MARGIN=date_margin, FUN=function(timevec) unique(timevec))
  if(!all.equal(unique_dates, names(date_table))) stop("couldn't fit given dates into matrix")
  
  # confirm that every day has the same modal timestep and put a value on that timestep
  obs_times <- time_by_date_matrix(as.numeric(data$solar.time - data$solar.time[1], units='days'))
  unique_timesteps <- unique(apply(obs_times, MARGIN=date_margin, FUN=function(timevec) unique(round(diff(timevec), digits=12)))) # 10 digits is 8/1000000 of a second. 14 digits exceeds machine precision for datetimes
  if(length(unique_timesteps) != 1) stop("could not determine a single timestep for all observations")
  timestep_days <- mean(apply(obs_times, MARGIN=date_margin, FUN=function(timevec) mean(diff(timevec))))
  
  # parse model name into features for deciding what data to include
  features <- mm_parse_name(model_name)
  
  # Format the data for JAGS/Stan. Stan disallows period-separated names, so
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
        # normalize light by the sum of light in the first 24 hours of the time window
        in_solar_day <- apply(obs_times, MARGIN=date_margin, FUN=function(timevec) {timevec - timevec[1] <= 1} )
        if(engine == 'jags') in_solar_day <- t(in_solar_day)
        mat_light <- time_by_date_matrix(data$light)
        sum_by_date <- switch(engine, jags=rowSums, stan=colSums)
        sweep(mat_light, MARGIN=date_margin, STATS=sum_by_date(mat_light*in_solar_day), FUN=`/`)
      },
      frac_ER  = time_by_date_matrix(timestep_days),
      frac_D   = time_by_date_matrix(timestep_days), # the yackulic shortcut models rely on this being constant over time
      KO2_conv = time_by_date_matrix(convert_k600_to_kGAS(k600=1, temperature=data$temp.water, gas="O2")),
      depth    = time_by_date_matrix(data$depth),
      DO_sat   = time_by_date_matrix(data$DO.sat),
      DO_obs   = time_by_date_matrix(data$DO.obs)
    ),
    
    specs[c(
      # Hyperparameters - this section should be identical to the
      # hyperparameters section of specs|bayes
      c('GPP_daily_mu','GPP_daily_sigma','ER_daily_mu','ER_daily_sigma'),
      switch(
        features$pool_K600,
        none=c('K600_daily_mu', 'K600_daily_sigma'),
        normal=c('K600_daily_mu_mu', 'K600_daily_mu_sigma', 'K600_daily_sigma_shape', 'K600_daily_sigma_rate'),
        linear=c('K600_daily_beta_mu', 'K600_daily_beta_sigma', 'K600_daily_sigma_shape', 'K600_daily_sigma_rate'),
        binned=c('K600_daily_beta_mu', 'K600_daily_beta_sigma', 'K600_daily_sigma_shape', 'K600_daily_sigma_rate')),
      if(features$err_obs_iid) c('err_obs_iid_sigma_shape', 'err_obs_iid_sigma_rate'),
      if(features$err_proc_acor) c('err_proc_acor_phi_shape', 'err_proc_acor_phi_rate', 'err_proc_acor_sigma_shape', 'err_proc_acor_sigma_rate'),
      if(features$err_proc_iid) c('err_proc_iid_sigma_shape', 'err_proc_iid_sigma_rate')
    )]
  )
  if(priors) {
    switch(
      engine,
      jags={ data_list <- data_list[-which(names(data_list)=="DO_obs")] },
      stan={ stop("sorry, Stan doesn't allow NAs in data, so priors can't be TRUE") })
  }
  
  data_list
}

#' Run an MCMC simulation on a formatted data ply
#' 
#' @param data_list a formatted list of inputs to the JAGS model
#' @param engine character string indicating which software to use
#' @param model_path the JAGS model file to use, as a full file path
#' @param params_out a character vector of parameters whose values in the MCMC 
#'   runs should be recorded and summarized
#' @param keep_mcmc logical. If TRUE, the Jags or Stan output object will be 
#'   saved. Be careful; these can be big, and a run with many models might 
#'   overwhelm R's memory.
#' @param n_chains the number of chains to run
#' @param n_cores the number of cores to apply to this run
#' @param adapt_steps the number of steps per chain to use in adapting the model
#' @param burnin_steps the number of steps per chain to run and ignore before
#'   starting to collect MCMC 'data'
#' @param saved_steps the number of MCMC steps per chain to save
#' @param thin_steps the number of steps to move before saving another step. 1
#'   means save all steps.
#' @param verbose logical. give status messages?
#' @return a data.frame of outputs
#' @import parallel
#' @keywords internal
mcmc_bayes <- function(data_list, engine=c('stan','jags'), model_path, params_out, split_dates, keep_mcmc=FALSE, n_chains=4, n_cores=4, adapt_steps=1000, burnin_steps=4000, saved_steps=40000, thin_steps=1, verbose=FALSE) {
  engine <- match.arg(engine)
  bayes_function <- switch(engine, jags = runjags_bayes, stan = runstan_bayes)
  
  tot_cores <- detectCores()
  if (!is.finite(tot_cores)) { tot_cores <- 1 } 
  n_cores <- min(tot_cores, n_cores)
  message(paste0("MCMC (",engine,"): requesting ",n_chains," chains on ",n_cores," of ",tot_cores," available cores\n"))
  
  bayes_function(
    data_list=data_list, model_path=model_path, params_out=params_out, split_dates=split_dates, keep_mcmc=keep_mcmc, n_chains=n_chains, n_cores=n_cores, 
    adapt_steps=adapt_steps, burnin_steps=burnin_steps, saved_steps=saved_steps, thin_steps=thin_steps, verbose=verbose)
}

#' Run JAGS on a formatted data ply
#' 
#' Seems to need to import rjags but does not, for now, because I can't get 
#' rjags to install on the Condor cluster. Including an import rjags line here 
#' allowed runjags to do its job last time I tried.
#' 
#' @inheritParams mcmc_bayes
#' @param ... args passed to other runxx_bayes functions but ignored here
#' @import dplyr
#' @keywords internal
runjags_bayes <- function(data_list, model_path, params_out, split_dates, keep_mcmc=FALSE, n_chains=4, adapt_steps=1000, burnin_steps=4000, saved_steps=40000, thin_steps=1, verbose=FALSE, ...) {
  
  if(!requireNamespace("runjags", quietly = TRUE)) {
    stop("the runjags package is required for JAGS MCMC models")
  }
  
  inits_fun <- function(chain) {
    list(.RNG.name=
           c("base::Wichmann-Hill",
             "base::Marsaglia-Multicarry",
             "base::Super-Duper",
             "base::Mersenne-Twister")[chain])
    # Let JAGS initialize other parameters automatically
  }
  
  runjags_out <- runjags::run.jags(
    method=c("rjags","parallel","snow")[2],
    model=model_path,
    monitor=params_out,
    data=data_list,
    inits=inits_fun,
    n.chains=n_chains,
    adapt=adapt_steps,
    burnin=burnin_steps,
    sample=saved_steps,
    thin=thin_steps,
    summarise=TRUE,
    plots=FALSE,
    silent.jags=!verbose)
  
  # format output
  if(split_dates) {
    # for one-day models, use a 1-row data.frame. see ls('package:rstan')
    jags_mat <- cbind(runjags_out$summary$statistics[,c('Naive SE','Time-series SE')], 
                      runjags_out$summaries,
                      runjags_out$summary$quantiles) %>% as.matrix() # combine 2 matrices of statistics
    names_params <- rep(rownames(jags_mat), each=ncol(jags_mat)) # the GPP, ER, etc. part of the name
    names_stats <- rep(tolower(gsub(" |-", "_", gsub("%", "pct", colnames(jags_mat)))), times=nrow(jags_mat)) # add the mean, sd, etc. part of the name
    jags_out <- jags_mat %>% t %>% c %>% # get a 1D vector of GPP_daily_mean, GPP_sd, ..., ER_daily_mean, ER_daily_sd, ... etc
      t %>% as.data.frame() %>% # convert from 1D vector to 1-row data.frame
      setNames(paste0(names_params, "_", names_stats))
    
    # add the model object if requested
    if(keep_mcmc == TRUE) {
      jags_out <- mutate(jags_out, mcmcfit=list(runjags_out))
    }
  } else {
    # for multi-day or unsplit models, format output into a list of data.frames,
    # one per unique number of nodes sharing a variable name
    jags_mat <- cbind(runjags_out$summary$statistics[,c('Naive SE','Time-series SE')], 
                      runjags_out$summaries,
                      runjags_out$summary$quantiles) %>% as.matrix() # combine 2 matrices of statistics
    colnames(jags_mat) <- gsub("%", "pct", colnames(jags_mat))
    
    stat <- val <- . <- rowname <- variable <- index <- varstat <- '.dplyr_var'
    
    # determine how many unique nrows, & therefore data.frames, there should be
    var_table <- table(gsub("\\[[[:digit:]]\\]", "", rownames(jags_mat)))
    all_dims <- unique(var_table) %>% setNames(., .)
    names(all_dims)[all_dims == 1] <- "overall"
    names(all_dims)[all_dims == data_list$d] <- "daily" # overrides 'overall' if d==1. that's OK
    
    # for each unique nrows, create the data.frame with vars in columns and indices in rows
    jags_out <- lapply(all_dims, function(odim) {
      dim_params <- names(var_table[which(var_table == odim)])
      dim_rows <- sort(do.call(c, lapply(dim_params, function(dp) grep(paste0("^", dp, "(\\[|$)"), rownames(jags_mat)))))
      row_order <- names(sort(sapply(dim_params, function(dp) grep(paste0("^", dp, "(\\[|$)"), rownames(jags_mat))[1])))
      varstat_order <- paste0(rep(row_order, each=ncol(jags_mat)), '_', rep(colnames(jags_mat), times=length(row_order)))
      
      as.data.frame(jags_mat[dim_rows,]) %>%
        add_rownames() %>%
        gather(stat, value=val, 2:ncol(.)) %>%
        mutate(variable=gsub("\\[[[:digit:]]\\]", "", rowname),
               index=if(odim == 1) 1 else sapply(strsplit(rowname, "\\[|\\]"), `[[`, 2),
               varstat=ordered(paste0(variable, '_', stat), varstat_order)) %>%
        select(index, varstat, val) %>%
        spread(varstat, val)
    })
    
    jags_out <- c(jags_out, list(mcmcfit=runjags_out))
  }

  return(jags_out)
}

#' Run Stan on a formatted data ply
#' 
#' @inheritParams mcmc_bayes
#' @param ... args passed to other runxx_bayes functions but ignored here
#' @import parallel
#' @import dplyr
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
  
  runstan_out <- rstan::stan(
    file=model_path,
    data=data_list,
    pars=params_out,
    include=TRUE,
    chains=n_chains,
    warmup=burnin_steps,
    iter=saved_steps+burnin_steps,
    thin=thin_steps,
    init="random",
    save_dso=TRUE, # must be true if you're using more than one core
    verbose=verbose,
    open_progress=FALSE,
    cores=n_cores)

  # this is a good place for a breakpoint when running small numbers of models
  # manually (or keep_mcmc also helps with inspection)
  #   show(runstan_out)
  #   rstan::plot(runstan_out)
  #   pairs(runstan_out)
  #   traceplot(runstan_out)
  
  # format output
  if(split_dates) {
    # for one-day models, use a 1-row data.frame. see ls('package:rstan')
    stan_mat <- rstan::summary(runstan_out)$summary
    names_params <- rep(gsub("\\[1\\]", "", rownames(stan_mat)), each=ncol(stan_mat)) # the GPP, ER, etc. part of the name
    names_stats <- rep(gsub("%", "pct", colnames(stan_mat)), times=nrow(stan_mat)) # add the mean, sd, etc. part of the name
    stan_out <- stan_mat %>% t %>% c %>% # get a 1D vector of GPP_daily_mean, GPP_sd, ..., ER_daily_mean, ER_daily_sd, ... etc
      t %>% as.data.frame() %>% # convert from 1D vector to 1-row data.frame
      setNames(paste0(names_params, "_", names_stats))
  
    # add the model object if requested
    if(keep_mcmc == TRUE) {
      stan_out <- mutate(stan_out, mcmcfit=list(runstan_out))
    }
  } else {
    # for multi-day or unsplit models, format output into a list of data.frames,
    # one per unique number of nodes sharing a variable name
    stan_mat <- rstan::summary(runstan_out)$summary
    colnames(stan_mat) <- gsub("%", "pct", colnames(stan_mat))
    
    stat <- val <- . <- rowname <- variable <- index <- varstat <- '.dplyr_var'
    
    # determine how many unique nrows, & therefore data.frames, there should be
    var_table <- table(gsub("\\[[[:digit:]]+\\]", "", rownames(stan_mat)))
    all_dims <- unique(var_table) %>% setNames(., .)
    names(all_dims)[all_dims == 1] <- "overall"
    names(all_dims)[all_dims == data_list$d] <- "daily" # overrides 'overall' if d==1. that's OK

    # for each unique nrows, create the data.frame with vars in columns and indices in rows
    stan_out <- lapply(all_dims, function(odim) {
      dim_params <- names(var_table[which(var_table == odim)])
      dim_rows <- sort(do.call(c, lapply(dim_params, function(dp) grep(paste0("^", dp, "(\\[|$)"), rownames(stan_mat)))))
      row_order <- names(sort(sapply(dim_params, function(dp) grep(paste0("^", dp, "(\\[|$)"), rownames(stan_mat))[1])))
      varstat_order <- paste0(rep(row_order, each=ncol(stan_mat)), '_', rep(colnames(stan_mat), times=length(row_order)))
      
      as.data.frame(stan_mat[dim_rows,]) %>%
        add_rownames() %>%
        gather(stat, value=val, 2:ncol(.)) %>%
        mutate(variable=gsub("\\[[[:digit:]]+\\]", "", rowname),
               index=if(odim == 1) 1 else as.numeric(sapply(strsplit(rowname, "\\[|\\]"), `[[`, 2)),
               varstat=ordered(paste0(variable, '_', stat), varstat_order)) %>%
        select(index, varstat, val) %>%
        spread(varstat, val)
    })

    stan_out <- c(stan_out, list(mcmcfit=runstan_out))
  } 
  
  return(stan_out)
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
  slots=c(mcmc="ANY")
)

#' Extract any MCMC model objects that were stored with the model
#' 
#' A function specific to metab_bayes models. Returns an MCMC object or, for 
#' nopool models, a list of MCMC objects. These objects are not saved by 
#' default; see \code{keep_mcmcs} argument to \code{\link{specs}} for options.
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

#' Make metabolism predictions from a fitted metab_model.
#' 
#' Makes daily predictions of GPP, ER, and K600.
#' 
#' @inheritParams predict_metab
#' @return A data.frame of predictions, as for the generic 
#'   \code{\link{predict_metab}}.
#' @import dplyr
#' @export
#' @family predict_metab
predict_metab.metab_bayes <- function(metab_model, date_start=NA, date_end=NA, ...) {
  warn_strs <- metab_model@fit$warnings
  stop_strs <- metab_model@fit$errors
  metab_model <- metab_model(fit=mutate(
    metab_model@fit$daily, 
    warnings=if(length(warn_strs) == 0) NA else "see 'warnings' attribute", 
    errors=if(length(stop_strs) == 0) NA else "see 'errors' attribute"))
  preds <- NextMethod()
  attr(preds, 'warnings') <- warn_strs
  attr(preds, 'errors') <- stop_strs
  preds
}