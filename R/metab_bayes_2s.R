#' @include metab_model-class.R metab_bayes.R
NULL

# Suppress R CMD CHECK NOTEs for column names used as unbound globals in
# dplyr NSE calls (e.g. mutate, summarise). These are data frame column
# names resolved at runtime, not missing variable declarations.
utils::globalVariables(c(".", "metab_50pct", "DO.mod.down"))

#' Two-station Bayesian metabolism model fitting function
#'
#' Fits a two-station (upstream/downstream, Variable Flow Two-Station) Bayesian
#' model to estimate GPP, ER, and K600 from paired upstream and downstream DO,
#' temperature, light, and travel-time data, using the single fixed Stan model
#' in \code{inst/models/b2_np_oi_tr_plrckm.stan}. See \code{\link{mm_name}} to
#' choose a Bayesian model and \code{\link{specs}} for relevant options for the
#' \code{specs} argument.
#'
#' Unlike \code{\link{metab_bayes}}, which supports many model structures via
#' \code{split_dates}/\code{pool_K600}/etc., \code{metab_bayes_2s} always
#' fits every date jointly in a single Stan call (\code{specs$split_dates} is
#' forced to \code{FALSE} by \code{\link{specs}}), because the
#' upstream-downstream lag shift ties each date's first modeled rows to the
#' previous date's last rows.
#'
#' @inheritParams metab
#' @return A metab_bayes_2s object containing the fitted model. This object
#'   can be inspected with the functions in the
#'   \code{\link{metab_model_interface}} and also \code{\link{get_mcmc}}.
#'
#' @section Two-station data requirements: In addition to the checks
#'   performed by \code{\link{mm_validate_data}}, \code{data$travel.time} (the
#'   reach travel time between the upstream and downstream stations, in days)
#'   must be strictly positive and less than 1 day (values >= 1 usually
#'   indicate travel time was supplied in the wrong units). There must also be
#'   enough lead-in observations of upstream DO before the first row of
#'   \code{data} to cover the longest travel time in the dataset, given the
#'   (median) timestep of \code{data$solar.time}.
#'
#' @export
#' @family metab_model
#' @importFrom utils modifyList
metab_bayes_2s <- function(
  specs=specs(mm_name('bayes_2s')),
  data=mm_data(solar.time, DO.obs.up, DO.sat.up, DO.obs.down, DO.sat.down,
               light, depth, temp.water, travel.time),
  data_daily=mm_data(date, optional='all'),
  info=NULL
) {

  stanfit <- NULL
  fitting_time <- system.time({
    # Check data for correct column names, units, and travel.time bounds
    # (mm_validate_data()), then check lead-in coverage
    # (mm_validate_data_2station()), before any data prep begins
    dat_list <- mm_validate_data(data, data_daily, 'metab_bayes_2s')
    mm_validate_data_2station(dat_list$data)

    data_v <- v(dat_list$data)
    travel_time <- data_v$travel.time
    solar_time <- data_v$solar.time

    # Reconstruct the same "modeled rows" (post-lead-in-trim) index set that
    # prepdata_bayes_2s() computes internally, using the identical
    # timestep_days/lag/max_lag formula, so that Stan's date_index/time_index
    # can be mapped back to actual dates and solar.times for the daily and
    # instantaneous results below. (Duplicated here rather than exposed by
    # prepdata_bayes_2s(), which returns only the Stan-ready matrices.)
    timestep_days <- stats::median(as.numeric(diff(solar_time), units='days'))
    max_lag <- max(round(travel_time / timestep_days))
    n_total <- nrow(dat_list$data)
    keep <- seq.int(max_lag + 1, n_total)
    modeled_solar_time <- solar_time[keep]
    modeled_dates <- as.Date(modeled_solar_time)
    date_df <- tibble::tibble(date=unique(modeled_dates), date_index=seq_along(unique(modeled_dates)))
    n_days <- nrow(date_df)
    if(length(keep) %% n_days != 0) {
      stop(paste0(
        'dates have differing numbers of modeled rows after lead-in removal; ',
        'observations cannot be combined into a matrix: ',
        paste(sprintf('%s (%d rows)', names(table(modeled_dates)), table(modeled_dates)), collapse=', ')))
    }
    n_obs <- length(keep) / n_days
    obs_index_df <- tibble::tibble(
      solar.time=modeled_solar_time,
      DO.obs.down=data_v$DO.obs.down[keep],
      date_index=rep(date_df$date_index, each=n_obs),
      time_index=rep(seq_len(n_obs), times=n_days))

    # Prepare the Stan data list (matrices from data, plus scalar priors from
    # specs). modifyList (not c()) is used because prepdata_bayes_2s() already
    # supplies K600_lnorm_meanlog/K600_lnorm_sdlog (read from specs), and
    # those two names are also in specs$params_in; a plain c() would create
    # duplicate-named list elements instead of overriding
    data_list <- prepdata_bayes_2s(dat_list$data, specs=specs)
    data_list <- modifyList(data_list, specs[specs$params_in])

    # Check and parse model file path
    specs$model_path <- mm_locate_filename(specs$model_name)

    # determine how many cores to use, as in runstan_bayes()
    tot_cores <- parallel::detectCores()
    if(!is.finite(tot_cores)) tot_cores <- 1
    n_cores <- min(tot_cores, specs$n_cores)

    # Fit the model, collecting errors/warnings as strings rather than
    # letting a bad dataset halt execution without reporting anything back
    stop_strs <- character(0)
    warn_strs <- character(0)
    daily <- NULL
    inst <- NULL
    withCallingHandlers(
      tryCatch({
        if (!requireNamespace("rstan", quietly = TRUE)) stop("rstan is required but not installed. Install it with: install.packages('rstan')")

        consolelog <- utils::capture.output(
          stanfit <- rstan::stan(
            file=specs$model_path, data=data_list, pars=specs$params_out,
            chains=specs$n_chains, cores=n_cores,
            iter=specs$burnin_steps + specs$saved_steps, warmup=specs$burnin_steps,
            thin=specs$thin_steps, verbose=specs$verbose, open_progress=FALSE),
          split=specs$verbose)

        if(stanfit@mode == 2L) {
          stop(paste(utils::capture.output(print(stanfit)), collapse='\n'))
        }

        # format the Stan summary matrix into per-variable data.frames
        stan_mat <- rstan::summary(stanfit)$summary
        mcmc_out <- format_mcmc_mat_nosplit(
          stan_mat, data_list$n_days, data_list$n_obs, specs$model_name,
          keep_mcmc=isTRUE(specs$keep_mcmcs), stanfit)

        # daily GPP/ER/K600 estimates: join Stan's date_index back to dates
        date_index <- time_index <- index <- '.dplyr.var'
        daily <- mcmc_out$daily %>%
          dplyr::left_join(date_df, by='date_index') %>%
          dplyr::select(-date_index, -time_index, -index) %>%
          dplyr::select(date, dplyr::everything())

        # instantaneous DO.mod.down estimates come from the 'metab' Stan
        # transformed parameter (posterior median), which format_mcmc_mat_nosplit()
        # buckets by row count rather than by name since 'metab' isn't in its
        # par_homes lookup table; find that bucket by its column names instead
        is_metab_bucket <- sapply(mcmc_out, function(df) is.data.frame(df) && any(grepl('^metab_', names(df))))
        metab_bucket_name <- names(mcmc_out)[is_metab_bucket][1]
        if(is.na(metab_bucket_name)) {
          stop("could not find 'metab' in the Stan output; check that specs$params_out includes 'metab'")
        }
        inst <- mcmc_out[[metab_bucket_name]] %>%
          dplyr::select(date_index, time_index, DO.mod.down=metab_50pct) %>%
          dplyr::inner_join(obs_index_df, by=c('date_index','time_index')) %>%
          dplyr::select(solar.time, DO.obs.down, DO.mod.down) %>%
          dplyr::arrange(solar.time)

      }, error=function(err) {
        stop_strs <<- c(stop_strs, err$message)
      }), warning=function(war) {
        warn_strs <<- c(warn_strs, war$message)
        invokeRestart("muffleWarning")
      })

    # if fitting failed, fill in NA daily estimates (with real dates) so the
    # returned model at least reports which dates were attempted
    if(length(stop_strs) > 0 || is.null(daily)) {
      na_vec <- rep(as.numeric(NA), nrow(date_df))
      daily <- data.frame(
        date=date_df$date,
        GPP_daily_2.5pct=na_vec, GPP_daily_50pct=na_vec, GPP_daily_97.5pct=na_vec,
        ER_daily_2.5pct=na_vec, ER_daily_50pct=na_vec, ER_daily_97.5pct=na_vec,
        K600_daily_2.5pct=na_vec, K600_daily_50pct=na_vec, K600_daily_97.5pct=na_vec)
      inst <- NULL
    }
    daily <- dplyr::mutate(daily, valid_day=TRUE, warnings='', errors='')

    fit <- list(
      daily=daily, inst=inst,
      warnings=trimws(unique(warn_strs)), errors=trimws(unique(stop_strs)))
  })

  # Package and return results
  mm <- metab_model(
    model_class="metab_bayes_2s",
    info=info,
    fit=fit,
    log=NULL,
    mcmc=if(isTRUE(specs$keep_mcmcs)) stanfit else NULL,
    mcmc_data=if(isTRUE(specs$keep_mcmc_data)) data_list else NULL,
    fitting_time=fitting_time,
    compile_time=system.time({}), # rstan::stan() compiles & samples in one call; not timed separately
    specs=specs,
    data=dat_list$data, # keep the units if given
    data_daily=dat_list$data_daily)

  # Update data with DO predictions
  success <- !is.null(fit$inst) && length(fit$errors) == 0
  if(success) {
    mm@data <- predict_DO(mm)
  } else {
    warntxt <- paste0(
      'Modeling failed\n',
      if(length(fit$warnings) > 0) paste0('  Warnings:\n', paste0('    ', fit$warnings, collapse='\n')),
      if(length(fit$errors) > 0) paste0('  Errors:\n', paste0('    ', fit$errors, collapse='\n')))
    warning(warntxt)
  }

  # Return
  mm
}


#' Reshape long-format two-station data into the list expected by the
#' two-station Stan model
#'
#' Time-shifts the upstream DO series to match the travel time between
#' stations, then pivots the result into the \code{n_obs x n_days} matrices
#' expected by the \code{data} block of \code{inst/models/b2_np_oi_tr_plrckm.stan}
#' (see \code{\link{metab_bayes_2s}}).
#'
#' The upstream observation that "matches" a given downstream observation at
#' row \code{i} was recorded \code{lag[i] <- round(travel.time[i] /
#' timestep_days)} timesteps earlier, where \code{timestep_days} is the
#' median timestep of \code{data$solar.time}. This must be computed the same
#' way as in \code{\link{mm_validate_data_2station}}'s lead-in check, so that the
#' \code{max(lag)} computed here always agrees with the lead-in requirement
#' already validated there. The first \code{max(lag)} rows of \code{data} are
#' lead-in rows: they supply upstream DO for the shift but are never
#' themselves treated as modeled (downstream) observations.
#'
#' @param data data.frame as validated by \code{\link{mm_validate_data}} for
#'   \code{\link{metab_bayes_2s}}: must contain \code{solar.time},
#'   \code{DO.obs.up}, \code{DO.sat.up}, \code{DO.obs.down},
#'   \code{DO.sat.down}, \code{light}, \code{depth}, \code{temp.water},
#'   \code{travel.time}, sorted ascending by \code{solar.time}, and must
#'   include the lead-in rows required to cover the longest travel time (see
#'   \code{\link{metab_bayes_2s}}).
#' @param specs a list of model specs (see \code{\link{specs}}), expected to
#'   already contain \code{K600_lnorm_meanlog} and \code{K600_lnorm_sdlog}
#'   -- e.g., the object returned by \code{specs(mm_name('bayes_2s'))}, which
#'   populates them with sensible defaults. This function does not supply
#'   its own fallback values; if \code{specs} is omitted or missing these
#'   fields, the resulting Stan data will contain NULL/missing values for
#'   them.
#' @return a named list with all variables in the Stan model's data block:
#'   \code{n_obs}, \code{n_days}, \code{DO_obs_up}, \code{DO_sat_up},
#'   \code{DO_obs_down}, \code{DO_sat_down}, \code{light}, \code{depth},
#'   \code{temp_water}, \code{travel_time} (each an \code{n_obs x n_days}
#'   matrix, unitless), and \code{K600_lnorm_meanlog}/\code{K600_lnorm_sdlog}
#' @importFrom unitted v
#' @keywords internal
prepdata_bayes_2s <- function(data, specs=NULL) {

  # strip units; Stan cannot handle unitted vectors/matrices
  data <- v(data)

  # timestep_days must match mm_validate_data_2station()'s lead-in check
  # exactly (median timestep, in days), so that max_lag here agrees with
  # what was validated there
  timestep_days <- stats::median(as.numeric(diff(data$solar.time), units='days'))

  # lag, in timesteps, that the upstream series must be shifted by to line up
  # with each row's downstream observation
  lag <- round(data$travel.time / timestep_days)
  max_lag <- max(lag)
  n_total <- nrow(data)

  # the first max_lag rows are lead-in only (upstream DO used for the shift,
  # but never modeled themselves). every row i > max_lag is guaranteed to
  # have a valid shift target (i - lag[i] >= 1) because lag[i] <= max_lag
  keep <- seq.int(max_lag + 1, n_total)
  shift_idx <- keep - lag[keep]
  if(any(shift_idx < 1)) {
    # should be unreachable given max_lag's definition; guards against
    # programming errors rather than expected user input
    stop('internal error: upstream shift index falls before the first row of data')
  }

  modeled <- data.frame(
    solar.time = data$solar.time[keep],
    DO_obs_up = data$DO.obs.up[shift_idx],
    DO_sat_up = data$DO.sat.up[shift_idx],
    DO_obs_down = data$DO.obs.down[keep],
    DO_sat_down = data$DO.sat.down[keep],
    light = data$light[keep],
    depth = data$depth[keep],
    temp_water = data$temp.water[keep],
    travel_time = data$travel.time[keep]
  )

  # pivot into n_obs x n_days matrices, one column per unique date, following
  # the same time_by_date_matrix approach used for the 1-station Stan models
  # (see prepdata_bayes() in metab_bayes.R)
  date_vec <- as.character(as.Date(modeled$solar.time))
  date_table <- table(date_vec)
  n_days <- length(date_table)
  n_obs_per_day <- unique(unname(date_table))
  if(length(n_obs_per_day) > 1) {
    stop(
      'dates have differing numbers of modeled rows after lead-in removal; ',
      'observations cannot be combined into a matrix: ',
      paste(sprintf('%s (%d rows)', names(date_table), date_table), collapse=', '))
  }
  n_obs <- n_obs_per_day

  to_matrix <- function(vec) matrix(vec, nrow=n_obs, ncol=n_days, byrow=FALSE)

  # confirm each date occupies a contiguous block of rows, i.e., that data
  # was sorted by solar.time; otherwise the matrix pivot below would silently
  # scramble which rows belong to which date
  date_mat <- to_matrix(date_vec)
  unique_dates_per_col <- apply(date_mat, MARGIN=2, FUN=unique)
  if(is.list(unique_dates_per_col) || !isTRUE(all.equal(unname(unique_dates_per_col), names(date_table)))) {
    stop('data must be sorted by solar.time so that each date occupies a contiguous block of rows')
  }

  list(
    n_obs = n_obs,
    n_days = n_days,
    DO_obs_up = to_matrix(modeled$DO_obs_up),
    DO_sat_up = to_matrix(modeled$DO_sat_up),
    DO_obs_down = to_matrix(modeled$DO_obs_down),
    DO_sat_down = to_matrix(modeled$DO_sat_down),
    light = to_matrix(modeled$light),
    depth = to_matrix(modeled$depth),
    temp_water = to_matrix(modeled$temp_water),
    travel_time = to_matrix(modeled$travel_time),
    K600_lnorm_meanlog = specs$K600_lnorm_meanlog,
    K600_lnorm_sdlog = specs$K600_lnorm_sdlog
  )
}


#### metab_bayes_2s class ####

#' Metabolism model fitted by two-station (Variable Flow Two-Station) Bayesian
#' MCMC
#'
#' \code{metab_bayes_2s} models use Bayesian MCMC methods to fit values of
#' GPP, ER, and K600 from paired upstream/downstream DO curves. This class
#' inherits from \code{metab_bayes} (same \code{log}/\code{mcmc}/
#' \code{mcmc_data}/\code{compile_time} slots, and therefore the same
#' \code{\link{get_mcmc}}, \code{\link{get_mcmc_data}}, and
#' \code{\link{get_log}} methods), but \code{predict_metab} and
#' \code{predict_DO} are overridden below because two-station's fitted-value
#' structure and output columns differ from one-station's.
#'
#' @exportClass metab_bayes_2s
#' @family metab.model.classes
setClass("metab_bayes_2s", contains="metab_bayes")


#' @describeIn get_params Does the same Stan-output-to-streamMetabolizer
#'   renaming as \code{get_params.metab_bayes}, but (unlike that method) does
#'   not delegate the rest of the work to \code{get_params.metab_model} via
#'   \code{NextMethod()}: that generic implementation looks up parameter
#'   names via \code{get_param_names()}, which builds an ODE-based dDOdt
#'   function using streamMetabolizer's one-station instantaneous-rate
#'   framework (\code{ode_method}/\code{GPP_fun}/\code{ER_fun}/
#'   \code{deficit_src}) -- machinery that doesn't apply to the two-station
#'   steady-state model's daily GPP/ER/K600 parameters. \code{fixed}
#'   column/star annotations (relevant only to models that can take fixed
#'   daily parameters from \code{data_daily}) are not supported here.
#' @export
#' @import dplyr
get_params.metab_bayes_2s <- function(
  metab_model, date_start=NA, date_end=NA, uncertainty=c('sd','ci','none'), messages=TRUE, ...) {

  uncertainty <- match.arg(uncertainty)

  fit <- metab_model@fit$daily
  if(is.null(fit)) return(NULL)

  # Stan prohibits '.' in variable names, so convert back from '_' to '.',
  # as in get_params.metab_bayes
  parnames <- setNames(gsub('_', '\\.', metab_model@specs$params_out), metab_model@specs$params_out)
  parnames <- parnames[order(nchar(parnames), decreasing=TRUE)]
  for(i in seq_along(parnames)) {
    names(fit) <- gsub(names(parnames[i]), parnames[[i]], names(fit))
  }
  names(fit) <- gsub('_mean$', '', names(fit))
  names(fit) <- gsub('_sd$', '.sd', names(fit))
  names(fit) <- gsub('_50pct$', '.median', names(fit))
  names(fit) <- gsub('_2.5pct$', '.lower', names(fit))
  names(fit) <- gsub('_97.5pct$', '.upper', names(fit))

  fit <- mm_filter_dates(fit, date_start=date_start, date_end=date_end)

  metab.vars <- c('GPP.daily', 'ER.daily', 'K600.daily')
  for(mv in metab.vars) {
    if(paste0(mv, '.median') %in% names(fit)) fit[[mv]] <- fit[[paste0(mv, '.median')]]
  }
  keep.cols <- c('date', unlist(lapply(metab.vars, function(mv) grep(paste0('^', mv, '($|\\.)'), names(fit), value=TRUE))))
  params <- fit[intersect(keep.cols, names(fit))]

  params <- switch(
    uncertainty,
    'none' = params[!grepl('\\.median$|\\.sd$|\\.lower$|\\.upper$', names(params))],
    'sd'   = params[!grepl('\\.median$|\\.lower$|\\.upper$', names(params))],
    'ci'   = params[!grepl('\\.median$|\\.sd$', names(params))])

  # attach raw warnings/errors columns (not yet compressed into a single
  # column); show()'s pretty_print_ddat()/compress_msgs() does that
  # compression itself at print time, as in get_params.metab_model
  if(messages && exists('date', fit) && any(c('warnings','errors') %in% names(fit))) {
    msgs <- fit[c('date','warnings','errors') %>% { .[. %in% names(fit)] }]
    params <- left_join(params, msgs, by='date', copy=TRUE)
  }

  params
}


#' @describeIn predict_metab Pulls daily GPP, ER, and K600 estimates out of
#'   the two-station Stan model results.
#' @export
#' @import dplyr
predict_metab.metab_bayes_2s <- function(metab_model, date_start=NA, date_end=NA, ...) {

  Var1 <- Var2 <- '.dplyr.var'
  fit.names <- expand.grid(c('50pct','2.5pct','97.5pct'), c('GPP_daily','ER_daily','K600_daily'), stringsAsFactors=FALSE) %>%
    select(Var2, Var1) %>%
    apply(MARGIN=1, FUN=function(row) do.call(paste, c(as.list(row), list(sep='_'))))
  metab.names <- expand.grid(c('','.lower','.upper'), c('GPP','ER','K600'), stringsAsFactors=FALSE) %>%
    select(Var2, Var1) %>%
    apply(MARGIN=1, FUN=function(row) do.call(paste0, as.list(row)))

  fit <- metab_model@fit$daily %>%
    mm_filter_dates(date_start=date_start, date_end=date_end)
  if(is.null(fit) || !all(fit.names %in% names(fit))) {
    stop('could not find GPP_daily, ER_daily, and K600_daily estimates in the model fit')
  }
  preds <- fit[c('date', fit.names)] %>%
    setNames(c('date', metab.names))

  # add date-specific fitting warnings/errors, as in predict_metab.metab_bayes
  warnings <- errors <- '.dplyr.var'
  if(!is.null(fit) && all(c('date','warnings','errors') %in% names(fit))) {
    messages <- fit %>%
      select(date, warnings, errors) %>%
      compress_msgs('msgs.fit', warnings.overall=metab_model@fit$warnings, errors.overall=metab_model@fit$errors)
    preds <- full_join(preds, messages, by='date', copy=TRUE)
  } else {
    preds <- mutate(preds, msgs.fit=NA)
  }

  preds <- mutate(
    preds,
    warnings=if(length(metab_model@fit$errors) > 0) NA else '',
    errors=if(length(metab_model@fit$errors) > 0) NA else '')

  preds
}


#' @describeIn predict_DO Two-station (Variable Flow Two-Station) models.
#'   Returns a data.frame with columns \code{solar.time}, \code{DO.obs.down}
#'   (the observed downstream DO from the input data), and \code{DO.mod.down}
#'   (the posterior median of the two-station Stan model's fitted downstream DO)
#'   -- unlike the one-station \code{predict_DO} methods, which return
#'   \code{DO.obs}/\code{DO.mod}. The values are those computed once at fitting
#'   time (see \code{\link{metab_bayes_2s}}); \code{use_saved=FALSE} (on-demand
#'   recomputation from the fitted daily GPP/ER/K600 medians) is not
#'   implemented.
#' @export
predict_DO.metab_bayes_2s <- function(metab_model, date_start=NA, date_end=NA, ..., use_saved=TRUE) {

  if(!isTRUE(use_saved)) {
    stop("predict_DO(use_saved=FALSE) is not implemented for metab_bayes_2s; only the fitted-time DO.mod.down values are available")
  }

  inst <- metab_model@fit$inst
  if(is.null(inst)) {
    stop("no DO.mod.down predictions are available; the model fit may have failed (see get_fit(metab_model))")
  }

  # NOTE: R/plot_DO_preds.R and tests/testthat/helper-rmse_DO.R both
  # hard-code the one-station DO.obs/DO.mod column names; they still need to
  # branch on (or be parameterized for) DO.obs.down/DO.mod.down before those
  # tools will work with two-station predictions -- deferred, out of scope
  # for this PR.
  mm_filter_dates(inst, date_start=date_start, date_end=date_end)
}
