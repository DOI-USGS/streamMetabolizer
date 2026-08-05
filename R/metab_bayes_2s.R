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
#' @section Two-station day window: Two-station days begin at 06:00 and run a
#'   full 24 hours, to 06:00 the next day. This convention is unrelated to the
#'   4 AM / 28-hour \code{day_start}/\code{day_end} window used by the
#'   one-station models, which describes an overlapping diel window rather
#'   than a partition of the time series; the two are not interchangeable.
#'   Days that do not fill the window -- at the edges of a dataset whose
#'   bounds don't fall on 06:00, or where observations are missing -- are
#'   dropped with a message.
#'
#' @section Two-station data requirements: In addition to the checks
#'   performed by \code{\link{mm_validate_data}}, \code{data$travel.time} (the
#'   reach travel time between the upstream and downstream stations, in days)
#'   must be strictly positive, and at least one row must have enough
#'   preceding observations of upstream DO to cover its own travel time. Rows
#'   at the start of \code{data} that lack that lead-in are not an error: they
#'   serve as lead-in only, supplying upstream DO for later rows without being
#'   modeled themselves.
#'
#'   Travel time is separately subject to a ceiling, \code{specs$
#'   max_travel_time_hours} (10 hours by default, configurable up to 12).
#'   Beyond it, a day's upstream parcel originates before the day's own 06:00
#'   start, where the day-normalized light it experienced no longer has a
#'   well-defined day total to be a proportion of. Days exceeding the ceiling
#'   are dropped with a message rather than treated as an error, since the
#'   remaining days are unaffected. A travel time far above the ceiling
#'   usually means the column was supplied in the wrong units -- days are
#'   expected, not minutes or hours.
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

    # Determine the same modeled-row index set that prepdata_bayes_2s()
    # applies internally, so that Stan's date_index/time_index can be mapped
    # back to actual dates and solar.times for the daily and instantaneous
    # results below. Both call mm_align_2s() rather than reimplementing the
    # lag/day-window math, so the two can't disagree about which rows are
    # modeled (prepdata_bayes_2s() returns only the Stan-ready matrices, so
    # the alignment itself isn't available to read back off its result).
    aln <- mm_align_2s(data_v, max_travel_time_hours=specs$max_travel_time_hours)
    keep <- aln$keep
    modeled_solar_time <- data_v$solar.time[keep]
    date_df <- tibble::tibble(date=unique(aln$date), date_index=seq_len(aln$n_days))
    obs_index_df <- tibble::tibble(
      solar.time=modeled_solar_time,
      DO.obs.down=data_v$DO.obs.down[keep],
      date_index=rep(date_df$date_index, each=aln$n_obs),
      time_index=rep(seq_len(aln$n_obs), times=aln$n_days))

    # Prepare the Stan data list (matrices from data, plus scalar priors from
    # specs). modifyList (not c()) is used because prepdata_bayes_2s() already
    # supplies K600_lnorm_meanlog/K600_lnorm_sdlog (read from specs), and
    # those two names are also in specs$params_in; a plain c() would create
    # duplicate-named list elements instead of overriding
    data_list <- prepdata_bayes_2s(dat_list$data, specs=specs, aln=aln)
    data_list <- modifyList(data_list, specs[specs$params_in])

    # Check and parse model file path
    specs$model_path <- mm_locate_filename(specs$model_name)

    # determine how many cores to use, as in runstan_bayes()
    n_cores <- mm_determine_cores(specs$n_cores, n_chains=specs$n_chains, verbose=specs$verbose)

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
          # mirror runstan_bayes()'s warn-and-continue pattern for a failed
          # run: report the diagnostic as a warning (caught below) and skip
          # the post-processing steps, which all assume a successful fit
          warning(paste(utils::capture.output(print(stanfit)), collapse='\n'))
        } else {
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
        }

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
#' The alignment itself -- the per-row lag, the per-row lead-in test, the
#' 06:00 day window, the travel-time ceiling, and the whole-day completeness
#' requirement -- is computed by \code{mm_align_2s} (see \code{mm_lag_2s.R}),
#' which \code{\link{metab_bayes_2s}} and \code{\link{mm_validate_data_2station}}
#' also route through. This function only applies the resulting indices and
#' pivots the result into matrices.
#'
#' Rows at the start of \code{data} whose upstream shift would reach before
#' the first row are lead-in rows: they supply upstream DO for the shift but
#' are never themselves treated as modeled (downstream) observations.
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
#' @param aln optional, the alignment already computed for \code{data} by
#'   \code{mm_align_2s} (see \code{mm_lag_2s.R}). \code{\link{metab_bayes_2s}}
#'   needs the same alignment to map Stan's indices back to dates, so it
#'   computes it once and passes it here; leave \code{NULL} to compute it.
#'   Supplying an alignment that doesn't correspond to \code{data} will
#'   silently produce wrong matrices.
#' @return a named list with all variables in the Stan model's data block:
#'   \code{n_obs}, \code{n_days}, \code{DO_obs_up}, \code{DO_sat_up},
#'   \code{DO_obs_down}, \code{DO_sat_down}, \code{light}, \code{depth},
#'   \code{temp_water}, \code{travel_time} (each an \code{n_obs x n_days}
#'   matrix, unitless), and \code{K600_lnorm_meanlog}/\code{K600_lnorm_sdlog}
#' @importFrom unitted v
#' @keywords internal
prepdata_bayes_2s <- function(data, specs=NULL, aln=NULL) {

  # strip units; Stan cannot handle unitted vectors/matrices
  data <- v(data)

  # per-row lag, per-row lead-in, 06:00 day window, travel-time ceiling, and
  # whole-day completeness -- all owned by mm_align_2s() so that this
  # function, metab_bayes_2s(), and mm_validate_data_2station() share one
  # definition. metab_bayes_2s() needs the same alignment to map Stan's
  # indices back to dates, so it computes it once and passes it in; recompute
  # it here only when called directly. specs may be NULL or lack the ceiling,
  # in which case mm_align_2s()'s own default applies
  if(is.null(aln)) {
    aln <- if(is.null(specs$max_travel_time_hours)) {
      mm_align_2s(data)
    } else {
      mm_align_2s(data, max_travel_time_hours=specs$max_travel_time_hours)
    }
  }
  keep <- aln$keep
  shift_idx <- aln$shift_idx

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

  # defensive check: every day mm_align_2s() considers complete must have
  # non-NA light. Stan rejects NA data outright, and since two-station fits
  # every date jointly, one bad day would fail the whole multi-day fit with
  # an opaque diagnostic that never mentions light or day boundaries -- name
  # the day here instead
  na_light <- is.na(modeled$light)
  if(any(na_light)) {
    bad_dates <- unique(aln$date[na_light])
    stop(paste0(
      'light is NA for ', length(bad_dates), ' day(s) that mm_align_2s() otherwise ',
      'considers complete: ', paste(bad_dates, collapse=', '), '. ',
      'Whatever computed the light column is using a different day window than ',
      'mm_align_2s(), or introduced other NAs -- check that before fitting; Stan ',
      'does not accept NA data.'), call.=FALSE)
  }

  # pivot into n_obs x n_days matrices, one column per two-station day, using
  # the same mm_time_by_date_matrix()/mm_check_dates_contiguous() helpers
  # shared with prepdata_bayes() (see mm_time_by_date_matrix.R). mm_align_2s()
  # has already guaranteed every day holds exactly n_obs rows
  date_vec <- as.character(aln$date)
  date_table <- table(date_vec)
  n_obs <- aln$n_obs
  n_days <- aln$n_days

  to_matrix <- mm_time_by_date_matrix(n_obs, n_days)

  # confirm each date occupies a contiguous block of rows, i.e., that data
  # was sorted by solar.time; otherwise the matrix pivot below would silently
  # scramble which rows belong to which date
  mm_check_dates_contiguous(
    to_matrix(date_vec), date_table,
    'data must be sorted by solar.time so that each date occupies a contiguous block of rows')

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


#' Compute the within-day light proportion for two-station light lagging
#'
#' Converts a single, already-combined light value per timestep into the
#' travel-time-weighted, within-day light proportion required by
#' \code{\link{metab_bayes_2s}}'s \code{light} column: Bishop et al. (2026)
#' Eq. 2, light summed over the upstream travel-time window ending at each
#' timestep, divided by that day's total light sum. See
#' \code{\link{two_station_example}} for the quantity this must match.
#'
#' Reuses \code{\link{mm_lag_2s}} for the per-row lag, window-start index,
#' and lead-in test, so this function's window can't drift from
#' \code{\link{mm_align_2s}}'s definition of the same lag. Rows lacking
#' lead-in (\code{shift_idx < 1}) get \code{NA} rather than a value computed
#' from a truncated window.
#'
#' This function expects a single, already-combined light value per
#' timestep.
#'
#' @section Day-sum denominator: the daily total in the denominator uses the
#'   same \code{\link{mm_date_2s}} 06:00-06:00 day window that
#'   \code{\link{mm_align_2s}}/\code{\link{metab_bayes_2s}} use elsewhere in
#'   the two-station pipeline, so any 06:00-day \code{\link{mm_align_2s}}
#'   counts as full is guaranteed to have every one of its rows counted as
#'   full here too -- its light proportions are never \code{NA}. (An earlier
#'   version used calendar midnight-to-midnight days instead, which could
#'   disagree with \code{\link{mm_align_2s}} at the edges of
#'   \code{solar.time}; see \code{test-metab_bayes_2s.R} for the regression
#'   test.)
#'
#'   A day that doesn't hold a full \code{round(1/timestep_days)} rows --
#'   at either end of \code{solar.time}, or from missing sensor data -- has
#'   no well-defined total to divide by, so every row in that day gets
#'   \code{NA} rather than a proportion computed from a partial sum.
#'
#' @section Regular-timestep requirement (temporary limitation): this
#'   function assumes every row of \code{solar.time} is exactly one nominal
#'   timestep apart, so that \code{\link{mm_lag_2s}}'s window-start indices
#'   can be used directly as row offsets. Real deployment data with gaps or
#'   multiple, mutually phase-shifted deployments breaks that assumption.
#'   TODO(#475-item7): snap-to-bin rounding would lift this restriction. For
#'   now, this function checks timestep regularity up front and errors
#'   rather than silently computing wrong windows on irregular input.
#'
#' @param solar.time POSIXct vector of timestamps, in UTC, sorted ascending.
#' @param light numeric vector, the same length as \code{solar.time}, of a
#'   single combined light value per timestep.
#' @param travel.time numeric vector of reach travel times, in days, the
#'   same length as \code{solar.time}.
#' @return a numeric vector, the same length as \code{solar.time}, of the
#'   within-day light proportion, or \code{NA} where \code{solar.time} lacks
#'   lead-in or falls in an incomplete 06:00-06:00 day.
#' @keywords internal
mm_lag_light_2s <- function(solar.time, light, travel.time) {

  # regularity guard: mm_lag_2s()'s shift_idx is a row offset, which only
  # matches a time offset when every row is exactly one nominal timestep
  # apart (see the roxygen section above)
  timesteps <- mm_get_timestep(solar.time, format='unique')
  if(length(timesteps) != 1) {
    stop(
      'mm_lag_light_2s() requires a regular timestep grid; solar.time has ',
      'gaps or multiple, phase-shifted deployments, which this function ',
      'does not yet support (see issue #475)', call.=FALSE)
  }

  lagged <- mm_lag_2s(solar.time, travel.time)
  n_total <- length(light)

  # per-row window sum; window width varies row to row with travel.time, so
  # this can't be a fixed-width rolling sum
  window_sum <- vapply(seq_len(n_total), function(i) {
    if(!lagged$has_leadin[i]) return(NA_real_)
    sum(light[lagged$shift_idx[i]:i], na.rm=TRUE)
  }, numeric(1))

  # 06:00-day total (mm_date_2s(), same day window as mm_align_2s()), NA'd
  # out for any day that doesn't fill the nominal timestep count (see the
  # roxygen section on the day-sum denominator)
  day_total <- '.dplyr.var'
  day <- mm_date_2s(solar.time)
  expected_n <- round(1 / lagged$timestep_days)
  day_totals <- tibble::tibble(day=day, light=light) %>%
    group_by(day) %>%
    mutate(day_total=ifelse(n() == expected_n, sum(light, na.rm=TRUE), NA_real_)) %>%
    ungroup() %>%
    pull(day_total)

  window_sum / day_totals
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
