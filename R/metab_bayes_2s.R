#' @include metab_model-class.R metab_bayes.R
#' @include mm_modeled_rows_2s.R mm_filter_valid_days_2s.R
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
#' @section Fitting all dates at once, or one at a time: With
#'   \code{specs$split_dates=FALSE} (the default) every date is fitted
#'   together in a single Stan call; with \code{TRUE}, each date is fitted in
#'   its own call. The difference that matters is the observation-error
#'   \code{sigma}: the Stan model carries one \code{sigma} shared across
#'   whatever dates it is given, so fitting dates one at a time estimates a
#'   separate \code{sigma} per date instead of one pooled across the record
#'   -- the same tradeoff one-station's \code{split_dates=TRUE} makes on a
#'   structurally identical model. Joint fitting is the default because the
#'   pooled estimate uses the whole record; splitting is useful when dates
#'   should not share an error estimate, or to keep one unfittable date from
#'   costing the whole run, since each date's errors are collected and
#'   reported against that date alone.
#'
#'   Fitting one date at a time does not mean slicing the data by date: each
#'   date's upstream observations are drawn from earlier rows, across the
#'   date boundary, exactly as in the joint fit.
#'
#' @inheritParams metab
#' @return A metab_bayes_2s object containing the fitted model. This object
#'   can be inspected with the functions in the
#'   \code{\link{metab_model_interface}} and also \code{\link{get_mcmc}}.
#'
#' @section Two-station day window: Two-station days run 06:00-06:00 (24
#'   hours) -- unrelated to the one-station models' overlapping 4 AM/28-hour
#'   \code{day_start}/\code{day_end} window, which is not a partition of the
#'   time series; the two are not interchangeable. Days that do not fill the
#'   window -- at the edges of a dataset whose bounds don't fall on 06:00, or
#'   where observations are missing and the gap was too long to bridge (see
#'   the next section) -- are dropped with a message and reported as invalid
#'   days in the results.
#'
#' @section Gap filling: Brief interruptions in the record are bridged by
#'   linear interpolation before day completeness is assessed, so a day
#'   marred by a short sensor dropout can be modeled rather than discarded.
#'   Missing timesteps and \code{NA} values are treated alike, and a gap is
#'   measured by the time it spans rather than by how many rows are missing
#'   across it. Runs longer than \code{specs$max_gap_hours} (1 hour by
#'   default, configurable up to 2) are left untouched, and the days holding
#'   them are dropped: the policy is deliberately two-tier, because the Stan
#'   model's fixed-shape matrices cannot represent a partially observed day.
#'   Gaps at the very start or end of the record are never filled, since
#'   interpolation there would mean extrapolating from one side.
#'
#'   Two consequences follow. First, the returned model's data -- what
#'   \code{\link{get_data}} and \code{\link{predict_DO}} report -- covers the
#'   filled timesteps too, so it can hold rows that were not in the
#'   \code{data} argument, and \code{predict_DO} returns predictions at those
#'   timesteps accordingly. Second, if \code{data} was prepared by
#'   \code{\link{mm_format_data_2s}} its gaps are already filled and nothing
#'   happens here; but data formatted by hand already holds \code{light} as
#'   the within-day proportion rather than raw light, and interpolating an
#'   already-normalized value is an approximation, since the day total it
#'   was divided by is itself short by the missing terms. Supplying raw
#'   light to \code{\link{mm_format_data_2s}} and letting it do the
#'   conversion avoids that.
#'
#' @section Two-station data requirements: In addition to the checks
#'   performed by \code{\link{mm_validate_data}}, \code{data$travel.time}
#'   (the reach travel time between stations, in days) must be strictly
#'   positive, and at least one row must have enough preceding observations
#'   of upstream DO to cover its own travel time. Rows at the start of
#'   \code{data} that lack that lead-in are not an error: they serve as
#'   lead-in only, supplying upstream DO for later rows without being
#'   modeled themselves.
#'
#'   Travel time is also subject to a ceiling, \code{specs$
#'   max_travel_time_hours} (10 hours by default, configurable up to 12):
#'   beyond it, a day's upstream parcel originates before the day's own
#'   06:00 start, where the light it experienced no longer has a
#'   well-defined day total to be a proportion of. Days exceeding the
#'   ceiling are dropped with a message rather than treated as an error,
#'   since the remaining days are unaffected. A travel time far above the
#'   ceiling usually means the column was supplied in the wrong units --
#'   days are expected, not minutes or hours.
#'
#' @section Day validity tests: the checks above concern a day's structure --
#'   does it fill its window, is its travel time usable. \code{specs$day_tests}
#'   covers the data itself, and defaults to
#'   \code{c('complete_data', 'pos_depth')}: a day is dropped if any value it
#'   would be modeled from is \code{NA}, or if its depth is not positive.
#'   Neither can be modeled around, and one bad day fails \emph{every} date of
#'   a joint fit with a diagnostic naming neither the column nor the day.
#'
#'   The values tested are the ones the day is modeled from, not the rows
#'   falling inside it: a day's upstream DO comes from one travel time
#'   earlier, routinely from the preceding day's rows.
#'
#'   Only a subset of the one-station tests is accepted: \code{full_day} and
#'   \code{even_timesteps} are rejected, each rejecting nearly every
#'   two-station day.
#'
#' @section Dropped days in the results: a day dropped for any reason above
#'   never reaches Stan, but is not omitted from the results either. It appears
#'   among the daily estimates as a \code{valid_day=FALSE} row with \code{NA}
#'   values and the reason in its \code{errors} column, so that every date
#'   supplied in \code{data} is accounted for. Instantaneous predictions cover
#'   only the modeled timesteps, a dropped day having none to report.
#'
#' @export
#' @family metab_model
#' @importFrom utils modifyList
metab_bayes_2s <- function(
  specs=specs(mm_name('bayes_2s')),
  data={
    d <- mm_data(solar.time, DO.obs.up, DO.sat.up, DO.obs.down, DO.sat.down,
                 light, depth, temp.water, travel.time)
    # light here is the within-day proportion computed by
    # mm_lag_light_2s(), not the raw PAR that mm_data()'s shared 'light'
    # template describes -- unitless, matching two_station_example.
    d$light <- u(v(d$light), NA)
    d
  },
  data_daily=mm_data(date, optional='all'),
  info=NULL
) {

  if(missing(specs)) {
    # if specs is left to the default, it gets confused about whether specs() is
    # the argument or the function. tell it which:
    specs <- streamMetabolizer::specs(mm_name('bayes_2s'))
  }

  stanfit <- NULL
  mcmc_data <- NULL
  bayes_log <- NULL
  compile_time <- system.time({})
  fitting_time <- system.time({
    # Check data for correct column names, units, and travel.time bounds
    # (mm_validate_data()), then check lead-in coverage
    # (mm_validate_data_2station()), before any data prep begins
    dat_list <- mm_validate_data(data, data_daily, 'metab_bayes_2s')

    # Bridge short gaps before anything measures day completeness. This
    # reassigns dat_list$data, and every consumer downstream -- the lead-in
    # check, mm_align_2s(), the fitting functions, and the @data slot that
    # predict_DO() reads back -- uses the filled frame, so the row indices
    # mm_align_2s() returns keep referring to the frame they were computed
    # from. Data prepared by mm_format_data_2s() has already been filled, so
    # this is a no-op for it; hand-formatted data is filled here or nowhere
    dat_list$data <- mm_fill_gaps_2s(dat_list$data, max_gap_hours=specs$max_gap_hours)

    mm_validate_data_2station(dat_list$data)

    # Determine the same modeled-row index set that prepdata_bayes_2s()
    # applies internally, so that Stan's date_index/time_index can be mapped
    # back to actual dates and solar.times for the daily and instantaneous
    # results. Both call mm_align_2s() rather than reimplementing the
    # lag/day-window math, so the two can't disagree about which rows are
    # modeled (prepdata_bayes_2s() returns only the Stan-ready matrices, so
    # the alignment itself isn't available to read back off its result).
    aln <- mm_align_2s(v(dat_list$data), max_travel_time_hours=specs$max_travel_time_hours)

    # Reject days whose modeled values fail specs$day_tests. Validation
    # upstream checks structure only, so this is the first look at the values
    filtered <- mm_filter_valid_days_2s(dat_list$data, aln, day_tests=specs$day_tests)

    # from both stages, so no day that failed to reach Stan goes unreported
    removed <- rbind(aln$removed, filtered$removed)
    removed <- removed[order(removed$date), , drop=FALSE]
    rownames(removed) <- NULL
    aln <- filtered$aln

    if(aln$n_days == 0) {
      stop(paste0(
        'no days remain after day_tests (',
        paste(specs$day_tests, collapse=', '), '): all ', nrow(filtered$removed),
        ' remaining day(s) were dropped. First few: ',
        paste(sprintf('%s (%s)',
                      as.character(utils::head(filtered$removed$date, 3)),
                      utils::head(filtered$removed$errors, 3)), collapse='; ')), call.=FALSE)
    }

    # Check and parse model file path
    specs$model_path <- mm_locate_filename(specs$model_name)

    # check the format of keep_mcmcs/keep_mcmc_data, coercing a date vector
    # to Date
    specs <- mm_check_keep_mcmc_specs(specs)

    if(isTRUE(specs$split_dates)) {
      # one date at a time. Each date's data preparation happens inside its
      # own error-collecting handler, so a single bad date is reported in that
      # date's errors column rather than halting the rest of the run
      perday <- bayes_perday_2s(dat_list$data, specs=specs, aln=aln)

      stanfit <- perday$mcmcs
      mcmc_data <- perday$mcmc_datas
      compile_time <- perday$compile_time
      bayes_log <- perday$log
      fit <- list(
        daily=perday$daily, inst=perday$inst,
        warnings=perday$warnings, errors=perday$errors)

    } else {
      # every date jointly, in a single Stan call.
      #
      # Prepare the Stan data list (matrices from data, plus scalar priors
      # from specs). modifyList (not c()) is used because prepdata_bayes_2s()
      # already supplies K600_lnorm_meanlog/K600_lnorm_sdlog (read from
      # specs), and those two names are also in specs$params_in; a plain c()
      # would create duplicate-named list elements instead of overriding.
      #
      # Prepared here, outside bayes_1fit_2s(), so that data problems
      # prepdata_bayes_2s() detects (NA light in an otherwise complete day,
      # non-contiguous dates) halt the joint fit with the original error
      # rather than being collected and reported as a failed fit
      data_list <- prepdata_bayes_2s(dat_list$data, specs=specs, aln=aln)
      data_list <- modifyList(data_list, specs[specs$params_in])

      # errors/warnings are collected as strings rather than letting a bad
      # dataset halt execution without reporting anything back
      fit1 <- bayes_1fit_2s(
        dat_list$data, aln=aln, specs=specs, data_list=data_list,
        keep_mcmc=isTRUE(specs$keep_mcmcs))
      stanfit <- fit1$stanfit
      mcmc_data <- if(isTRUE(specs$keep_mcmc_data)) data_list else NULL
      compile_time <- fit1$compile_time
      # mirrors metab_bayes()'s split_dates==FALSE naming (bayes_allply's
      # 'Compilation'/'MCMC_All_Days'), for the same single joint Stan call
      bayes_log <- setNames(list(fit1$compile_log, fit1$log), c('Compilation', 'MCMC_All_Days'))
      bayes_log <- bayes_log[!sapply(bayes_log, is.null)]

      # if fitting failed, fill in NA daily estimates (with real dates) so the
      # returned model at least reports which dates were attempted
      daily <- fit1$daily
      inst <- fit1$inst
      if(length(fit1$errors) > 0 || is.null(daily)) {
        daily <- mm_na_daily(fit1$date_df$date)
        inst <- NULL
      }
      daily <- dplyr::mutate(daily, valid_day=TRUE, warnings='', errors='')

      fit <- list(
        daily=daily, inst=inst,
        warnings=fit1$warnings, errors=fit1$errors)
    }

    # same in both modes: the days dropped before fitting come back as
    # valid_day=FALSE rows carrying the reason they were dropped
    fit$daily <- mm_rejoin_removed_days_2s(fit$daily, removed)
  })

  # Package and return results
  mm <- metab_model(
    model_class="metab_bayes_2s",
    info=info,
    fit=fit,
    log=bayes_log,
    # a single stanfit in joint mode, a list of them named by date in per-day
    # mode -- the same two shapes metab_bayes() already returns for
    # split_dates=FALSE/TRUE
    mcmc=stanfit,
    mcmc_data=mcmc_data,
    fitting_time=fitting_time - compile_time,
    compile_time=compile_time,
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


#### fitting helpers ####

#' Add pre-fit dropped days back into a two-station daily results frame
#'
#' Called once, where both fitting modes converge, so that neither can report
#' a dropped day differently from the other.
#'
#' @param daily the fitted daily results, already carrying
#'   \code{valid_day}/\code{warnings}/\code{errors} columns.
#' @param removed data.frame of \code{date} and \code{errors} for the dropped days.
#' @return \code{daily} with one added row per dropped day, sorted by date;
#'   columns present only in \code{daily} are \code{NA} on the added rows.
#' @keywords internal
mm_rejoin_removed_days_2s <- function(daily, removed) {

  if(nrow(removed) == 0) return(daily)

  # bind_rows, not rbind: the added rows carry only the GPP/ER/K600 quantile
  # columns, and the rest fill with NA rather than being enumerated here
  na_rows <- mm_na_daily(removed$date)
  na_rows$valid_day <- FALSE
  na_rows$warnings <- ''
  na_rows$errors <- removed$errors

  date <- '.dplyr.var'
  dplyr::arrange(dplyr::bind_rows(daily, na_rows), date)
}


#' Run one two-station Stan fit over one alignment
#'
#' The single Stan call plus output formatting shared by
#' \code{\link{metab_bayes_2s}}'s joint (all-dates-at-once) fit and
#' \code{bayes_perday_2s}'s per-day loop. Nothing here depends on
#' \code{aln} covering more than one day, so the same code path serves both:
#' the joint fit passes the full alignment from \code{mm_align_2s}, and the
#' per-day loop passes a one-day slice of it.
#'
#' The Stan call itself -- compile (or load from the \code{.stanrds} cache),
#' sample, and format the summary matrix into per-variable data.frames -- is
#' delegated to \code{\link{runstan_bayes}}, which handles two-station data
#' unchanged. What remains here is the two-station-specific part: collecting
#' errors, joining Stan's indices back to dates and solar.times, and pulling
#' the modeled downstream DO out of the \code{metab} parameter.
#'
#' Errors and warnings are collected as strings rather than propagated, in
#' the same shape \code{\link{bayes_1ply}} uses for one-station, so that a
#' caller looping over days can record a bad day and continue.
#'
#' @param data data.frame as validated by \code{\link{mm_validate_data}} for
#'   \code{\link{metab_bayes_2s}}, units optional. This is the \emph{full}
#'   dataset in both modes, never a per-day slice: \code{aln}'s \code{keep}
#'   and \code{shift_idx} index into it, and \code{shift_idx} routinely
#'   reaches back into the previous day's rows for upstream values.
#' @param aln an alignment as returned by \code{mm_align_2s} (see
#'   \code{mm_lag_2s.R}), either the whole thing or a single-day slice of it
#'   (see \code{bayes_perday_2s}).
#' @param specs a list of model specs (see \code{\link{specs}}), including
#'   \code{model_path} as resolved by \code{\link{mm_locate_filename}}. Passed
#'   on to \code{\link{runstan_bayes}}, which reads the chain/step/core
#'   settings from it.
#' @param data_list optional, the Stan data list already prepared for
#'   \code{data}/\code{aln} by \code{\link{prepdata_bayes_2s}} (plus the
#'   \code{specs$params_in} scalars). Supply it to keep data-preparation
#'   errors propagating to the caller; leave \code{NULL} to prepare it here,
#'   inside the error-collecting handler, so that preparation failures are
#'   reported as a failed fit instead. Supplying a list that doesn't
#'   correspond to \code{data}/\code{aln} will silently produce wrong results.
#' @param keep_mcmc logical. Retain the \code{stanfit} object in the result?
#'   Resolved by the caller, which knows whether \code{specs$keep_mcmcs} names
#'   this particular date.
#' @return a list with \code{daily} (data.frame of per-date estimates, or
#'   \code{NULL} if the fit failed), \code{inst} (data.frame of
#'   \code{solar.time}/\code{DO.obs.down}/\code{DO.mod.down}, or \code{NULL}),
#'   \code{stanfit} (\code{NULL} if the run failed or \code{keep_mcmc} is
#'   \code{FALSE}), \code{data_list},
#'   \code{compile_time}, \code{log}/\code{compile_log} (as returned by
#'   \code{\link{runstan_bayes}}; collected into the model object's log by the
#'   caller), \code{date_df} (the date/date_index lookup, available even on failure),
#'   and \code{warnings}/\code{errors} character vectors
#' @importFrom utils modifyList
#' @keywords internal
bayes_1fit_2s <- function(data, aln, specs, data_list=NULL, keep_mcmc=TRUE) {

  data_v <- v(data)
  keep <- aln$keep

  # map Stan's date_index/time_index back to real dates and solar.times.
  # mm_align_2s() guarantees each day holds exactly n_obs rows, in ascending
  # solar.time order, so the indices tile day-major over keep
  date_df <- tibble::tibble(date=unique(aln$date), date_index=seq_len(aln$n_days))
  obs_index_df <- tibble::tibble(
    solar.time=data_v$solar.time[keep],
    DO.obs.down=data_v$DO.obs.down[keep],
    date_index=rep(date_df$date_index, each=aln$n_obs),
    time_index=rep(seq_len(aln$n_obs), times=aln$n_days))

  stop_strs <- character(0)
  warn_strs <- character(0)
  daily <- NULL
  inst <- NULL
  stanfit <- NULL
  stan_log <- NULL
  compile_log <- NULL
  compile_time <- system.time({})

  withCallingHandlers(
    tryCatch({
      if (!requireNamespace("rstan", quietly = TRUE)) stop("rstan is required but not installed. Install it with: install.packages('rstan')")

      if(is.null(data_list)) {
        data_list <- prepdata_bayes_2s(data, specs=specs, aln=aln)
        data_list <- modifyList(data_list, specs[specs$params_in])
      }

      # Compile, sample, and format via one-station's runstan_bayes(), which
      # handles all three identically for two-station data: it takes the
      # nosplit branch (specs$split_dates is always FALSE here) and its
      # mm_compile_stan_model() call is the same .stanrds cache two-station
      # already used. runstan_bayes() also emits the mode==2L
      # failed-run warning, which the handler below collects.
      #
      # CAUTION: this delegation is only safe because
      # format_mcmc_mat_nosplit()'s data_list_d/data_list_n parameters are
      # declared but never referenced in its body -- runstan_bayes() passes
      # data_list$d/$n, which do not exist in a two-station data list (whose
      # Stan data block names them n_days/n_obs) and so arrive as NULL. If
      # anyone ever wires those two parameters up to actually be used,
      # two-station breaks here, silently and without an error: check this
      # call site before doing that.
      #
      # split_dates is forced FALSE for this call regardless of the user's
      # setting, because it selects runstan_bayes()'s output *formatter*, not
      # its looping behavior: every two-station fit -- the joint one, and each
      # date of bayes_perday_2s()'s loop -- is a single Stan call whose summary
      # must go through format_mcmc_mat_nosplit(). format_mcmc_mat_split()
      # collapses the summary to one row and strips a '[1]' suffix, which would
      # mangle metab[i,t] and leave no 'daily' element for the caller to find,
      # silently turning every fit into an apparent failure.
      specs$split_dates <- FALSE
      specs$keep_mcmc <- keep_mcmc
      stan_out <- do.call(runstan_bayes, c(list(data_list=data_list), specs))

      compile_time <- stan_out$compile_time
      stan_log <- stan_out$log
      compile_log <- stan_out$compile_log
      stanfit <- stan_out$mcmcfit

      # on a failed run runstan_bayes() warns (collected below) and returns
      # no parameter data.frames, so there is nothing to reshape
      if(!is.null(stan_out$daily)) {

        # daily GPP/ER/K600 estimates: join Stan's date_index back to dates
        date_index <- time_index <- index <- '.dplyr.var'
        daily <- stan_out$daily %>%
          dplyr::left_join(date_df, by='date_index') %>%
          dplyr::select(-date_index, -time_index, -index) %>%
          dplyr::select(date, dplyr::everything())

        # instantaneous DO.mod.down estimates come from the 'metab' Stan
        # transformed parameter (posterior median), which format_mcmc_mat_nosplit()
        # buckets by row count rather than by name since 'metab' isn't in its
        # par_homes lookup table; find that bucket by its column names instead
        is_metab_bucket <- sapply(stan_out, function(df) is.data.frame(df) && any(grepl('^metab_', names(df))))
        metab_bucket_name <- names(stan_out)[is_metab_bucket][1]
        if(is.na(metab_bucket_name)) {
          stop("could not find 'metab' in the Stan output; check that specs$params_out includes 'metab'")
        }
        inst <- stan_out[[metab_bucket_name]] %>%
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

  list(
    daily=daily, inst=inst, stanfit=stanfit, data_list=data_list,
    compile_time=compile_time, log=stan_log, compile_log=compile_log,
    date_df=date_df,
    warnings=trimws(unique(warn_strs)), errors=trimws(unique(stop_strs)))
}

#' Fit the two-station Stan model one date at a time
#'
#' Loops over the dates in a two-station alignment, fitting each in its own
#' Stan call, and reassembles the per-date results into the same
#' \code{daily}/\code{inst} shape \code{\link{metab_bayes_2s}}'s joint fit
#' produces -- so \code{\link{predict_metab}}, \code{\link{predict_DO}}, and
#' \code{\link{get_params}} see a consistent structure either way.
#'
#' A failed date does not abort the run: each date's fit collects its own
#' errors and warnings (see \code{bayes_1fit_2s}), which are recorded in that
#' date's \code{warnings}/\code{errors} columns while the remaining dates
#' proceed. This mirrors \code{\link{bayes_1ply}}'s per-day behavior for
#' one-station models. Run-level \code{warnings}/\code{errors} are
#' deliberately left empty for date-specific failures, because
#' \code{\link{predict_metab}} blanks out \emph{every} date's estimates when
#' the run-level slots are non-empty.
#'
#' The dates come from \code{mm_align_2s}'s non-overlapping 06:00-06:00
#' partition, already computed and carried on \code{aln}; day membership is
#' not re-derived here. In particular this does not route through
#' \code{\link{mm_model_by_ply}}, whose overlapping \code{day_start}/
#' \code{day_end} diel window is a different partition of the same rows (see
#' the two-station day window section of \code{\link{metab_bayes_2s}}).
#'
#' @section Per-date sigma: fitting one date at a time estimates a separate
#'   observation-error \code{sigma} per date rather than one pooled across
#'   the record -- the same tradeoff described in
#'   \code{\link{metab_bayes_2s}}'s "Fitting all dates at once, or one at a
#'   time" section.
#'
#' @inheritParams bayes_1fit_2s
#' @param aln optional, the alignment already computed for \code{data} by
#'   \code{mm_align_2s}. Leave \code{NULL} to compute it here.
#' @return a list with \code{daily} (one row per date, with
#'   \code{valid_day}/\code{warnings}/\code{errors} columns), \code{inst}
#'   (all dates' instantaneous predictions, ordered by \code{solar.time}, or
#'   \code{NULL} if no date succeeded), \code{mcmcs} and \code{mcmc_datas}
#'   (named by date, honoring \code{specs$keep_mcmcs}/\code{keep_mcmc_data};
#'   \code{NULL} if no date was selected), \code{compile_time}, \code{log}
#'   (compile and per-date MCMC logs, named as in
#'   \code{\link{metab_bayes}}'s \code{split_dates=TRUE} mode),
#'   \code{dates_fit}, \code{dates_failed}, and empty run-level
#'   \code{warnings}/\code{errors}
#' @keywords internal
bayes_perday_2s <- function(data, specs, aln=NULL) {

  if(is.null(aln)) {
    aln <- mm_align_2s(v(data), max_travel_time_hours=specs$max_travel_time_hours)
  }

  # resolved once rather than per date: it doesn't depend on which date is
  # being fit. Core count is not hoisted the same way -- runstan_bayes()
  # resolves it per call, as it already does for each one-station ply
  if(is.null(specs$model_path)) specs$model_path <- mm_locate_filename(specs$model_name)

  dates <- unique(aln$date)

  per_date <- lapply(dates, function(dt) {

    # slice the alignment, not the data: keep/shift_idx index into the full
    # dataset, and shift_idx routinely points at the previous date's rows for
    # this date's upstream values
    rows <- which(aln$date == dt)
    aln_1 <- list(
      keep=aln$keep[rows],
      shift_idx=aln$shift_idx[rows],
      date=aln$date[rows],
      n_obs=aln$n_obs,
      n_days=1L,
      timestep_days=aln$timestep_days)

    # keep_mcmcs/keep_mcmc_data are each either a single logical or a vector
    # of dates naming which fits to retain, as in bayes_1ply()
    dt <- as.Date(dt)
    keep_mcmc <- if(is.logical(specs$keep_mcmcs)) {
      isTRUE(specs$keep_mcmcs)
    } else {
      isTRUE(dt %in% specs$keep_mcmcs)
    }
    keep_mcmc_dat <- if(is.logical(specs$keep_mcmc_data)) {
      isTRUE(specs$keep_mcmc_data)
    } else {
      isTRUE(dt %in% specs$keep_mcmc_data)
    }

    fit1 <- bayes_1fit_2s(data, aln=aln_1, specs=specs, keep_mcmc=keep_mcmc)

    failed <- length(fit1$errors) > 0 || is.null(fit1$daily)
    daily <- if(failed) mm_na_daily(dt) else fit1$daily

    # valid_day refers to whether mm_align_2s() accepted the date, not to
    # whether its fit converged -- every date reaching this loop passed the
    # day-window and travel-time checks, so it stays TRUE and a failure is
    # reported through the errors column, as in the joint fit
    daily <- dplyr::mutate(
      daily,
      valid_day=TRUE,
      warnings=paste0(fit1$warnings, collapse='; '),
      errors=paste0(fit1$errors, collapse='; '))

    list(
      daily=daily, inst=if(failed) NULL else fit1$inst,
      stanfit=fit1$stanfit,
      mcmc_data=if(keep_mcmc_dat) fit1$data_list else NULL,
      compile_time=fit1$compile_time, failed=failed,
      log=fit1$log, compile_log=fit1$compile_log)
  })
  names(per_date) <- as.character(dates)

  failed <- vapply(per_date, `[[`, logical(1), 'failed')

  # bind_rows (not rbind) because a failed date contributes only the nine
  # NA quantile columns while a successful one carries the full set of Stan
  # summary statistics; the missing columns fill with NA
  daily <- dplyr::bind_rows(lapply(per_date, `[[`, 'daily'))

  insts <- lapply(per_date, `[[`, 'inst')
  insts <- insts[!vapply(insts, is.null, logical(1))]
  inst <- if(length(insts) > 0) {
    dplyr::arrange(dplyr::bind_rows(insts), solar.time)
  } else NULL

  # only the first date should pay a compile cost; the rest load the
  # .stanrds cache and contribute a zero compile_time
  compile_time <- Reduce(`+`, lapply(per_date, `[[`, 'compile_time'), system.time({}))

  if(any(failed)) {
    message(paste0(
      'fit ', sum(!failed), ' of ', length(dates), ' date(s); ',
      sum(failed), ' failed: ', paste(as.character(dates[failed]), collapse=', ')))
  }

  # named lists keyed by date, with a NULL entry for any date not selected by
  # keep_mcmcs/keep_mcmc_data; NULL overall when no date was selected, so that
  # get_mcmc()/get_mcmc_data() report "nothing kept" the same way the joint
  # fit does
  drop_if_all_null <- function(x) if(all(vapply(x, is.null, logical(1)))) NULL else x
  mcmcs <- drop_if_all_null(lapply(per_date, `[[`, 'stanfit'))
  mcmc_datas <- drop_if_all_null(lapply(per_date, `[[`, 'mcmc_data'))

  # mirrors metab_bayes()'s split_dates==TRUE naming: only the date(s) that
  # triggered a compile contribute a 'Compilation' entry, and each date's Stan
  # log is kept as 'MCMC_<date>'
  compile_log <- lapply(per_date, `[[`, 'compile_log')
  compile_log <- compile_log[!vapply(compile_log, is.null, logical(1))]
  if(length(compile_log) > 0) names(compile_log) <- rep('Compilation', length(compile_log))
  log <- lapply(per_date, `[[`, 'log')
  log <- log[!vapply(log, is.null, logical(1))]
  log <- setNames(log, paste0('MCMC_', names(log)))
  bayes_log <- c(compile_log, log)

  list(
    daily=daily,
    inst=inst,
    mcmcs=mcmcs,
    mcmc_datas=mcmc_datas,
    compile_time=compile_time,
    log=bayes_log,
    dates_fit=dates[!failed],
    dates_failed=dates[failed],
    warnings=character(0),
    errors=character(0))
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
  modeled <- mm_modeled_rows_2s(data, aln)

  # defensive, for direct calls only: a fit routed through metab_bayes_2s()
  # has already dropped these days. Stan's own diagnostic for NA data names
  # neither column nor day, so name them here instead
  if(anyNA(modeled)) {
    bad_cols <- names(modeled)[vapply(modeled, anyNA, logical(1))]
    bad_dates <- unique(aln$date[!stats::complete.cases(modeled)])
    stop(paste0(
      'NAs in ', paste(bad_cols, collapse=', '), ' for ', length(bad_dates),
      ' day(s) that mm_align_2s() otherwise considers complete: ',
      paste(bad_dates, collapse=', '),
      '. Stan does not accept NA data; run mm_filter_valid_days_2s() first, or ',
      'fit via metab_bayes_2s(), which does.'), call.=FALSE)
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
    DO_obs_up = to_matrix(modeled$DO.obs.up),
    DO_sat_up = to_matrix(modeled$DO.sat.up),
    DO_obs_down = to_matrix(modeled$DO.obs.down),
    DO_sat_down = to_matrix(modeled$DO.sat.down),
    light = to_matrix(modeled$light),
    depth = to_matrix(modeled$depth),
    temp_water = to_matrix(modeled$temp.water),
    travel_time = to_matrix(modeled$travel.time),
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
#' @section Day-sum denominator: the daily total in the denominator uses the
#'   same \code{\link{mm_date_2s}} 06:00-06:00 day window that
#'   \code{\link{mm_align_2s}}/\code{\link{metab_bayes_2s}} use elsewhere in
#'   the two-station pipeline, so any 06:00-day \code{\link{mm_align_2s}}
#'   counts as full is guaranteed to have every one of its rows counted as
#'   full here too -- its light proportions are never \code{NA}.
#'
#'   A day that doesn't hold a full \code{round(1/timestep_days)} rows --
#'   at either end of \code{solar.time}, or from missing sensor data -- has
#'   no well-defined total to divide by, so every row in that day gets
#'   \code{NA} rather than a proportion computed from a partial sum.
#'
#' @section Snap-to-bin requirement: \code{\link{mm_lag_2s}} matches rows by
#'   timestep bin, so a missing bin inside a light-sum window simply
#'   contributes no term, and a missing bin at a window's start yields
#'   \code{NA} via \code{has_leadin} -- both without special-casing here.
#'   This function does not itself check that \code{solar.time} sits on a
#'   bin grid; \code{\link{mm_lag_2s}} does.
#'
#' @param solar.time POSIXct vector of timestamps, in UTC, sorted ascending,
#'   already snapped to a single nominal timestep via
#'   \code{\link{mm_snap_to_bin_2s}}.
#' @param light numeric vector, the same length as \code{solar.time}, of a
#'   single combined light value per timestep.
#' @param travel.time numeric vector of reach travel times, in days, the
#'   same length as \code{solar.time}.
#' @return a numeric vector, the same length as \code{solar.time}, of the
#'   within-day light proportion, or \code{NA} where \code{solar.time} lacks
#'   lead-in or falls in an incomplete 06:00-06:00 day.
#' @keywords internal
mm_lag_light_2s <- function(solar.time, light, travel.time) {

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
#'   renaming as \code{get_params.metab_bayes}, but does not delegate to
#'   \code{get_params.metab_model} via \code{NextMethod()}: the two-station
#'   steady-state model's daily GPP/ER/K600 parameters don't fit that
#'   generic's one-station, ODE-based parameter-name lookup. \code{fixed}
#'   column/star annotations (relevant only to models that can take fixed
#'   daily parameters from \code{data_daily}) are not supported here.
#' @export
#' @import dplyr
get_params.metab_bayes_2s <- function(
  metab_model, date_start=NA, date_end=NA, uncertainty=c('sd','ci','none'), messages=TRUE, ...) {

  # not delegated to get_params.metab_model via NextMethod(): that generic's
  # parameter-name lookup builds an ODE-based dDOdt function from one-station's
  # ode_method/GPP_fun/ER_fun/deficit_src specs, which don't exist for this
  # steady-state model
  uncertainty <- match.arg(uncertainty)

  fit <- metab_model@fit$daily
  if(is.null(fit)) return(NULL)

  fit <- mm_rename_stan_params(fit, metab_model@specs$params_out)

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
#' @importFrom lifecycle deprecated is_present
#' @importFrom unitted get_units u
predict_metab.metab_bayes_2s <- function(metab_model, date_start=NA, date_end=NA, ..., attach.units=deprecated()) {

  Var1 <- Var2 <- '.dplyr.var'

  # check units-related arguments
  if (lifecycle::is_present(attach.units)) {
    unitted_deprecate_warn("predict_metab(attach.units)")
  } else {
    attach.units <- FALSE
  }

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

  preds <- mm_attach_fit_msgs(
    preds, fit,
    warnings.overall=metab_model@fit$warnings, errors.overall=metab_model@fit$errors)

  # attach.units if requested
  if(attach.units) {
    pred.units <- get_units(mm_data())[sapply(names(preds), function(x) strsplit(x, '\\.')[[1]][1], USE.NAMES=FALSE)]
    preds <- u(preds, pred.units)
  }
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
#'   implemented. Rows cover every timestep that was modeled, which includes
#'   any that were interpolated to bridge a short gap in the input data (see
#'   \code{\link{metab_bayes_2s}}'s gap-filling section), so \code{DO.obs.down}
#'   may hold an interpolated rather than an observed value at those
#'   timesteps.
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
