
# # NB: these lines work within testthat() calls:
# skip()
# skip_on_cran()
# skip_on_travis()
# skip_on_appveyor()
# skip_if_not_installed('deSolve')


# prepdata_bayes() ---------------------------------------------------------
#
# Direct, fast (no Stan) unit tests for prepdata_bayes()'s matrix-pivot and
# contiguous-sort-check behavior. Written to pin current behavior ahead of
# extracting shared logic with prepdata_bayes_2s() (see mm_time_by_date_matrix.R).

# Build a minimal, traceable two-day data.frame: 4 rows/day, DO.obs ==
# original row index so the pivoted matrix's values can be checked exactly.
make_bayes_prepdata <- function(n_per_day=4, n_days=2, shuffle=FALSE) {
  n_total <- n_per_day * n_days
  solar.time <- as.POSIXct("2050-06-01 00:00:00", tz="UTC") +
    as.difftime(
      rep((seq_len(n_per_day) - 1) * (24 / n_per_day), n_days) +
        rep((seq_len(n_days) - 1) * 24, each=n_per_day),
      units="hours")
  dat <- data.frame(
    solar.time = solar.time,
    date = as.Date(solar.time),
    DO.obs = seq_len(n_total),
    DO.sat = rep(10, n_total),
    depth = rep(1, n_total),
    temp.water = rep(20, n_total),
    light = rep(100, n_total)
  )
  if(shuffle) {
    # interleave day1/day2 rows to break contiguity while keeping an equal
    # row count per date (so the earlier date_table-based check still passes
    # and only the contiguous-sort check is exercised)
    dat <- dat[order(rep(seq_len(n_per_day), n_days)), ]
  }
  dat
}

test_that("prepdata_bayes() pivots data into the expected n x d matrix (dims and values)", {
  sp <- specs(mm_name('bayes'))
  dat <- make_bayes_prepdata(n_per_day=4, n_days=2)

  out <- prepdata_bayes(data=dat, data_daily=NULL, ply_date=NA, specs=sp)

  expect_equal(out$d, 2)
  expect_equal(out$n, 4)
  expect_equal(dim(out$DO_obs), c(4, 2))
  expect_equal(dim(out$DO_sat), c(4, 2))
  expect_equal(dim(out$depth), c(4, 2))
  # values: matrix(vec, nrow=4, ncol=2, byrow=FALSE) fills column-wise, so
  # day 1 = original rows 1:4, day 2 = original rows 5:8
  expect_equal(out$DO_obs[,1], as.numeric(1:4))
  expect_equal(out$DO_obs[,2], as.numeric(5:8))
})

test_that("prepdata_bayes()'s contiguous-sort check passes for properly sorted data", {
  sp <- specs(mm_name('bayes'))
  dat <- make_bayes_prepdata(n_per_day=4, n_days=2)

  expect_silent(prepdata_bayes(data=dat, data_daily=NULL, ply_date=NA, specs=sp))
})

test_that("prepdata_bayes() errors for data that isn't sorted by date", {
  # HISTORY: prior to the mm_time_by_date_matrix()/mm_check_dates_contiguous()
  # extraction, this check was `if(!all.equal(unique_dates, names(date_table)))
  # stop("couldn't fit given dates into matrix")`, unguarded by isTRUE().
  # Confirmed empirically pre-refactor: whenever the dates were truly
  # non-contiguous, all.equal() returned a character vector (not TRUE), and
  # `!` on a character vector always errors with "invalid argument type" in
  # base R -- so the intended "couldn't fit given dates into matrix" message
  # was actually unreachable. The shared mm_check_dates_contiguous() helper
  # fixes this as a side effect of the extraction, so this test now
  # asserts the originally-intended message.
  sp <- specs(mm_name('bayes'))
  dat <- make_bayes_prepdata(n_per_day=4, n_days=2, shuffle=TRUE)

  expect_error(prepdata_bayes(data=dat, data_daily=NULL, ply_date=NA, specs=sp), "couldn't fit given dates into matrix")
})

test_that("format_mcmc_mat_nosplit()'s data_list_d/data_list_n parameters stay unused", {
  # bayes_1fit_2s() (two-station) delegates to runstan_bayes(), which calls
  # format_mcmc_mat_nosplit(stan_mat, data_list$d, data_list$n, ...).
  # Two-station's data_list has no $d/$n -- it uses n_days/n_obs instead --
  # so data_list_d/data_list_n always arrive here as NULL for a two-station
  # fit. That's only safe because format_mcmc_mat_nosplit() never actually
  # reads either parameter. If either is ever wired into the function body,
  # two-station's delegation breaks silently (NULL flowing into logic that
  # expects a real value) rather than erroring. This pins the invariant
  # directly rather than waiting to notice a silent break downstream.
  used_vars <- all.vars(body(format_mcmc_mat_nosplit))
  expect_false('data_list_d' %in% used_vars)
  expect_false('data_list_n' %in% used_vars)
})


manual_test4 <- function() {

  library(streamMetabolizer)
  library(testthat)
  library(dplyr)
  library(tidyr)
  library(gridExtra)
  library(ggplot2)
  source('tests/testthat/helper-rmse_DO.R')

  # test data
  dat <- data_metab('1', res='30')

  # create a list of all models to run
  opts <- expand.grid(
    type='bayes',
    pool_K600='none',#c('none','normal','linear','binned'),
    err_obs_iid=c(TRUE, FALSE),
    err_proc_acor=F,#c(FALSE, TRUE),
    err_proc_iid=c(FALSE, TRUE),
    ode_method=c('trapezoid','euler'),
    GPP_fun='linlight',
    ER_fun='constant',
    deficit_src=c('DO_mod','DO_obs'),
    engine='stan',
    check_validity=FALSE,
    stringsAsFactors=FALSE)
  stanfiles <- opts %>%
    rowwise %>% do(tibble::tibble(model_name=do.call(mm_name, .))) %>%
    unlist(use.names=FALSE) %>% sort %>% {.[!grepl('__', .)]} %>%
    .[. %in% mm_valid_names('bayes')]
  stanfiles

  mms <- lapply(setNames(nm=stanfiles[7]), function(sf) {
    message(sf)
    sp <- revise(specs(sf, burnin_steps=100, saved_steps=50), params_out=union(params_out, c('err_obs_iid','err_proc_iid')))
    metab(sp, dat)
  })

  bind_rows(lapply(mms, get_params))
  sapply(mms, get_fitting_time)

  bind_rows(lapply(mms, function(mm) {
    get_params(mm) %>%
      mutate(model=get_specs(mm)$model_name) %>%
      select(model, everything()) %>%
      bind_cols(select(get_fit(mm)$daily, ends_with('Rhat')))
  }))

  # if models were run on a cluster, pull and query them here. get in the right directory.
  mm_files <- dir(path='../metab_tests/explore/007_stan_tests/', pattern='16.*', full.names = TRUE)
  mm_run <- substring(gsub('\\.Rds', '.stan', basename(mm_files)), 8)
  mm_unrun <- sort(setdiff(stanfiles, mm_run))
  mms <- lapply(mm_files, readRDS)
  mm_failures <- mm_run[sapply(mms, function(mm) is.null(get_mcmc(mm)))]
  setdiff(mm_run, mm_failures)

  library(gridExtra)
  do.call(grid.arrange, c(lapply(mms[c(1,4,2,3)], function(mm)
    plot_DO_preds(mm) + ggtitle(mm@specs$model_name)), list(nrow=2, ncol=2)))
  do.call(grid.arrange, c(lapply(mms[c(1,4,2,3)], function(mm)
    plot_metab_preds(mm) + ggtitle(mm@specs$model_name)), list(nrow=2, ncol=2)))
  do.call(grid.arrange, c(lapply(mms[c(1,5,9,2,6,10,3,7,11,4,8,12)], function(mm)
    plot_DO_preds(mm) + ggtitle(mm@specs$model_name)), list(nrow=4, ncol=3)))
  do.call(grid.arrange, c(lapply(mms[c(1,5,9,2,6,10,3,7,11,4,8,12)], function(mm)
    traceplot(get_mcmc(mm), 'K600_daily_mu') + ggtitle(mm@specs$model_name)), list(nrow=4, ncol=3)))
}

manual_test1 <- function() {

  library(streamMetabolizer)
  library(testthat)
  library(dplyr)
  source('tests/testthat/helper-rmse_DO.R')

  test_that("lots of bayesian models available", {
    expect_lt(42, length(mm_valid_names('bayes')))
  })

  test_that("simple bayesian models run and implement the interface", {

    # simple test data (except for being light saturating)
    dat <- data_metab('1', res='30')

    # 1-core model
    mm <- mm_name('bayes', err_proc_acor=FALSE, err_proc_iid=FALSE) %>%
      specs(n_chains=1, n_cores=1, burnin_steps=300, saved_steps=100) %>%
      metab(data=dat)

    # run the model through its interface paces
    expect_equal(1, length(grep("SAMPLING FOR MODEL 'b_np_oi_tr_plrckm' NOW", get_log(mm)$MCMC_All_Days)))
    expect_lt(get_fitting_time(mm)['elapsed'], 120)
    expect_lt(rmse_DO(predict_DO(mm)), 0.2)
    plot_metab_preds(mm)
    plot_DO_preds(mm)
    traceplot(get_mcmc(mm), pars='GPP_daily')


    # 4-core model
    mm <- mm_name('bayes', err_proc_acor=FALSE, err_proc_iid=FALSE) %>%
      specs(n_chains=2, n_cores=4, burnin_steps=300, saved_steps=100) %>%
      metab(data=dat)

    # run the model through its interface paces
    expect_equal(2, length(grep("SAMPLING FOR MODEL 'b_np_oi_tr_plrckm' NOW", get_log(mm)$MCMC_All_Days)))
    expect_lt(get_fitting_time(mm)['elapsed'], 120)
    expect_lt(rmse_DO(predict_DO(mm)), 0.2)
    plot_metab_preds(mm)
    plot_DO_preds(mm)
    traceplot(get_mcmc(mm), pars='GPP_daily')

  })

  # sp() applies to next two test_that calls
  sp <- function(split_dates) { replace(
    specs(mm_name('bayes', err_proc_iid=FALSE),
          n_cores=3, n_chains=3, burnin_steps=300, saved_steps=200, verbose=FALSE),
    'split_dates', split_dates
  ) }

  test_that("error-free models can be run with split or combined dates", {
    dat <- data_metab('1', res='30')
    nosplit <- metab(sp(split_dates=FALSE), dat)
    split <- metab(sp(split_dates=TRUE), dat)
    # expect the same fitted parameter dimensions and similar estimates by either method
    expect_true(all(dim(get_params(nosplit)) == dim(get_params(split))))
    expect_true(max(abs(get_params(nosplit)[c('GPP.daily','ER.daily','K600.daily')] / get_params(split)[c('GPP.daily','ER.daily','K600.daily')] - 1)) < 0.2)

    dat <- data_metab('3', res='30')
    nosplit <- metab(sp(split_dates=FALSE), dat)
    split <- metab(sp(split_dates=TRUE), dat)
    # expect the same fitted parameter dimensions and similar estimates by either method
    expect_true(all(dim(get_params(nosplit)) == dim(get_params(split))))
    expect_true(max(abs(get_params(nosplit)[c('GPP.daily','ER.daily','K600.daily')] / get_params(split)[c('GPP.daily','ER.daily','K600.daily')] - 1)) < 0.2)

  })

  test_that("error and warning messages are printed with the mm object if present", {
    dat <- data_metab('1', res='30', flaws=c('missing start'))
    expect_warning(metab(sp(FALSE), dat), "Modeling failed: no valid days of data")
    expect_warning(metab(sp(TRUE), dat), "Modeling failed: no valid days of data")

    dat <- data_metab('3', res='30', flaws=c('missing middle'))
    expect_equal(get_params(metab(sp(FALSE), dat))$errors, c('','uneven timesteps',''))
    expect_equal(get_params(metab(sp(TRUE), dat))$errors, c('','uneven timesteps',''))
  })
}

manual_test2 <- function() {

  library(streamMetabolizer)
  library(testthat)
  library(dplyr)
  source('tests/testthat/helper-rmse_DO.R')

  # Make sure many combinations of models run
  dat <- data_metab('1','30')

  mm <- metab(revise(specs("b_np_oi_tr_plrckm.stan"), params_out=c(params_out, 'DO_mod')), dat)
  mm <- metab(revise(specs("b_np_oi_tr_plrcko.stan"), params_out=c(params_out, 'DO_mod')), dat)

}

manual_test3 <- function() {

  # as of 9/15/2016, the 'old' models are actually the newest, having been
  # rewritten since these tests

  library(streamMetabolizer)
  library(dplyr)
  # faster stan Kl_pcpi_ko model?
  dat <- mutate(data_metab('10', res='10'), discharge=3)
  sp <- specs("b_Kl_pcpi_tr_plrcko.stan", n_chains=3, n_cores=3, burnin_steps=300, saved_steps=100, keep_mcmcs=TRUE)
  mm_old <- metab(specs=sp, data=dat) # used to be wronger. 20 sec
  get_fitting_time(mm_old)
  plot_metab_preds(mm_old)
  plot_DO_preds(mm_old)
  pairs(get_mcmc(mm_old), pars=c("GPP_daily[7]", "ER_daily[7]", "K600_daily[7]"))
  sp <- specs("b_Kl_pcpi_tr_plrcko.stan", n_chains=3, n_cores=3, burnin_steps=300, saved_steps=100, keep_mcmcs=TRUE,
              K600_daily_beta_mu=c(intercept=1, slope=2.3), K600_daily_beta_sigma=c(intercept=0.3, slope=0.3)) %>%
    revise(model_name='inst/models/b_Kl_pcpi_pm_plrcko_sfs2loglog.stan',
           K600_daily_sigma_rate=2, err_proc_acor_phi_shape=1, err_proc_acor_phi_rate=50000, err_proc_acor_sigma_rate=0.001, err_proc_iid_sigma_rate=0.02,
           params_in=c('GPP_daily_mu','GPP_daily_sigma','ER_daily_mu','ER_daily_sigma','K600_daily_beta_mu','K600_daily_beta_sigma','K600_daily_sigma_rate',
                       'err_proc_acor_phi_shape','err_proc_acor_phi_rate','err_proc_acor_sigma_rate','err_proc_iid_sigma_rate'))
  mm_new <- metab(specs=sp, data=dat) # 1:54, 1:51 with default K600_daily_beta_mu and _sigma. 0:40 with better ones
  get_fitting_time(mm_new)
  predict_metab(mm_new)
  plot_metab_preds(mm_new)
  plot_DO_preds(mm_new)
  traceplot(get_mcmc(mm_new), pars=c("GPP_daily[2]", "ER_daily[2]", "K600_daily[2]", "GPP_daily[7]", "ER_daily[7]", "K600_daily[7]"))
  traceplot(get_mcmc(mm_new), pars=c("err_proc_acor_phi", "err_proc_acor_sigma", "err_proc_iid_sigma"))
  traceplot(get_mcmc(mm_new), pars=c("K600_daily_beta","K600_daily_sigma"))
  pairs(get_mcmc(mm_new), pars=c("GPP_daily[7]", "ER_daily[7]", "K600_daily[7]"))
  pairs(get_mcmc(mm_new), pars=c("GPP_daily[2]", "ER_daily[2]", "K600_daily[2]"))
  pairs(get_mcmc(mm_new), pars=c("err_proc_acor_phi", "err_proc_acor_sigma", "err_proc_iid_sigma"))

  # faster oipi_km model (state space)
  dat <- mutate(data_metab('10', res='10'), discharge=3)
  sp <- specs("b_Kl_oipi_tr_plrckm.stan", n_chains=3, n_cores=3, burnin_steps=300, saved_steps=100, keep_mcmcs=TRUE)
  sp <- sp %>% revise(params_out=sp$params_out[-which(sp$params_out == 'err_proc_iid')])
  mm_old <- metab(specs=sp, data=dat) # 0:18 but magnitudes are all off by about 10 or 15
  sp <- specs("b_Kl_oipi_tr_plrckm.stan", n_chains=3, n_cores=3, burnin_steps=300, saved_steps=100, keep_mcmcs=TRUE, keep_mcmc_data=TRUE,
              K600_daily_beta_mu=c(intercept=1, slope=2.3), K600_daily_beta_sigma=c(intercept=0.3, slope=0.3)) %>%
    revise(model_name='inst/models/b_Kl_oipi_pm_plrckm_sfs2loglog.stan',
           K600_daily_sigma_rate=0.2, err_obs_iid_sigma_rate=0.05, err_proc_iid_sigma_rate=0.01,
           params_in=c('GPP_daily_mu','GPP_daily_sigma','ER_daily_mu','ER_daily_sigma','K600_daily_beta_mu','K600_daily_beta_sigma','K600_daily_sigma_rate','err_obs_iid_sigma_rate','err_proc_iid_sigma_rate'),
           params_out=setdiff(params_out, 'err_proc_iid')) # a little faster f run without 'err_proc_iid' in params_out (but just 1 Mb)
  mm_new <- metab(specs=sp, data=dat) # 2:37
  predict_metab(mm_new)
  plot_metab_preds(mm_new)
  plot_DO_preds(mm_new)
  traceplot(get_mcmc(mm_new), pars=c("GPP_daily[2]", "ER_daily[2]", "K600_daily[2]", "GPP_daily[7]", "ER_daily[7]", "K600_daily[7]"))
  traceplot(get_mcmc(mm_new), pars=c("K600_daily_beta","K600_daily_sigma","err_obs_iid_sigma", "err_proc_iid_sigma"), inc_warmup=TRUE)
  pairs(get_mcmc(mm_new), pars=c("GPP_daily[7]", "ER_daily[7]", "K600_daily[7]"))
  pairs(get_mcmc(mm_new), pars=c("GPP_daily[2]", "ER_daily[2]", "K600_daily[2]"))
  pairs(get_mcmc(mm_new), pars=c("err_proc_iid_sigma", "err_obs_iid_sigma"))
  # now with bob's reworking
  sp <- specs("b_Kl_oipi_tr_plrckm.stan", n_chains=3, n_cores=3, burnin_steps=300, saved_steps=100, keep_mcmcs=TRUE, keep_mcmc_data=TRUE,
              K600_daily_beta_mu=c(intercept=1, slope=2.3), K600_daily_beta_sigma=c(intercept=0.3, slope=0.3)) %>%
    revise(model_name='inst/models/b_Kl_oipi_pm_plrckm_sfs3loglog.stan',
           K600_daily_sigma_rate=0.2, err_obs_iid_sigma_rate=0.05, err_proc_iid_sigma_rate=0.01,
           params_in=c('GPP_daily_mu','GPP_daily_sigma','ER_daily_mu','ER_daily_sigma','K600_daily_beta_mu','K600_daily_beta_sigma','K600_daily_sigma_rate','err_obs_iid_sigma_rate','err_proc_iid_sigma_rate'),
           params_out=setdiff(params_out, 'err_proc_iid'))
  mm_new2 <- metab(specs=sp, data=dat) # 17, 20, 30 seconds for 300/100 steps, or 1:21 min for 500/500
  get_log(mm_new2)
  predict_metab(mm_new2)
  plot_metab_preds(mm_new2)
  plot_DO_preds(mm_new2)
  traceplot(get_mcmc(mm_new2), pars=c("GPP_daily[2]", "ER_daily[2]", "K600_daily[2]", "GPP_daily[7]", "ER_daily[7]", "K600_daily[7]"), inc_warmup=F)
  traceplot(get_mcmc(mm_new2), pars=c("K600_daily_beta","K600_daily_sigma"), inc_warmup=F)
  pairs(get_mcmc(mm_new2), pars=c("GPP_daily[7]", "ER_daily[7]", "K600_daily[7]"))
  pairs(get_mcmc(mm_new2), pars=c("GPP_daily[2]", "ER_daily[2]", "K600_daily[2]"))
  pairs(get_mcmc(mm_new2), pars=c("err_obs_iid_sigma", "err_proc_iid_sigma"))

  # compare to process error model
  dat <- mutate(data_metab('10', res='10'), discharge=3)
  sp <- specs("b_Kl_pi_tr_plrcko.stan", n_chains=3, n_cores=3, burnin_steps=200, saved_steps=100, keep_mcmcs=TRUE)
  mm_old_cim <- metab(specs=sp, data=dat) # 11 seconds for 200/100
  get_fit(mm_old_cim)$daily %>% select(ends_with('Rhat'))
  predict_metab(mm_old_cim)
  plot_metab_preds(mm_old_cim)
  plot_DO_preds(mm_old_cim)
  traceplot(get_mcmc(mm_old_cim), pars=c("GPP_daily[2]", "ER_daily[2]", "K600_daily[2]", "GPP_daily[7]", "ER_daily[7]", "K600_daily[7]"))
  traceplot(get_mcmc(mm_old_cim), pars=c("err_proc_iid_sigma"))
  traceplot(get_mcmc(mm_old_cim), pars=c("K600_daily_beta","K600_daily_sigma"))
  pairs(get_mcmc(mm_old_cim), pars=c("GPP_daily[7]", "ER_daily[7]", "K600_daily[7]"))
  pairs(get_mcmc(mm_old_cim), pars=c("GPP_daily[2]", "ER_daily[2]", "K600_daily[2]"))
  sp <- specs("b_Kl_pi_pm_plrcko_sfs.stan", n_chains=3, n_cores=3, burnin_steps=200, saved_steps=100, keep_mcmcs=TRUE, keep_mcmc_data=FALSE,
              K600_daily_beta_mu=c(intercept=1, slope=2.3), K600_daily_beta_sigma=c(intercept=0.3, slope=0.3)) %>%
    revise(K600_daily_sigma_rate=1, err_proc_iid_sigma_rate=0.03,
           params_in=c('GPP_daily_mu','GPP_daily_sigma','ER_daily_mu','ER_daily_sigma','K600_daily_beta_mu','K600_daily_beta_sigma',
                       'K600_daily_sigma_rate','err_proc_iid_sigma_rate'),
           delete=c('K600_daily_sigma_scale','err_proc_iid_sigma_scale'))
  mm_new <- metab(specs=sp, data=dat) # 8-12 sec for 200/100
  plot_metab_preds(mm_new)
  plot_DO_preds(mm_new)
  traceplot(get_mcmc(mm_new), pars=c("GPP_daily[2]", "ER_daily[2]", "K600_daily[2]", "GPP_daily[7]", "ER_daily[7]", "K600_daily[7]"))

  # faster stan oi model
  dat <- data_metab('10', res='10')
  sp <- specs('b_np_oi_pm_plrckm.stan', n_chains=3, n_cores=3, burnin_steps=300, saved_steps=100, verbose=FALSE, keep_mcmcs=TRUE,
              GPP_daily_sigma=4, ER_daily_sigma=4, K600_daily_sigma=4)
  mm_slow <- metab(specs=sp, data=dat) # 9 sec
  # SFS version
  mm_fast <- sp %>%
    revise(err_obs_iid_sigma_rate=0.2, model_name='b_np_oi_pm_plrckm_faster.stan', delete=c('err_obs_iid_sigma_scale'),
           params_in=setdiff(c(params_in, 'err_obs_iid_sigma_rate'),c('err_obs_iid_sigma_scale'))) %>%
    metab(data=dat) # 51 sec w/ new compilation, 16-17 sec w/ err_obs_sigma_rate at 0.2 or 2 or 0.02
  mm <- mm_fast
  get_fitting_time(mm)
  plot_DO_preds(predict_DO(mm))
  plot_metab_preds(predict_metab(mm))
  select(get_fit(mm)$daily, date, ends_with('Rhat')) # 400 iterations is enough!
  traceplot(get_mcmc(mm), pars=c('GPP_daily[6]','ER_daily[6]','K600_daily[6]','err_obs_iid_sigma'), inc_warmup=TRUE)

  # faster stan Kl model?
  dat <- mutate(data_metab('10', res='10'), discharge=3)
  sp <- specs("b_Kl_oi_tr_plrckm.stan", n_chains=3, n_cores=3, burnin_steps=300, saved_steps=100, verbose=FALSE, keep_mcmcs=TRUE)
  mm_old <- metab(
    specs=revise(sp, model_name='inst/models/b_Kl_oi_tr_plrckm.stan'),
    data=dat) # 20 sec - baseline after june 2016 speed-ups
  sp2 <- sp %>% revise(
    params_in=c('GPP_daily_mu','GPP_daily_sigma','ER_daily_mu','ER_daily_sigma',
                'K600_daily_beta_mu','K600_daily_beta_sigma','K600_daily_sigma_rate','err_obs_iid_sigma_rate'),
    K600_daily_beta_mu=c(intercept=1, slope=2.3),
    K600_daily_beta_sigma=c(intercept=0.3, slope=0.3),
    K600_daily_sigma_rate=1,err_obs_iid_sigma_rate=0.1,
    delete=c("K600_daily_sigma_location","K600_daily_sigma_scale","err_obs_iid_sigma_scale"))
  sp3 <- sp2 %>% revise(
    K600_daily_beta_mu=c(intercept=10, slope=3),
    K600_daily_beta_sigma=c(intercept=8, slope=2))
  mm_new2linlog <- metab(revise(sp2, model_name='inst/models/b_Kl_oi_pm_plrckm_sfs2linlog.stan'), data=dat) # 18 sec
  mm_new2loglog <- metab(revise(sp2, model_name='inst/models/b_Kl_oi_pm_plrckm_sfs2loglog.stan'), data=dat) # 18 sec
  mm_new3linlog <- metab(revise(sp3, model_name='inst/models/b_Kl_oi_pm_plrckm_sfs2linlog.stan'), data=dat) # 21 sec
  mm_new3loglog <- metab(revise(sp3, model_name='inst/models/b_Kl_oi_pm_plrckm_sfs2loglog.stan'), data=dat) # 53 sec. good priors matter!
}


# takes too long to do all the time. also, saving and reloading a stan model
# doesn't work!
test_that("test that metab_models can be saved & reloaded (see helper-save_load.R)", {

  skip("NB: saving and reloading a stan model doesn't work!")

  source('tests/testthat/helper-save_load.R')

  # fit model
  dat <- data_metab('1', res='30')
  mmb <- mm_name('bayes', err_proc_acor=FALSE, err_proc_iid=FALSE, engine='stan') %>%
    specs(n_chains=1, n_cores=1, burnin_steps=300, saved_steps=100) %>%
    metab(data=dat)
  mm <- mmb

  # see if saveRDS with gzfile, compression=9 works well
  rdstimes <- save_load_timing(mm, reps=1)
  expect_true('gz6' %in% rdstimes$typelevel[1:3], info="gz6 is reasonably efficient for saveRDS")
  plot_save_load_timing(rdstimes)

  # save and load the mm, make sure it stays the same
  mmls <- test_save_load_recovery(mm) # fails!
  expect_equal(get_mcmc(mmls$original), get_mcmc(mmls$reloaded)) # fails!

})

## Code you can run after fitting any MCMC model
useful_code <- function() {

  # when things go bad
  debug_metab <- function(specs) {
    mcmc_args <- c('engine','model_path','params_out','n_chains','n_cores','burnin_steps','saved_steps','thin_steps','verbose')
    specs$model_path <- system.file(paste0("models/", specs$model_name), package="streamMetabolizer")
    data_list <- streamMetabolizer:::prepdata_bayes(
      data=vfrench1day, data_daily=NULL, ply_date="2012-08-24",
      specs=specs, engine=specs$engine, model_name=specs$model_name, priors=specs$priors)
    system.time({
      suppressWarnings(
        {mcmc_out <- do.call(streamMetabolizer:::mcmc_bayes, c(list(data_list=data_list, keep_mcmc=TRUE), specs[mcmc_args]))})
    })
    mcmc_out
  }
  mcmc <- debug_metab(specs)

  # when things didn't break
  expect_accurate <- function(mm) {
    mfile <- mm@specs$model_name
    expect_silent(metab <- predict_metab(mm))
    expect_silent(DO_preds <- predict_DO(mm))
    expect_lt(rmse_DO(DO_preds), 0.3) #, info=mfile)
  }
  expect_is(get_fitting_time(mm), "proc_time")
  expect_equal(names(mm@mcmc)[2], "2012-08-24")
  expect_accurate(mm)
  expect_equal(predict_metab(mm)$GPP.lower, get_fit(mm)$GPP_daily_2.5pct)

  # useful things to report
  plot_DO_preds(predict_DO(mm))

  ## Code you can run after fitting any Stan model
  rstan::traceplot(get_mcmc(mm)[[1]])
  rstan::traceplot(get_mcmc(mm)[[1]], inc_warmup=TRUE)
  expect_is(get_mcmc(mm)[[2]], "stanfit")
  get_fit(mm)[grep("Rhat", names(get_fit(mm)))]

}
