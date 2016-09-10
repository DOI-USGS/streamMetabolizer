context("metab_bayes")

manual_test4 <- function() {
  # these three lines work within testthat() calls:
  skip_on_cran()
  skip_on_travis()
  skip_on_appveyor()
  skip_if_not_installed('deSolve')
  
  library(streamMetabolizer)
  library(dplyr)
  dat <- mutate(data_metab('3', res='30'), discharge=3)
  sp <- specs("b_Kl_oipi_eu_plrckm.stan", n_chains=3, n_cores=3, burnin_steps=200, saved_steps=100, verbose=TRUE, keep_mcmcs=TRUE)
  mm <- metab(specs=sp, data=dat)
  show_log(mm)
  plot_metab_preds(mm)
  plot_DO_preds(mm)
  traceplot(get_mcmc(mm), pars='GPP_daily')
}

manual_test5 <- function() {
  test_that("error-free models can be run with split or combined dates", {
    # this test doesn't actually do any test-that checks
    sp <- function(split_dates, engine) { replace(
      specs(mm_name('bayes', err_proc_iid=FALSE, engine=engine),
            n_cores=3, n_chains=3, burnin_steps=300, saved_steps=100, verbose=FALSE),
      'split_dates', split_dates
    ) }
    dat <- data_metab('1', res='30')
    # new compilation of any Stan model gives deprecation warnings as of 
    # 7/12/16; THESE ARE OKAY. subsequent runs of the compiled Stan model are 
    # quieter
    mm <- metab(sp(FALSE,'stan'), dat)
    expect_true(grepl("SAMPLING FOR MODEL 'b_np_oi_tr_plrckm' NOW", get_log(mm)))
    mm <- metab(sp(TRUE, 'stan'), dat)
    dat <- data_metab('3', res='30')
    mm <- metab(sp(FALSE,'stan'), dat)
    mm <- metab(sp(TRUE, 'stan'), dat)
    
    # ran into an interesting issue - shouldn't be bayes specific, but this is where I found it.
    # the error message (occurs during predict_DO, on about the 3rd call to metab()):
    # > mm <- metab(sp(FALSE,'stan'), dat)
    # Error in .Call("call_rkFixed", as.double(y), as.double(times), Func, Initfunc,  : 
    #   "call_rkFixed" not resolved from current namespace (deSolve)
    # Error in .C("unlock_solver") : 
    #   "unlock_solver" not resolved from current namespace (deSolve)
    
    # the solution may be to call library(deSolve) before loading/running
    # streamMetabolizer (see
    # http://r.789695.n4.nabble.com/deSolve-unresolved-namespace-error-td4679888.html)
    
  })
}


manual_test6 <- function() {
  test_that("error and warning messages are printed with the mm object if present", {
    # this test doesn't actually do any test-that checks
    sp <- function(split_dates, engine) { replace(
      specs(mm_name('bayes', err_proc_iid=FALSE, engine=engine),
            n_cores=3, n_chains=3, burnin_steps=300, saved_steps=100, verbose=FALSE),
      'split_dates', split_dates
    ) }
    dat <- data_metab('1', res='30', flaws=c('missing start'))
    mm <- metab(sp(FALSE,'stan'), dat) # Modeling failed: argument is of length zero
    mm <- metab(sp(TRUE, 'stan'), dat) # Modeling failed: no valid days of data
    dat <- data_metab('3', res='30', flaws=c('missing middle'))
    mm <- metab(sp(FALSE,'stan'), dat) # no problem
    mm <- metab(sp(TRUE, 'stan'), dat) # no problem
  })
}

manual_test3 <- function() {
  library(streamMetabolizer)
  library(dplyr)
  # faster stan Kl_pcpi_ko model?
  dat <- mutate(data_metab('10', res='10'), discharge=3)
  sp <- specs("b_Kl_pcpi_tr_plrcko.stan", n_chains=3, n_cores=3, burnin_steps=300, saved_steps=100, verbose=TRUE, keep_mcmcs=TRUE)
  mm_old <- metab(specs=sp, data=dat) # 12 sec, but hugely wrong
  plot_metab_preds(mm_old)
  plot_DO_preds(mm_old)
  sp <- specs("b_Kl_pcpi_tr_plrcko.stan", n_chains=3, n_cores=3, burnin_steps=300, saved_steps=100, verbose=TRUE, keep_mcmcs=TRUE,
              err_proc_acor_phi_shape=1, err_proc_acor_phi_rate=50000, err_proc_acor_sigma_rate=0.001, err_proc_iid_sigma_rate=0.02,
              K600_daily_beta_mu=c(intercept=1, slope=2.3), K600_daily_beta_sigma=c(intercept=0.3, slope=0.3)) %>%
    replace('model_name', 'inst/models/b_Kl_pcpi_tr_plrcko_sfs2loglog.stan')
  mm_new <- metab(specs=sp, data=dat) # 1:54, 1:51 with default K600_daily_beta_mu and _sigma. 0:40 with better ones
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
  sp <- specs("b_Kl_oipi_tr_plrckm.stan", n_chains=3, n_cores=3, burnin_steps=300, saved_steps=100, verbose=TRUE, keep_mcmcs=TRUE)
  sp <- sp %>% replace('params_out', list(sp$params_out[-which(sp$params_out == 'err_proc_iid')]))
  mm_old <- metab(specs=sp, data=dat) # 3:47 with compile, not much convergence; 4:33, 4:02 thereafter. this is a lin-log K model
  sp <- specs("b_Kl_oipi_tr_plrckm.stan", n_chains=3, n_cores=3, burnin_steps=300, saved_steps=100, verbose=TRUE, keep_mcmcs=TRUE, keep_mcmc_data=TRUE,
              K600_daily_sigma_rate=0.2, err_obs_iid_sigma_rate=0.05, err_proc_iid_sigma_rate=0.01,
              K600_daily_beta_mu=c(intercept=1, slope=2.3), K600_daily_beta_sigma=c(intercept=0.3, slope=0.3)) %>%
    replace('model_name', 'inst/models/b_Kl_oipi_tr_plrckm_sfs2loglog.stan')
  sp <- sp %>% replace('params_out', list(sp$params_out[-which(sp$params_out == 'err_proc_iid')]))
  mm_new <- metab(specs=sp, data=dat) # 5:28 with compile, 4:50 thereafter (14.7 Mb), 4:42 if run without 'err_proc_iid' in params_out (but just 1 Mb)
  predict_metab(mm_new)
  plot_metab_preds(mm_new)
  plot_DO_preds(mm_new)
  traceplot(get_mcmc(mm_new), pars=c("GPP_daily[2]", "ER_daily[2]", "K600_daily[2]", "GPP_daily[7]", "ER_daily[7]", "K600_daily[7]"))
  traceplot(get_mcmc(mm_new), pars=c("K600_daily_beta","K600_daily_sigma","err_obs_iid_sigma", "err_proc_iid_sigma"), inc_warmup=TRUE)
  pairs(get_mcmc(mm_new), pars=c("GPP_daily[7]", "ER_daily[7]", "K600_daily[7]"))
  pairs(get_mcmc(mm_new), pars=c("GPP_daily[2]", "ER_daily[2]", "K600_daily[2]"))
  pairs(get_mcmc(mm_new), pars=c("err_proc_acor_phi", "err_proc_acor_sigma", "err_proc_iid_sigma"))
  # now with bob's reworking
  sp <- specs("b_Kl_oipi_tr_plrckm.stan", n_chains=3, n_cores=3, burnin_steps=500, saved_steps=500, verbose=TRUE, keep_mcmcs=TRUE, keep_mcmc_data=TRUE,
              K600_daily_sigma_rate=0.2, err_obs_iid_sigma_rate=0.05, err_proc_iid_sigma_rate=0.01,
              K600_daily_beta_mu=c(intercept=1, slope=2.3), K600_daily_beta_sigma=c(intercept=0.3, slope=0.3)) %>%
    replace('model_name', 'inst/models/b_Kl_oipi_tr_plrckm_sfs3loglog.stan')
  sp <- sp %>% replace('params_out', list(sp$params_out[-which(sp$params_out == 'err_proc_iid')]))
  mm_new2 <- metab(specs=sp, data=dat) # 17, 20, 30 seconds for 300/100 steps, or 1:21 min for 500/500
  show_log(mm_new2)
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
  sp <- specs("b_Kl_pi_tr_plrcko.stan", n_chains=3, n_cores=3, burnin_steps=200, saved_steps=100, verbose=TRUE, keep_mcmcs=TRUE)
  mm_old_cim <- metab(specs=sp, data=dat) # 1:12 with compile for 300/300
  get_fit(mm_old_cim)$daily %>% select(ends_with('Rhat'))
  predict_metab(mm_old_cim)
  plot_metab_preds(mm_old_cim)
  plot_DO_preds(mm_old_cim)
  traceplot(get_mcmc(mm_old_cim), pars=c("GPP_daily[2]", "ER_daily[2]", "K600_daily[2]", "GPP_daily[7]", "ER_daily[7]", "K600_daily[7]"))
  traceplot(get_mcmc(mm_old_cim), pars=c("err_proc_iid_sigma"))
  traceplot(get_mcmc(mm_old_cim), pars=c("K600_daily_beta","K600_daily_sigma"))
  pairs(get_mcmc(mm_old_cim), pars=c("GPP_daily[7]", "ER_daily[7]", "K600_daily[7]"))
  pairs(get_mcmc(mm_old_cim), pars=c("GPP_daily[2]", "ER_daily[2]", "K600_daily[2]"))
  sp <- specs(
    "b_Kl_pi_tr_plrcko_sfs.stan", n_chains=3, n_cores=3, burnin_steps=200, saved_steps=100, verbose=TRUE, keep_mcmcs=TRUE, keep_mcmc_data=FALSE,
    K600_daily_sigma_rate=1, err_proc_iid_sigma_rate=0.03,
    K600_daily_beta_mu=c(intercept=1, slope=2.3), K600_daily_beta_sigma=c(intercept=0.3, slope=0.3))
  mm_new <- metab(specs=sp, data=dat) # 8-12 sec for 200/100
  plot_metab_preds(mm_new)
  plot_DO_preds(mm_new)
  traceplot(get_mcmc(mm_new), pars=c("GPP_daily[2]", "ER_daily[2]", "K600_daily[2]", "GPP_daily[7]", "ER_daily[7]", "K600_daily[7]"))
  
  # faster stan oi model
  dat <- data_metab('10', res='10')
  mmb <- mm_name('bayes', err_proc_acor=FALSE, err_proc_iid=FALSE, engine='stan')
  sp <- specs(mmb, n_chains=3, n_cores=3, burnin_steps=300, saved_steps=100, verbose=TRUE, keep_mcmcs=TRUE, 
              GPP_daily_sigma=4, ER_daily_sigma=4, K600_daily_sigma=4)
  #GPP_daily_sigma=10, ER_daily_sigma=10, K600_daily_sigma=10)
  mm_slow <- metab(specs=replace(sp, 'err_obs_iid_sigma_rate', 10), data=dat) # 75 sec w/ new compilation, 36-42 sec w/ daily sigmas at 4, 38-40 sec w/ daily sigmas at 10, 48 w/ err_obs_iid_sigma_rate=1 or 1000, 35 w/ err_obs_iid_sigma_rate=100 or 10
  sp$err_obs_iid_sigma_rate <- 0.2 # need to adjust because we're using the 'rate' as a lognormal scaling parameter now
  sp$model_name <- 'inst/models/b_np_oi_tr_plrckm_faster.stan'
  mm_fast <- metab(specs=replace(sp, 'err_obs_iid_sigma_rate', 0.02), data=dat) # 51 sec w/ new compilation, 16-17 sec w/ err_obs_sigma_rate at 0.2 or 2 or 0.02
  
  mm <- mm_fast
  get_fitting_time(mm)
  plot_DO_preds(predict_DO(mm))
  plot_metab_preds(predict_metab(mm))
  select(get_fit(mm)$daily, date, ends_with('Rhat')) # 400 iterations is enough!
  traceplot(get_mcmc(mm), pars=c('GPP_daily','ER_daily','K600_daily','err_obs_iid_sigma'), inc_warmup=TRUE)
  
  # faster stan Kl model?
  dat <- mutate(data_metab('10', res='10'), discharge=3)
  sp <- specs("b_Kl_oi_tr_plrckm.stan", n_chains=3, n_cores=3, burnin_steps=300, saved_steps=100, verbose=TRUE, keep_mcmcs=TRUE)
  mm_old <- metab(
    specs=replace(sp, 'model_name', 'inst/models/b_Kl_oi_tr_plrckm.stan'), 
    data=dat) # 48,41,53 sec - baseline
  mm_new <- metab(
    specs=sp %>%
      replace('model_name', 'inst/models/b_Kl_oi_tr_plrckm_sfs.stan') %>%
      replace('K600_daily_beta_mu', list(c(intercept=1, slope=2.3))) %>%
      replace('K600_daily_beta_sigma', list(c(intercept=0.3, slope=0.3))), 
    data=dat) # 44,31,40,36 sec - faster than no scaling but by precious little
  mm_newb <- metab(
    specs=sp %>%
      replace('model_name', 'inst/models/b_Kl_oi_tr_plrckm_sfs.stan') %>%
      replace('K600_daily_beta_mu', list(c(intercept=10, slope=3))) %>%
      replace('K600_daily_beta_sigma', list(c(intercept=8, slope=2))), 
    data=dat) # 53,49,52,41 sec - slightly slower than above b/c of parameters
  mm_loglog <- metab(
    specs=sp %>%
      replace('model_name', 'inst/models/b_Kl_oi_tr_plrckm_sfsloglog.stan') %>%
      replace('K600_daily_beta_mu', list(c(intercept=1, slope=2.3))) %>%
      replace('K600_daily_beta_sigma', list(c(intercept=0.3, slope=0.3))), 
    data=dat) # 1:19, 1:22, 1:32 sec - log-log is definitely slow
  mm_new2 <- metab(
    specs=sp %>%
      replace('model_name', 'inst/models/b_Kl_oi_tr_plrckm_sfs2.stan') %>%
      replace('K600_daily_beta_mu', list(c(intercept=10, slope=3))) %>%
      replace('K600_daily_beta_sigma', list(c(intercept=8, slope=2))), 
    data=dat) # 22,21 sec. where have you been all this time?
  mm_new2loglog <- metab(
    specs=sp %>%
      replace('model_name', 'inst/models/b_Kl_oi_tr_plrckm_sfs2loglog.stan') %>%
      replace('K600_daily_beta_mu', list(c(intercept=1, slope=2.3))) %>%
      replace('K600_daily_beta_sigma', list(c(intercept=0.3, slope=0.3))),
    data=dat) # 18,20 sec. so log-log isn't slow after all??
  mm_new2loglogb <- metab(
    specs=sp %>%
      replace('model_name', 'inst/models/b_Kl_oi_tr_plrckm_sfs2loglog.stan') %>%
      replace('K600_daily_beta_mu', list(c(intercept=10, slope=3))) %>%
      replace('K600_daily_beta_sigma', list(c(intercept=8, slope=2))),
    data=dat) # 45 sec. so good priors do matter
  
  # stan pi model
  dat <- data_metab('10', res='30')
  mm <- metab(specs=sp, data=dat)
  # elapsed time with no compilation: 85.39;  after pre-compilation: 61.42
  get_fitting_time(mm)
  plot_DO_preds(predict_DO(mm))
  plot_metab_preds(predict_metab(mm))
  traceplot(get_mcmc(mm), pars=c('GPP_daily[7]','ER_daily[7]','K600_daily[7]','err_obs_iid_sigma'), inc_warmup=TRUE)
  traceplot(get_mcmc(mm), pars=c('GPP_daily[7]','ER_daily[7]','K600_daily[7]','err_obs_iid_sigma'), inc_warmup=FALSE)
}

# The tests below cannot touch on all possible Bayesian models; that's
# work for a computing cluster. Instead we'll just inspect key features and a
# couple of simple models.

# even these simple tests take too long to do all the time.
manual_test1 <- function() {
  test_that("simple bayesian models run correctly", {
    
    # lots of bayesian models available
    expect_lt(42, length(mm_valid_names('bayes')))
    
    # get simplest possible data
    dat <- data_metab('1', res='30')
    
    # 1-day model in stan
    mmb <- mm_name('bayes', err_proc_acor=FALSE, err_proc_iid=FALSE, engine='stan') %>%
      specs(n_chains=1, n_cores=1, burnin_steps=300, saved_steps=100) %>%
      metab(data=dat)
    expect_lt(get_fitting_time(mmb)['elapsed'], 120)
    expect_lt(rmse_DO(predict_DO(mmb)), 0.2) #, info='stan')
    # plot_DO_preds(predict_DO(mmb))
    
  })
}

# takes too long to do all the time. also, saving and reloading a stan model
# doesn't work!
manual_test2 <- function() {
  testthat("test that metab_models can be saved & reloaded (see helper-save_load.R)", {
    
    # fit model
    dat <- data_metab('1', res='30')
    mmb <- mm_name('bayes', err_proc_acor=FALSE, err_proc_iid=FALSE, engine='stan') %>%
      specs(n_chains=1, n_cores=1, burnin_steps=300, saved_steps=100) %>%
      metab(data=dat)
    mm <- mmb
    
    # see if saveRDS with gzfile, compression=9 works well
    rdstimes <- save_load_timing(mm, reps=1) # autoloaded b/c script begins with 'helper' and is in this directory
    expect_true('gz6' %in% rdstimes$typelevel[1:3], info="gz6 is reasonably efficient for saveRDS")
    plot_save_load_timing(rdstimes)
    
    # save and load the mm, make sure it stays the same
    mmls <- test_save_load_recovery(mm) # fails!
    expect_equal(get_mcmc(mmls$original), get_mcmc(mmls$reloaded)) # fails!
    
  })
}

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
