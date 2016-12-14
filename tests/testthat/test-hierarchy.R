context("hierarchy")

## There are currently no automated tests here because everything is slow.
## Longer tests, where we can investigate parameter estimates & convergence, are
## in manual_tests()

manual_tests2 <- function() {
  
  library(streamMetabolizer)
  library(testthat)
  library(dplyr)
  library(tidyr)
  library(gridExtra)
  library(ggplot2)
  source('tests/testthat/helper-rmse_DO.R')
  
  # test data
  library(mda.streams)
  site <- 'nwis_03259757'
  datall <- get_metab_data(
    site, disch=choose_data_source('disch','nwis_03259757'), model='metab_bayes') %>%
    filter(complete.cases(.))
  dat <- streamMetabolizer:::mm_filter_dates(datall, date_start='2014-04-03', date_end='2014-04-12')
  
  # create a list of all models to run
  opts <- expand.grid(
    type='bayes',
    pool_K600=c('none','normal','linear','binned'),
    err_obs_iid=T,#c(TRUE, FALSE),
    err_proc_acor=F,#c(FALSE, TRUE),
    err_proc_iid=F,#c(FALSE, TRUE),
    ode_method='trapezoid',#c('trapezoid','euler'),
    GPP_fun='linlight',
    ER_fun='constant',
    deficit_src='DO_mod',#c('DO_mod','DO_obs'),
    engine='stan',
    check_validity=FALSE,
    stringsAsFactors=FALSE)
  stanfiles <- opts %>% 
    rowwise %>% do(data_frame(model_name=do.call(mm_name, .))) %>% 
    unlist(use.names=FALSE) %>% sort %>% {.[!grepl('__', .)]} %>%
    .[. %in% mm_valid_names('bayes')]
  
  mms <- lapply(setNames(nm=stanfiles), function(sf) {
    message(sf)
    msp <- if(mm_parse_name(sf)$pool_K600 %in% c('none')) specs(sf) else 
      specs(sf, K600_daily_sdlog_scale=0.0001, err_obs_iid_sigma_scale=3)
    mdat <- if(mm_parse_name(sf)$pool_K600 %in% c('linear','binned')) dat else select(dat, -discharge)
    #metab(revise(msp, burnin_steps=10, saved_steps=10), mdat) # just to compile model
    metab(revise(msp, burnin_steps=200, saved_steps=200), mdat)
  })
  
  # see issue #291 for periodic output reports
  
  # fitting times
  sapply(mms, function(mm) get_fitting_time(mm)[['elapsed']])
  
  # K~Q relationship estimates
  kqfits <- lapply(mms, function(mm) {
    mname <- get_specs(mm)$model_name
    pars <- get_params(mm, unc='ci')
    fit <- get_fit(mm)
    min.K.lower = min(c(pars$K600.daily.lower[pars$K600.daily.lower>0], pars$K600.daily[pars$K600.daily>0]))/2
    dailymean <- function(data_ply, data_daily_ply, day_start, day_end, ply_date, ply_validity, timestep_days, ...) {
      data.frame(dischdaily = if(isTRUE(ply_validity[1])) mean(data_ply$discharge) else NA)
    }
    dischdaily <- mm_model_by_ply(model_fun=dailymean, data=dat, day_start=mms[[1]]@specs$day_start, day_end=mms[[1]]@specs$day_end)
    
    switch(
      mm_parse_name(mname)$pool_K600,
      'none'={
        ggplot(pars, aes(x=date, y=K600.daily)) + 
          geom_abline(intercept=0, color='darkgrey') +
          geom_point() + geom_errorbar(aes(ymin=pmax(0, K600.daily.lower), ymax=K600.daily.upper)) +
          scale_y_log10() +
          ylab('K600') + xlab('date')
      },
      'normal'={
        K <- list(
          mean=fit$overall$K600_daily_predlog_50pct,
          sd=fit$overall$K600_daily_sdlog_50pct)
        ggplot(pars, aes(x=log(K600.daily))) +
          geom_density(fill='lightgrey') +
          geom_rug(sides='t') +
          geom_point(x=K$mean, y=0, color='red') +
          geom_errorbarh(aes(x=K$mean, y=0, xmin=K$mean-K$sd, xmax=K$mean+K$sd), color='red', height=0.05) + 
          stat_function(fun=dnorm, args=list(mean=K$mean, sd=K$sd), color="red") +
          xlab('ln(K600)') + ylab('density')
      },
      'linear'={
        parslnQ <- mutate(
          pars, 
          lnQold=get_mcmc_data(mm)$lnQ_daily,
          lnQ=log(dischdaily$dischdaily),
          lnK.pred=fit$daily$K600_daily_predlog_50pct,
          lnK=log(K600.daily),
          lnK.lower=log(pmax(min.K.lower, K600.daily.lower)),
          lnK.upper=log(K600.daily.upper))
        K <- list(
          icpt=fit$overall$lnK600_lnQ_intercept_50pct,
          slope=fit$overall$lnK600_lnQ_slope_50pct,
          sd=fit$overall$K600_daily_sdlog_50pct)
        ggplot(parslnQ, aes(x=lnQ, y=lnK, ymin=lnK.lower, ymax=lnK.upper)) +
          geom_abline(intercept=K$icpt, slope=K$slope, col='red') +
          geom_ribbon(aes(ymin=lnK.pred-K$sd, ymax=lnK.pred+K$sd), color=NA, fill='red', alpha=0.2) +
          geom_point() + geom_errorbar() +
          ylab('ln(K600)') + xlab('ln(Q)')
      },
      'binned'={
        K <- data_frame(
          lnQ=get_specs(mm)$K600_lnQ_nodes_centers,
          lnKbin=fit[[as.character(length(lnQ))]]$lnK600_lnQ_nodes_50pct,
          lnKbin.lower=fit[[as.character(length(lnQ))]]$lnK600_lnQ_nodes_2.5pct,
          lnKbin.upper=fit[[as.character(length(lnQ))]]$lnK600_lnQ_nodes_97.5pct)
        lnQdat <- 
          data_frame(
            lnQ_bin1=get_mcmc_data(mm)$lnQ_bins[1,],
            lnQ_bin2=get_mcmc_data(mm)$lnQ_bins[2,],
            lnQ_bin_weights1=get_mcmc_data(mm)$lnQ_bin_weights[1,],
            lnQ_bin_weights2=get_mcmc_data(mm)$lnQ_bin_weights[2,]) %>%
          mutate(
            lnQ_daily=K$lnQ[lnQ_bin1]*lnQ_bin_weights1 + K$lnQ[lnQ_bin2]*lnQ_bin_weights2,
            lnK_pred=K$lnKbin[lnQ_bin1]*lnQ_bin_weights1 + K$lnKbin[lnQ_bin2]*lnQ_bin_weights2)
        parslnQ <- mutate(
          pars, 
          lnQ=log(dischdaily$dischdaily), #lnQdat$lnQ_daily truncates at last bins
          lnK.pred=fit$daily$K600_daily_predlog_50pct,
          lnK.eq=lnQdat$lnK_pred,
          lnK=log(K600.daily),
          lnK.lower=log(pmax(min.K.lower, K600.daily.lower)),
          lnK.upper=log(K600.daily.upper))
        ggplot(parslnQ, aes(x=lnQ)) +
          geom_line(data=K, aes(y=lnKbin), col='red') +
          geom_ribbon(data=K, aes(ymin=lnKbin.lower, ymax=lnKbin.upper), color=NA, fill='red', alpha=0.2) +
          geom_point(data=K, aes(y=lnKbin), col='red') +
          geom_point(aes(y=lnK)) + geom_errorbar(aes(ymin=lnK.lower, ymax=lnK.upper)) +
          ylab('ln(K600)') + xlab('ln(Q)')
      }
    ) + theme_classic() + ggtitle(mname)
  })
  do.call(grid.arrange, kqfits)
  
  # daily parameter estimates
  pars <- 
    bind_rows(lapply(mms, function(mm) {
      get_params(mm) %>%
        mutate(model=get_specs(mm)$model_name) %>%
        select(model, everything()) %>%
        bind_cols(select(get_fit(mm)$daily, ends_with('Rhat')))
    })) %>%
    select(model, date, GPP.daily, ER.daily, K600.daily) %>%
    arrange(date, model)
  pars
  pars %>%
    gather(param, estimate, GPP.daily, ER.daily, K600.daily) %>%
    mutate(param=ordered(param, levels=c('GPP.daily', 'ER.daily', 'K600.daily'))) %>%
    ggplot(aes(x=date, y=estimate, color=model)) + 
    geom_abline(intercept=0, slope=0, color='darkgrey') +
    geom_point() + geom_line() + theme_classic() + 
    facet_grid(param ~ ., scales='free_y')
  
  # daily metabolism preds
  preds <- 
    bind_rows(lapply(mms, function(mm) {
      predict_metab(mm) %>%
        mutate(model=get_specs(mm)$model_name) %>%
        select(model, everything())
    })) %>%
    select(model, date, GPP, ER) %>%
    arrange(date, model)
  preds
  preds %>%
    gather(pred, estimate, GPP, ER) %>%
    mutate(pred=ordered(pred, levels=c('GPP', 'ER'))) %>%
    ggplot(aes(x=date, y=estimate, color=model)) + 
    geom_abline(intercept=0, slope=0, color='darkgrey') +
    geom_point() + geom_line() + theme_classic() + 
    facet_grid(pred ~ ., scales='free_y')
  
  # DO preds
  do.call(grid.arrange, c(lapply(mms, function(mm) 
    plot_DO_preds(mm, y_var='pctsat') + ggtitle(mm@specs$model_name)), list(nrow=2, ncol=2)))
  
}

manual_tests <- function() {
  
  # devtools::install_github("stan-dev/shinystan")
  # library(shinystan)
  library(streamMetabolizer)
  library(dplyr)
  
  # pool_K600="none"
  dat <- data_metab('1', res='30')
  oi <- metab(specs(mm_name('bayes', err_obs_iid=TRUE, err_proc_iid=FALSE, pool_K600='none', ode_method='euler')), dat)
  plot_distribs(oi, 'err_obs_iid_sigma')
  stan_trace(get_mcmc(oi), pars='err_obs_iid_sigma')
  stan_hist(get_mcmc(oi), pars='err_obs_iid_sigma')
  pi <- metab(specs(mm_name('bayes', err_obs_iid=FALSE, err_proc_iid=TRUE, pool_K600='none', deficit_src='DO_obs'), err_proc_iid_sigma_scale=0.001), dat) # DO_obs currently required for pi models
  plot_distribs(pi, 'err_proc_iid_sigma')
  stan_trace(get_mcmc(pi), pars='err_proc_iid_sigma')
  stan_hist(get_mcmc(pi), pars='err_proc_iid_sigma')
  oipi <- metab(specs(mm_name('bayes', err_obs_iid=TRUE, err_proc_iid=TRUE, pool_K600='none')), dat)
  plot_distribs(oipi, 'err_obs_iid_sigma')
  plot_distribs(oipi, 'err_proc_iid_sigma')
  stan_trace(get_mcmc(oipi), pars=c('err_obs_iid_sigma','err_proc_iid_sigma'))
  stan_hist(get_mcmc(oipi), pars=c('err_obs_iid_sigma','err_proc_iid_sigma'))
  
  lapply(list(oi,pi,oipi), get_fitting_time)
  lapply(list(oi,pi,oipi), get_params)
  lapply(list(oi,pi,oipi), predict_metab)
  lapply(list(oi,pi,oipi), plot_DO_preds)
  
  # pool_K600="normal"
  dat <- data_metab('10', res='30')
  sp <- specs(mm_name('bayes', pool_K600='normal'))
  
  # pool_K600="linear"
  dat <- data_metab('10', res='30')
  sp <- specs(mm_name('bayes', pool_K600='linear'))
  
  # pool_K600="binned"
  dat <- data_metab('10', res='30')
  sp <- specs(mm_name('bayes', pool_K600='binned'))
  
  
  
  #### normal hierarchical model ####
  x <- seq(0,1,by=0.01)
  plot(x=x, y=dgamma(x, shape=1, rate=10), type='l', ylim=c(0,4))
  qgamma(p=c(0.001,0.9), shape=1, rate=100)
  sp <- mm_name('bayes', pool_K600='normal') %>%
    specs(K600_daily_mu_mu=30, K600_daily_mu_sigma=30,
          K600_daily_sigma_shape=1, K600_daily_sigma_rate=0.5,
          err_obs_iid_sigma_shape=1, err_obs_iid_sigma_rate=1,
          err_proc_iid_sigma_shape=1, err_proc_iid_sigma_rate=1,
          n_chains=1, n_cores=4, burnin_steps=1500, saved_steps=1500)
  mma <- metab(sp, data=dat)
  get_fitting_time(mma) # 877 sec
  #mms <- replace(sp, 'split_dates', TRUE) %>% metab(data=dat)
  #get_fitting_time(mms) # 73 sec
  
  
  mn <- mm_name("bayes", err_proc_acor=FALSE, err_proc_iid=FALSE, ode_method="trapezoid", deficit_src="DO_mod", engine="stan")
  ms <- specs(mn, burnin_steps=300, saved_steps=200, split_dates=FALSE)
  fitone <- metab_bayes(specs=replace(ms, 'split_dates', TRUE), vfrenchmedium) # one day at a time
  get_fitting_time(fitone) # 126 sec
  fitall <- metab_bayes(vfrenchmedium, specs=replace(ms, 'split_dates', FALSE)) # all three days at once
  get_fitting_time(fitall) # 87 sec
  
  #### medium-hard model ####
  mn <- mm_name("bayes", err_proc_acor=FALSE, err_proc_iid=TRUE, ode_method="trapezoid", deficit_src="DO_mod", engine="stan")
  ms <- specs(mn, burnin_steps=1000, saved_steps=300, n_chains=3, n_cores=4)
  # takes a really really long time! but maybe less now that all the sigmas have
  # gamma distributions? was 2.76 hours for 3 days worth of data (burning=3000,
  # saved=1000)! also didn't converge (1.9 < Rhat < 18.3 for all params)
  fitall <- metab_bayes(vfrenchmedium, specs=replace(ms, 'split_dates', FALSE))
  get_fitting_time(fitall)
  # after switching to gammas for the sigmas, 2149 secs (0.6 hr) for 1000 
  # burnin, 300 saved, 3 days. Rhats 2.1 to 39. dunno if this is better. GPP,
  # ER, and K are still generally too close to 0, though more varied (after 1300
  # iters)
  
  #### really hard model ####
  ms <- mm_name("bayes", pool_K600="none", err_proc_acor=TRUE, ode_method='Euler', engine="stan")
  ms <- specs(ms, n_chains=1, n_cores=1, burnin_steps=300, saved_steps=100, verbose=FALSE)
  system.time(fitall <- metab_bayes(vfrenchmedium, specs=replace(ms, 'split_dates', FALSE)))
  get_fitting_time(fitall) # 
  
  #### diagnostics code to run on any model ####
  get_fitting_time(fitall)
  fit <- get_mcmc(fitall)
  traceplot(fit, pars=c('GPP_daily[1]','ER_daily[1]','K600_daily[1]','err_obs_iid_sigma','err_proc_iid_sigma','lp__'))
  pairs(fit)
  predict_metab(fitall)
  plot_DO_preds(predict_DO(fitall))
  
}
