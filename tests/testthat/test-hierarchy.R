context("hierarchy")

## There are currently no automated tests here because everything is slow.
## Longer tests, where we can investigate parameter estimates & convergence, are
## in manual_tests()

manual_tests <- function() {
  
  # devtools::install_github("stan-dev/shinystan")
  # library(shinystan)
  library(dplyr)
  
  #### data prep ####
  dat <- data_metab('10', res='10')
  
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
  
  
  mn <- mm_name("bayes", err_proc_acor=FALSE, err_proc_iid=FALSE, ode_method="pairmeans", deficit_src="DO_mod", engine="stan")
  ms <- specs(mn, burnin_steps=300, saved_steps=200, split_dates=FALSE)
  fitone <- metab_bayes(specs=replace(ms, 'split_dates', TRUE), vfrenchmedium) # one day at a time
  get_fitting_time(fitone) # 126 sec
  fitall <- metab_bayes(vfrenchmedium, specs=replace(ms, 'split_dates', FALSE)) # all three days at once
  get_fitting_time(fitall) # 87 sec
  
  #### super simple jags ####
  # mn <- mm_name("bayes", err_proc_acor=FALSE, err_proc_iid=FALSE, ode_method="pairmeans", deficit_src="DO_mod", engine="jags")
  # ms <- specs(mn, burnin_steps=300, saved_steps=200, split_dates=FALSE)
  # fitall <- metab_bayes(vfrenchmedium, specs=ms) # syntax error on line 3
  
  #### medium-hard model ####
  mn <- mm_name("bayes", err_proc_acor=FALSE, err_proc_iid=TRUE, ode_method="pairmeans", deficit_src="DO_mod", engine="stan")
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
  ms <- specs(ms, n_chains=1, n_cores=1, burnin_steps=300, saved_steps=100, verbose=TRUE)
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