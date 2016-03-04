context("hierarchy")

## There are currently no automated tests here because everything is slow.
## Longer tests, where we can investigate parameter estimates & convergence, are
## in manual_tests()

manual_tests <- function() {
  
  # devtools::install_github("stan-dev/shinystan")
  # library(shinystan)
  
  #### data prep ####
  vfrench <- streamMetabolizer:::load_french_creek(attach.units=FALSE)
  vfrenchmedium <- vfrench[vfrench$solar.time >= as.POSIXct("2012-09-17 21:00:00", tz="UTC") & 
                             vfrench$solar.time <= as.POSIXct("2012-09-22 06:00:00", tz="UTC"), ]
  # vfrenchshort <- vfrench[vfrench$solar.time >= as.POSIXct("2012-08-23 00:00:00", tz="UTC") & 
  #                           vfrench$solar.time <= as.POSIXct("2012-08-26 00:00:00", tz="UTC"), ]
  # vspring <- streamMetabolizer:::load_spring_creek(attach.units=FALSE)
  
  #### super simple model by split_dates=c(T,F) ####
  mn <- mm_name("bayes", err_proc_acor=FALSE, err_proc_iid=FALSE, ode_method="pairmeans", deficit_src="DO_mod", engine="stan")
  ms <- specs(mn, burnin_steps=300, saved_steps=200, split_dates=FALSE)
  fitone <- metab_bayes(vfrenchmedium, model_specs=replace(ms, 'split_dates', TRUE)) # one day at a time
  get_fitting_time(fitone) # 126 sec
  fitall <- metab_bayes(vfrenchmedium, model_specs=replace(ms, 'split_dates', FALSE)) # all three days at once
  get_fitting_time(fitall) # 87 sec
  
  #### super simple jags ####
  # mn <- mm_name("bayes", err_proc_acor=FALSE, err_proc_iid=FALSE, ode_method="pairmeans", deficit_src="DO_mod", engine="jags")
  # ms <- specs(mn, burnin_steps=300, saved_steps=200, split_dates=FALSE)
  # fitall <- metab_bayes(vfrenchmedium, model_specs=ms) # syntax error on line 3
  
  #### medium-hard model ####
  mn <- mm_name("bayes", err_proc_acor=FALSE, err_proc_iid=TRUE, ode_method="pairmeans", deficit_src="DO_mod", engine="stan")
  ms <- specs(mn, burnin_steps=1000, saved_steps=100, err_proc_iid_sigma_max=0.00001, n_chains=1)
  # takes a really really long time!
  fitall <- metab_bayes(vfrenchmedium, model_specs=replace(ms, 'split_dates', FALSE))
  # user  system elapsed 
  # 10.19    2.81 9955.42
  get_fitting_time(fitall) # 2.76 hours for 3 days worth of data!
  # user  system elapsed 
  # 9.78    2.81 9955.01 
  # also doesn't converge (1.9 < Rhat < 18.3 for all params)
  
  #### really hard model ####
  ms <- mm_name("bayes", pool_K600="none", err_proc_acor=TRUE, ode_method='Euler', engine="stan")
  ms <- specs(ms, err_proc_acor_phi_max=0.9, # err_proc_acor_sigma_max=1, err_proc_iid_sigma_max=0.1,
              n_chains=1, n_cores=1, burnin_steps=3000, saved_steps=1000, verbose=TRUE)
  system.time(fitall <- metab_bayes(vfrenchmedium, model_specs=replace(ms, 'split_dates', FALSE)))
  get_fitting_time(fitall) # 
  
  #### diagnostics code to run on any model ####
  get_fitting_time(fitall)
  fit <- get_mcmc(fitall)
  traceplot(fit, pars=c('GPP_daily[1]','ER_daily[1]','K600_daily[1]','err_obs_iid_sigma','err_proc_iid_sigma','lp__'))
  pairs(fit)
  predict_metab(fitall)
  plot_DO_preds(predict_DO(fitall))
  
}