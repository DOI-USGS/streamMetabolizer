context("metab_bayes")

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
    
    # 1-day model in jags
    mmj <- mm_name('bayes', err_proc_acor=FALSE, err_proc_iid=FALSE, engine='jags') %>%
      specs(n_chains=1, n_cores=1, adapt_steps=250, burnin_steps=250, saved_steps=1000) %>%
      metab(data=dat)
    expect_lt(get_fitting_time(mmj)['elapsed'], 30)
    expect_lt(rmse_DO(predict_DO(mmj)), 0.2) #, info='jags')
    # plot_DO_preds(predict_DO(mmj))
    
    # jags & stan models w/ same options should reach very similar results
    expect_lt(sqrt(mean((predict_DO(mmb)$DO.mod - predict_DO(mmj)$DO.mod)^2)), 0.1) #, info='stan vs jags')
    
  })
}

# takes too long to do all the time. also, saving and reloading a stan model
# doesn't work! (it does seem to work for jags)
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
    jags_args <- c('engine','model_path','params_out','n_chains','n_cores','adapt_steps','burnin_steps','saved_steps','thin_steps','verbose')
    mcmc_args <- if(specs$engine == 'jags') jags_args else jags_args[-which(jags_args=='adapt_steps')]
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
  runjags::failed.jags()
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
  
  ## Code for JAGS models
  library(coda)
  par(mar=c(2,2,2,0.2))
  plot(as.mcmc.list(get_mcmc(mm)[[2]]), density=FALSE)
  expect_is(get_mcmc(mm)[[2]], "runjags")
  get_fit(mm)[grep("psrf", names(get_fit(mm)))]
  
  ## Code you can run after fitting any Stan model
  rstan::traceplot(get_mcmc(mm)[[1]])
  rstan::traceplot(get_mcmc(mm)[[1]], inc_warmup=TRUE)
  expect_is(get_mcmc(mm)[[2]], "stanfit")
  get_fit(mm)[grep("Rhat", names(get_fit(mm)))]
  
}
