context("metab_bayes")

test_that("bayes metab_spec functions work", {

  # ?specs_bayes should give useful help
  expect_is(specs_bayes_jags_nopool_oi(), "list")
  expect_is(specs_bayes_jags_nopool_oipc(), "list")
  
})

# get, format, & subset data
vfrench <- streamMetabolizer:::load_french_creek(attach.units=FALSE)
vfrenchshort <- vfrench[vfrench$local.time >= as.POSIXct("2012-08-23 00:00:00", tz="Etc/GMT+7") & 
                          vfrench$local.time <= as.POSIXct("2012-08-26 00:00:00", tz="Etc/GMT+7"), ]
vfrench1day <- vfrench[vfrench$local.time >= as.POSIXct("2012-08-24 04:00:00", tz="Etc/GMT+7") & 
                         vfrench$local.time <= as.POSIXct("2012-08-25 04:00:00", tz="Etc/GMT+7"), ]

test_that("prepdata_bayes and mcmc_bayes run with JAGS", {
  model_specs <- specs_bayes_jags_nopool_oi()
  model_specs$model_path <- system.file(paste0("models/bayes/", model_specs$model_file), package="streamMetabolizer") # usually added in metab_bayes
  
  # prepdata
  data_list <- streamMetabolizer:::prepdata_bayes(
    data=vfrench1day, data_daily=NULL, local_date="2012-08-24", model_specs=model_specs, priors=FALSE) 
  expect_equal(length(data_list$DO_obs), nrow(vfrench1day))
  expect_equal(length(data_list$frac_ER), nrow(vfrench1day))
  
  # mcmc
  mcmc_out <- do.call(streamMetabolizer:::mcmc_bayes, c(
    list(data_list=data_list),
    model_specs[c('bayes_software','model_path','params_out','n_chains','n_cores','adapt_steps','burnin_steps','num_saved_steps','thin_steps','verbose')]))
  expect_is(mcmc_out, 'data.frame')
  expect_true(all(c("GPP_daily_mean","GPP_daily_sd", "GPP_daily_2.5pct") %in% names(mcmc_out)))
  
  # in plys
  ply_out <- streamMetabolizer:::bayes_1ply(
    data_ply=vfrench1day, data_daily_ply=NULL, day_start=4, day_end=28, local_date="2012-08-24",
    tests=c('full_day', 'even_timesteps', 'complete_data'),
    model_specs=model_specs)
  expect_is(ply_out, 'data.frame')
  expect_true(all(c("GPP_daily_mean","GPP_daily_sd", "GPP_daily_2.5pct","warnings") %in% names(ply_out)))
  
})

dontrun <- function() {
  
  test_that("prepdata_bayes and mcmc_bayes run with Stan", {
    model_specs <- specs_bayes_stan_nopool_oi(n_cores=2)
    model_specs$model_path <- system.file(paste0("models/bayes/", model_specs$model_file), package="streamMetabolizer") # usually added in metab_bayes
    
    # prepdata
    data_list <- streamMetabolizer:::prepdata_bayes(
      data=vfrench1day, data_daily=NULL, local_date="2012-08-24", model_specs=model_specs, priors=FALSE) 
    expect_equal(length(data_list$DO_obs), nrow(vfrench1day))
    expect_equal(length(data_list$frac_ER), nrow(vfrench1day))
    
    # mcmc
    system.time({
      mcmc_out <- do.call(streamMetabolizer:::mcmc_bayes, c(
        list(data_list=data_list),
        model_specs[c('bayes_software','model_path','params_out','n_chains','n_cores','burnin_steps','num_saved_steps','thin_steps','verbose')]))
    })
    expect_is(mcmc_out, 'data.frame')
    expect_true(all(c("GPP_daily_mean","GPP_daily_sd", "GPP_daily_2.5pct") %in% names(mcmc_out)))
    
    # in plys
    ply_out <- streamMetabolizer:::bayes_1ply(
      data_ply=vfrench1day, data_daily_ply=NULL, day_start=4, day_end=28, local_date="2012-08-24",
      tests=c('full_day', 'even_timesteps', 'complete_data'),
      model_specs=model_specs)
    expect_is(ply_out, 'data.frame')
    expect_true(all(c("GPP_daily_mean","GPP_daily_sd", "GPP_daily_2.5pct","warnings") %in% names(ply_out)))
    
  })
  
  test_that("prepdata_bayes and mcmc_bayes run with Stan procobserr", {
    model_specs <- specs_bayes_stan_nopool_oipc(n_cores=4, burnin_steps=2000, num_saved_steps=1000, err_proc_phi_max=1)
    model_specs$model_path <- system.file(paste0("models/bayes/", model_specs$model_file), package="streamMetabolizer") # usually added in metab_bayes
    
    # prepdata
    data_list <- streamMetabolizer:::prepdata_bayes(
      data=vfrench1day, data_daily=NULL, local_date="2012-08-24", model_specs=model_specs, priors=FALSE) 
    expect_equal(length(data_list$DO_obs), nrow(vfrench1day))
    expect_equal(length(data_list$frac_ER), nrow(vfrench1day))
    
    # mcmc
    system.time({
      mcmc_out <- do.call(streamMetabolizer:::mcmc_bayes, c(
        list(data_list=data_list, keep_mcmc=TRUE),
        model_specs[c('bayes_software','model_path','params_out','n_chains','n_cores','burnin_steps','num_saved_steps','thin_steps','verbose')]))
    })
    expect_is(mcmc_out, 'data.frame')
    expect_true(all(c("GPP_daily_mean","GPP_daily_sd", "GPP_daily_2.5pct") %in% names(mcmc_out)))
    expect_output(show(mcmc_out$mcmcfit[[1]]), "Inference for Stan model")
    plot(mcmc_out$mcmcfit[[1]])
    pairs(mcmc_out$mcmcfit[[1]])
    traceplot(mcmc_out$mcmcfit[[1]])
    
    # in plys
    ply_out <- streamMetabolizer:::bayes_1ply(
      data_ply=vfrench1day, data_daily_ply=NULL, day_start=4, day_end=28, local_date="2012-08-24",
      tests=c('full_day', 'even_timesteps', 'complete_data'),
      model_specs=model_specs)
    expect_is(ply_out, 'data.frame')
    expect_true(all(c("GPP_daily_mean","GPP_daily_sd", "GPP_daily_2.5pct","warnings") %in% names(ply_out)))
    
  })
  
  test_that("yackulic model runs", {
    model_specs <- specs_bayes_stan_nopool_pcpi(n_cores=4, err_proc_acor_phi_min=0.99, err_proc_acor_sigma_max = 0.0001, burnin_steps = 3000, num_saved_steps = 1000)
    model_specs$model_path <- system.file(paste0("models/bayes/", model_specs$model_file), package="streamMetabolizer") # usually added in metab_bayes
    data_list <- streamMetabolizer:::prepdata_bayes(
      data=vfrench1day, data_daily=NULL, local_date="2012-08-24", model_specs=model_specs, priors=FALSE) 
    
    system.time({
      mcmc_out <- do.call(streamMetabolizer:::mcmc_bayes, c(
        list(data_list=data_list, keep_mcmc=TRUE),
        model_specs[c('bayes_software','model_path','params_out','n_chains','n_cores','burnin_steps','num_saved_steps','thin_steps','verbose')]))
    })
    show(mcmc_out$mcmcfit[[1]])
    plot(mcmc_out$mcmcfit[[1]])
    # pairs(mcmc_out$mcmcfit[[1]])
    traceplot(mcmc_out$mcmcfit[[1]])
    
  })
}

test_that("metab_bayes predictions (predict_metab, predict_DO) make sense", {
  
  # specs_bayes_jags_nopool_oi
  mm <- mmOP <- metab_bayes(data=vfrenchshort)
  metab <- predict_metab(mm)
  DO_preds <- predict_DO(mm)
  DO_preds_Aug24<- dplyr::filter(DO_preds, local.date == "2012-08-24")
  expect_true(all(abs(DO_preds_Aug24$DO.obs - DO_preds_Aug24$DO.mod) < 0.3), "DO.mod tracks DO.obs with not too much error")
  # now with Euler solution
  mm <- mmOE <- metab_bayes(data=vfrenchshort, model_specs=specs_bayes_jags_nopool_oi(
    model_file="nopool_oi_Euler.jags", num_saved_steps=500, GPP_daily_mu=2))
  metab <- predict_metab(mm)
  DO_preds <- predict_DO(mm)
  DO_preds_Aug24<- dplyr::filter(DO_preds, local.date == "2012-08-24")
  expect_true(all(abs(DO_preds_Aug24$DO.obs - DO_preds_Aug24$DO.mod) < 0.3), "DO.mod tracks DO.obs with not too much error")
  
  # specs_bayes_jags_nopool_oipc. you really have to crank down the err.proc.sigma.max or else the errors are huge
  mm <- mmOPP <- metab_bayes(data=vfrenchshort, model_specs=specs_bayes_jags_nopool_oipc(keep_mcmcs = TRUE))
  #plot(get_mcmc(mm)[[2]], "trace", vars=mm@args$model_specs$params_out, layout=c(2,3))
  metab <- predict_metab(mm)
  DO_preds <- predict_DO(mm)
  DO_preds_Aug24<- dplyr::filter(DO_preds, local.date == "2012-08-24")
  expect_true(all(abs(DO_preds_Aug24$DO.obs - DO_preds_Aug24$DO.mod) < 0.3), "DO.mod tracks DO.obs with not too much error")
  # now with Euler solution
  mm <- mmOPE <- metab_bayes(data=vfrenchshort, model_specs=specs_bayes_jags_nopool_oipc(
    model_file="nopool_oipc_Euler.jags", num_saved_steps=800, GPP_daily_mu=2, verbose=FALSE))
  metab <- predict_metab(mm)
  DO_preds <- predict_DO(mm)
  DO_preds_Aug24<- dplyr::filter(DO_preds, local.date == "2012-08-24")
  expect_true(all(abs(DO_preds_Aug24$DO.obs - DO_preds_Aug24$DO.mod) < 0.3), "DO.mod tracks DO.obs with not too much error")
  
  # these just take forever; keep them as manual tests for now
  dontrun <- function() {
    # specs_bayes_stan_nopool_oi
    mm <- mmSOP <- metab_bayes(data=vfrenchshort, model_specs=specs_bayes_stan_nopool_oi(n_cores=4))
    metab <- predict_metab(mm)
    expect_equal(metab$GPP.lower, get_fit(mm)$GPP_daily_2.5pct)
    DO_preds <- predict_DO(mm)
    DO_preds_Aug24<- dplyr::filter(DO_preds, local.date == "2012-08-24")
    expect_true(all(abs(DO_preds_Aug24$DO.obs - DO_preds_Aug24$DO.mod) < 0.3), "DO.mod tracks DO.obs with not too much error")
    
    # specs_bayes_stan_nopool_oipc
    mm <- mmSOPP <- metab_bayes(data=vfrenchshort, model_specs=specs_bayes_stan_nopool_oipc(n_cores=4, n_chains=4, keep_mcmcs="2012-08-24", burnin_steps = 3000, num_saved_steps = 1000))
    expect_is(get_fitting_time(mm), "proc_time")
    expect_is(get_mcmc(mm)[[2]], "stanfit")
    expect_equal(names(mm@mcmc)[2], "2012-08-24")
    predict_metab(mm)
    get_fit(mm)[grep("Rhat", names(get_fit(mm)))]
    expect_equal(metab$GPP.lower, get_fit(mm)$GPP_daily_2.5pct)
    DO_preds <- predict_DO(mm)
    DO_preds_Aug24<- dplyr::filter(DO_preds, local.date == "2012-08-24")
    expect_true(all(abs(DO_preds_Aug24$DO.obs - DO_preds_Aug24$DO.mod) < 0.3), "DO.mod tracks DO.obs with not too much error")
    plot_DO_preds(DO_preds)
    traceplot(get_mcmc(mm)$"2012-08-24")

    # plot_DO_preds(DO_preds)
  
    # compare models - since this includes mmSOP and mmSOPP, keep within the dontrun function
    library(dplyr); library(tidyr); library(ggplot2); library(gridExtra)
    preds <- bind_rows(lapply(list(mmOP, mmOE, mmOPP, mmOPE, mmSOP, mmSOPP), function(mm) 
      bind_cols(
        predict_metab(mm)[2,], 
        get_fit(mm)[2, grepl("Rhat|psrf", names(get_fit(mm)))] %>% setNames(gsub("Rhat|psrf", "rhat", names(.))), #potential scale reduction factor
        data_frame(file=get_args(mm)$model_specs$model_file))
    )) %>% mutate(model=LETTERS[1:nrow(.)])
    grid.arrange(
      ggplot(preds, aes(x=model, y=GPP, color=file)) + geom_point() + geom_errorbar(aes(ymin=GPP.lower, ymax=GPP.upper)) + theme_bw() + ylim(min(preds$GPP.lower), NA) + theme(legend.position="none"),
      ggplot(preds, aes(x=model, y=ER, color=file)) + geom_point() + geom_errorbar(aes(ymin=ER.lower, ymax=ER.upper)) + theme_bw() + ylim(NA, max(preds$ER.upper)) + theme(legend.position="none"),
      ggplot(preds, aes(x=model, y=K600, color=file)) + geom_point() + geom_errorbar(aes(ymin=K600.lower, ymax=K600.upper)) + theme_bw() + ylim(min(preds$K600.lower), NA),
      layout_matrix=matrix(c(1,2,3,3), byrow=TRUE, ncol=2))
    psrfs <- gather(select(preds, model, file, ends_with('rhat')), var, rhat, ends_with('rhat'))
    ggplot(psrfs, aes(y=var, x=rhat, color=file)) + geom_point(size=4) + theme_bw() + geom_vline(xintercept=1.1, color="blue")
  }
})
