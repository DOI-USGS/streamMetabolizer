context("metab_bayes")

test_that("bayes metab_spec functions work", {

  # ?specs_bayes should give useful help
  expect_is(specs_bayes_jags_nopool_obserr(), "list")
  expect_is(specs_bayes_jags_nopool_procobserr(), "list")
  
})

# get, format, & subset data
vfrench <- streamMetabolizer:::load_french_creek(attach.units=FALSE)
vfrenchshort <- vfrench[vfrench$local.time >= as.POSIXct("2012-08-23 00:00:00", tz="Etc/GMT+7") & 
                          vfrench$local.time <= as.POSIXct("2012-08-26 00:00:00", tz="Etc/GMT+7"), ]
vfrench1day <- vfrench[vfrench$local.time >= as.POSIXct("2012-08-24 04:00:00", tz="Etc/GMT+7") & 
                         vfrench$local.time <= as.POSIXct("2012-08-25 04:00:00", tz="Etc/GMT+7"), ]

test_that("prepdata_bayes and mcmc_bayes run with JAGS", {
  model_specs <- specs_bayes_jags_nopool_obserr()
  model_specs$model_path <- system.file(paste0("models/bayes/", model_specs$model_file), package="streamMetabolizer") # usually added in metab_bayes
  
  # prepdata
  data_list <- streamMetabolizer:::prepdata_bayes(
    data=vfrench1day, data_daily=NULL, local_date="2012-08-24", model_specs=model_specs, priors=FALSE) 
  expect_equal(length(data_list$DO_obs), nrow(vfrench1day))
  expect_equal(length(data_list$frac_ER), nrow(vfrench1day))
  
  # mcmc
  mcmc_out <- do.call(streamMetabolizer:::mcmc_bayes, c(
    list(data_list=data_list),
    model_specs[c('bayes_software','model_path','params_out','max_cores','adapt_steps','burnin_steps','num_saved_steps','thin_steps','verbose')]))
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

test_that("prepdata_bayes and mcmc_bayes run with Stan", {
  model_specs <- specs_bayes_stan_nopool_obserr()
  model_specs$model_path <- system.file(paste0("models/bayes/", model_specs$model_file), package="streamMetabolizer") # usually added in metab_bayes
  
  # prepdata
  data_list <- streamMetabolizer:::prepdata_bayes(
    data=vfrench1day, data_daily=NULL, local_date="2012-08-24", model_specs=model_specs, priors=FALSE) 
  expect_equal(length(data_list$DO_obs), nrow(vfrench1day))
  expect_equal(length(data_list$frac_ER), nrow(vfrench1day))
  
  # mcmc
  mcmc_out <- do.call(streamMetabolizer:::mcmc_bayes, c(
    list(data_list=data_list),
    model_specs[c('bayes_software','model_path','params_out','max_cores','adapt_steps','burnin_steps','num_saved_steps','thin_steps','verbose')]))
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


test_that("metab_bayes predictions (predict_metab, predict_DO) make sense", {
  
  # specs_bayes_jags_nopool_obserr
  mmOP <- mm <- metab_bayes(data=vfrenchshort)
  metab <- predict_metab(mm)
  DO_preds <- predict_DO(mm)
  DO_preds_Aug24<- dplyr::filter(DO_preds, local.date == "2012-08-24")
  expect_true(all(abs(DO_preds_Aug24$DO.obs - DO_preds_Aug24$DO.mod) < 0.3), "DO.mod tracks DO.obs with not too much error")
  # now with Euler solution
  mmOE <- mm <- metab_bayes(data=vfrenchshort, model_specs=specs_bayes_jags_nopool_obserr(
    model_file="nopool_obserr_Euler.jags", num_saved_steps=500, GPP_daily_mu=2))
  metab <- predict_metab(mm)
  DO_preds <- predict_DO(mm)
  DO_preds_Aug24<- dplyr::filter(DO_preds, local.date == "2012-08-24")
  expect_true(all(abs(DO_preds_Aug24$DO.obs - DO_preds_Aug24$DO.mod) < 0.3), "DO.mod tracks DO.obs with not too much error")
  
  # specs_bayes_jags_nopool_procobserr. you really have to crank down the err.proc.sigma.max or else the errors are huge
  mmOPP <- mm <- metab_bayes(data=vfrenchshort, model_specs=specs_bayes_jags_nopool_procobserr(err_proc_sigma_max = 0.0005))
  metab <- predict_metab(mm)
  DO_preds <- predict_DO(mm)
  DO_preds_Aug24<- dplyr::filter(DO_preds, local.date == "2012-08-24")
  expect_true(all(abs(DO_preds_Aug24$DO.obs - DO_preds_Aug24$DO.mod) < 0.3), "DO.mod tracks DO.obs with not too much error")
  # now with Euler solution
  mmOPE <- mm <- metab_bayes(data=vfrenchshort, model_specs=specs_bayes_jags_nopool_procobserr(
    model_file="nopool_procobserr_Euler.jags", num_saved_steps=800, GPP_daily_mu=2, verbose=FALSE))
  metab <- predict_metab(mm)
  DO_preds <- predict_DO(mm)
  DO_preds_Aug24<- dplyr::filter(DO_preds, local.date == "2012-08-24")
  expect_true(all(abs(DO_preds_Aug24$DO.obs - DO_preds_Aug24$DO.mod) < 0.3), "DO.mod tracks DO.obs with not too much error")
  
  # specs_bayes_stan_nopool_obserr
  mmSOP <- mm <- metab_bayes(data=vfrenchshort, model_specs=specs_bayes_stan_nopool_obserr())
  metab <- predict_metab(mm)
  expect_equal(metab$GPP.lower, get_fit(mm)$GPP_daily_2.5pct)
  DO_preds <- predict_DO(mm)
  DO_preds_Aug24<- dplyr::filter(DO_preds, local.date == "2012-08-24")
  expect_true(all(abs(DO_preds_Aug24$DO.obs - DO_preds_Aug24$DO.mod) < 0.3), "DO.mod tracks DO.obs with not too much error")
  
  # plot_DO_preds(DO_preds)
  
  # compare models
  # library(dplyr); library(ggplot2)
  # preds <- bind_rows(lapply(list(mmOP, mmOE, mmOPP, mmOPE), function(mm) 
  #   bind_cols(predict_metab(mm)[2,], data_frame(file=get_args(mm)$model_specs$model_file))
  # ))
  # ggplot(preds, aes(x=file, y=GPP)) + geom_point() + geom_errorbar(aes(ymin=GPP.lower, ymax=GPP.upper)) + theme_bw() + ylim(0, NA)
  # ggplot(preds, aes(x=file, y=ER)) + geom_point() + geom_errorbar(aes(ymin=ER.lower, ymax=ER.upper)) + theme_bw() + ylim(NA, 0)
  # ggplot(preds, aes(x=file, y=K600)) + geom_point() + geom_errorbar(aes(ymin=K600.lower, ymax=K600.upper)) + theme_bw() + ylim(0, NA)
  
})
