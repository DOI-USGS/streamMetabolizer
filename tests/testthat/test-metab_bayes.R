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

test_that("metab_bayes predictions (predict_metab, predict_DO) make sense", {
  
  # specs_bayes_jags_nopool_obserr
  mmOP <- mm <- metab_bayes(data=vfrenchshort)
  metab <- predict_metab(mm)
  DO_preds <- predict_DO(mm)
  DO_preds_Aug24<- dplyr::filter(DO_preds, local.date == "2012-08-24")
  expect_true(all(abs(DO_preds_Aug24$DO.obs - DO_preds_Aug24$DO.mod) < 0.3), "DO.mod tracks DO.obs with not too much error")
  # now with Euler solution
  mmOE <- mm <- metab_bayes(data=vfrenchshort, model_specs=specs_bayes_jags_nopool_obserr(
    model_file="jags/nopool_obserr_Euler.txt", num_saved_steps=500, GPP.daily.mu=2))
  metab <- predict_metab(mm)
  DO_preds <- predict_DO(mm)
  DO_preds_Aug24<- dplyr::filter(DO_preds, local.date == "2012-08-24")
  expect_true(all(abs(DO_preds_Aug24$DO.obs - DO_preds_Aug24$DO.mod) < 0.3), "DO.mod tracks DO.obs with not too much error")
  
  # specs_bayes_jags_nopool_procobserr. you really have to crank down the err.proc.sigma.max or else the errors are huge
  mmOPP <- mm <- metab_bayes(data=vfrenchshort, model_specs=specs_bayes_jags_nopool_procobserr(err.proc.sigma.max = 0.0005))
  metab <- predict_metab(mm)
  DO_preds <- predict_DO(mm)
  DO_preds_Aug24<- dplyr::filter(DO_preds, local.date == "2012-08-24")
  expect_true(all(abs(DO_preds_Aug24$DO.obs - DO_preds_Aug24$DO.mod) < 0.3), "DO.mod tracks DO.obs with not too much error")
  # now with Euler solution
  mmOPE <- mm <- metab_bayes(data=vfrenchshort, model_specs=specs_bayes_jags_nopool_procobserr(
    model_file="jags/nopool_procobserr_Euler.txt", num_saved_steps=4000, GPP.daily.mu=2))
  metab <- predict_metab(mm)
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
