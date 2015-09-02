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
  mm <- metab_bayes(data=vfrenchshort, model_specs=specs_bayes_jags_nopool_obserr())
  metab <- predict_metab(mm)
  DO_preds <- predict_DO(mm)
  DO_preds_Aug24<- dplyr::filter(DO_preds, local.date == "2012-08-24")
  expect_true(all(abs(DO_preds_Aug24$DO.obs - DO_preds_Aug24$DO.mod) < 0.3), "DO.mod tracks DO.obs with not too much error")
  # plot_DO_preds(DO_preds, plot_as="pctsat")
  
  # specs_bayes_jags_nopool_procobserr
  mm <- metab_bayes(data=vfrenchshort, model_specs=specs_bayes_jags_nopool_procobserr())
  metab <- predict_metab(mm)
  DO_preds <- predict_DO(mm)
  DO_preds_Aug24<- dplyr::filter(DO_preds, local.date == "2012-08-24")
  expect_true(all(abs(DO_preds_Aug24$DO.obs - DO_preds_Aug24$DO.mod) < 0.3), "DO.mod tracks DO.obs with not too much error")
  #plot_DO_preds(DO_preds, plot_as="pctsat")
  
})
