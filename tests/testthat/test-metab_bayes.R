context("metab_bayes")

test_that("bayes metab_spec functions work", {

  # ?specs_bayes should give useful help
  expect_is(specs_bayes_jags_nopool_obserr(), "list")
  expect_is(specs_bayes_jags_nopool_procobserr(), "list")
  
})

vfrench <- unitted::v(streamMetabolizer:::load_french_creek())

test_that("metab_bayes predictions (predict_metab, predict_DO) make sense", {
  
  # specs_bayes_jags_nopool_obserr
  mm <- metab_bayes(data=vfrench, model_specs=specs_bayes_jags_nopool_obserr())
  metab <- predict_metab(mm)
  DO_preds <- predict_DO(mm)
  DO_preds_Aug24<- dplyr::filter(DO_preds, local.date == "2012-08-24")
  expect_true(all(abs(DO_preds_Aug24$DO.obs - DO_preds_Aug24$DO.mod) < 0.2), "DO.mod tracks DO.obs with not too much error")
  #plot_DO_preds(DO_preds, plot_as="pctsat")
  
  # specs_bayes_jags_nopool_procobserr
  mm <- metab_bayes(data=vfrench, model_specs=specs_bayes_jags_nopool_procobserr())
  metab <- predict_metab(mm)
  DO_preds <- predict_DO(mm)
  DO_preds_Aug24<- dplyr::filter(DO_preds, local.date == "2012-08-24")
  expect_true(all(abs(DO_preds_Aug24$DO.obs - DO_preds_Aug24$DO.mod) < 0.2), "DO.mod tracks DO.obs with not too much error")
  #plot_DO_preds(DO_preds, plot_as="pctsat")
  
})
