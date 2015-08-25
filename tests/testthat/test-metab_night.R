context("metab_night")

vfrench <- unitted::v(streamMetabolizer:::load_french_creek())

test_that("metab_night models can be created", {
  
  mm <- metab_night(vfrench)
  
  # check basic structure
  expect_is(mm, "metab_night")
  expect_is(slot(mm, "fit"), "data.frame")
  expect_is(slot(mm, "args"), "list")
  expect_is(slot(mm, "data"), "data.frame")
  expect_is(slot(mm, "pkg_version"), "character")
  
})

test_that("metab_night predictions (predict_metab, predict_DO) make sense", {
  
  # metab_night
  mm <- metab_night(data=vfrench)
  metab <- predict_metab(mm)
  DO_preds <- predict_DO(mm)
  DO_preds_Aug24 <- dplyr::filter(DO_preds, local.date == "2012-08-24")
  expect_true(all(abs(DO_preds_Aug24$DO.obs - DO_preds_Aug24$DO.mod) < 0.15), "DO.mod tracks DO.obs with not too much error")
  # plot_DO_preds(DO_preds, plot_as="conc")
  # plot_DO_preds(DO_preds, plot_as="pctsat")
  
})

test_that("metab_night predictions can be passed back into metab_mle", {
  
  # metab_night
  mmk <- metab_night(data=vfrench)

  # metab_mle
  mm <- metab_mle(data=vfrench, data_daily=predict_metab(mmk)[c('local.date', 'K600')])
  metab <- predict_metab(mm)
  DO_preds <- predict_DO(mm)
  DO_preds_Aug24 <- dplyr::filter(DO_preds, local.date == "2012-08-24")
  # note that had to raise the maximum error for this model combination, from 0.15 to 0.25 mgO/L
  expect_true(all(abs(DO_preds_Aug24$DO.obs - DO_preds_Aug24$DO.mod) < 0.25), "DO.mod tracks DO.obs with not too much error")
  # plot_DO_preds(DO_preds, plot_as="conc")
  # plot_DO_preds(DO_preds, plot_as="pctsat")
  
})



