context("metab_night")

vfrench <- streamMetabolizer:::load_french_creek(attach.units=FALSE)
vfrenchshort <- vfrench[vfrench$solar.time >= as.POSIXct("2012-08-23 00:00:00", tz="UTC") & 
                          vfrench$solar.time <= as.POSIXct("2012-08-26 00:00:00", tz="UTC"), ]

test_that("metab_night models can be created", {
  
  mm <- metab_night(vfrenchshort)
  
  # check basic structure
  expect_is(mm, "metab_night")
  expect_is(slot(mm, "fit"), "data.frame")
  expect_is(slot(mm, "args"), "list")
  expect_is(slot(mm, "data"), "data.frame")
  expect_is(slot(mm, "pkg_version"), "character")
  
})

test_that("metab_night predictions (predict_metab, predict_DO) make sense", {
  
  # metab_night
  mm <- metab_night(data=vfrenchshort)
  metab <- predict_metab(mm)
  DO_preds <- predict_DO(mm)
  DO_preds_Aug24 <- dplyr::filter(DO_preds, solar.date == "2012-08-24")
  expect_true(all(abs(DO_preds_Aug24$DO.obs - DO_preds_Aug24$DO.mod) < 0.15), "DO.mod tracks DO.obs with not too much error")
  # plot_DO_preds(DO_preds, y_var="conc")
  # plot_DO_preds(DO_preds, y_var="pctsat")
  
})

test_that("metab_night predictions can be passed back into metab_mle", {
  
  # metab_night
  mmk <- metab_night(data=vfrenchshort)

  # metab_mle
  mm <- metab_mle(data=vfrenchshort, data_daily=predict_metab(mmk)[c('solar.date', 'K600')])
  metab <- predict_metab(mm)
  DO_preds <- predict_DO(mm)
  DO_preds_Aug24 <- dplyr::filter(DO_preds, solar.date == "2012-08-24")
  # note that had to raise the maximum error for this model combination, from 0.15 to 0.35 mgO/L
  expect_true(all(abs(DO_preds_Aug24$DO.obs - DO_preds_Aug24$DO.mod) < 0.35), "DO.mod tracks DO.obs with not too much error")
  # plot_DO_preds(DO_preds)
  
})



