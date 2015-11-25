context("metab_mle")

vfrench <- streamMetabolizer:::load_french_creek(attach.units=FALSE)
vfrenchshort <- vfrench[vfrench$local.time >= as.POSIXct("2012-08-23 00:00:00", tz="Etc/GMT+7") & 
                          vfrench$local.time <= as.POSIXct("2012-08-30 00:00:00", tz="Etc/GMT+7"), ]

test_that("metab_mle models can be created", {
  
  mm <- metab_mle(data=vfrench)
  
  # check basic structure
  expect_is(mm, "metab_mle")
  expect_is(slot(mm, "fit"), "data.frame")
  expect_is(slot(mm, "args"), "list")
  expect_is(slot(mm, "data"), "data.frame")
  expect_is(slot(mm, "pkg_version"), "character")
  
})

test_that("metab_mle predictions (predict_metab, predict_DO) make sense", {
  
  # metab_mle
  mm <- metab_mle(data=vfrenchshort, day_start=-1, day_end=23)
  metab <- predict_metab(mm)
  DO_preds <- predict_DO(mm)
  DO_preds_Aug24 <- dplyr::filter(DO_preds, local.date == "2012-08-24")
  expect_true(all(abs(DO_preds_Aug24$DO.obs - DO_preds_Aug24$DO.mod) < 0.25), "DO.mod tracks DO.obs with not too much error")
  # plot_DO_preds(DO_preds_Aug24)
  # plot_DO_preds(DO_preds)
  
  # fit with different ODE methods
  mmE <- metab_mle(data=vfrenchshort, day_start=-1, day_end=23, model_specs=specs('m_np_oi_eu_km.nlm'))
  mmP <- metab_mle(data=vfrenchshort, day_start=-1, day_end=23, model_specs=specs('m_np_oi_pm_km.nlm'))
  plot_DO_preds(predict_DO(mmE))
  plot_DO_preds(predict_DO(mmP))
  predict_metab(mmE) - predict_metab(mmP)

  # predict with different ODE methods
  mm <- metab_mle(data=vfrenchshort, day_start=-1, day_end=23)
  plot_DO_preds(predict_DO(mm, calc_DO_args=list(ODE_method="Euler")))
  plot_DO_preds(predict_DO(mm))
})

test_that("metab_mle models can be fit with K specified", {
  
  # metab_mle with K600
  K600 <- data.frame(local.date=unique(as.Date(vfrenchshort$local.time)), K600=c(NA, 30, NA, 50, 40)) # 2 matches, one mismatch w/ modeled data
  mm <- metab_mle(data=vfrenchshort, data_daily=K600, day_start=-1, day_end=23)
  metab <- predict_metab(mm)
  DO_preds <- predict_DO(mm)
  DO_preds_Aug24 <- dplyr::filter(DO_preds, local.date == "2012-08-24")
  expect_true(all(abs(DO_preds_Aug24$DO.obs - DO_preds_Aug24$DO.mod) < 0.30), "DO.mod tracks DO.obs with not too much error")
  # plot_DO_preds(DO_preds, y_var="pctsat")
  
})
