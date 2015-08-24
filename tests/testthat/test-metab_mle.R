context("metab_mle")

vfrench <- unitted::v(streamMetabolizer:::load_french_creek())

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
  mm <- metab_mle(data=vfrench, day_start=-1, day_end=23)
  metab <- predict_metab(mm)
  DO_preds <- predict_DO(mm)
  DO_preds_Aug24 <- dplyr::filter(DO_preds, local.date == "2012-08-24")
  expect_true(all(abs(DO_preds_Aug24$DO.obs - DO_preds_Aug24$DO.mod) < 0.15), "DO.mod tracks DO.obs with not too much error")
  # plot_DO_preds(DO_preds, plot_as="pctsat")
  
})

test_that("metab_mle models can be fit with K specified", {
  
  # metab_mle with K600
  K600 <- data.frame(local.date=unique(as.Date(vfrench$local.time)), K600=c(NA, 30, NA, 0, NA, 0, 50, 40)) # 2 matches, one mismatch w/ modeled data
  mm <- metab_mle(data=vfrench, data_daily=K600, day_start=-1, day_end=23)
  metab <- predict_metab(mm)
  DO_preds <- predict_DO(mm)
  DO_preds_Aug24 <- dplyr::filter(DO_preds, local.date == "2012-08-24")
  expect_true(all(abs(DO_preds_Aug24$DO.obs - DO_preds_Aug24$DO.mod) < 0.15), "DO.mod tracks DO.obs with not too much error")
  # plot_DO_preds(DO_preds, plot_as="pctsat")
  
})
