context("metab_mle")

test_that("metab_mle models can be created", {
  
  mm <- metab_mle(data=data_metab('1', res='30'))
  
  # check basic structure
  expect_is(mm, "metab_mle")
  expect_is(slot(mm, "fit"), "data.frame")
  expect_is(slot(mm, "specs"), "list")
  expect_is(slot(mm, "data"), "data.frame")
  expect_is(slot(mm, "pkg_version"), "character")
  
})

test_that("metab_mle predictions (predict_metab, predict_DO) make sense", {
  
  # 1 day
  mm <- metab_mle(specs(mm_name('mle'), day_start=-1, day_end=23), data=data_metab('1', day_start=-1, day_end=23))
  metab <- predict_metab(mm)
  DO_preds <- predict_DO(mm)
  expect_true(rmse_DO(DO_preds) < 0.2, info="DO.mod tracks DO.obs with not too much error")
  # plot_DO_preds(DO_preds)
  
  # 10 days
  mm <- metab_mle(specs(mm_name('mle'), day_start=2, day_end=26), data=data_metab('10', day_start=2, day_end=26))
  metab <- predict_metab(mm)
  DO_preds <- predict_DO(mm)
  expect_true(rmse_DO(DO_preds) < 0.2, info="DO.mod tracks DO.obs with not too much error")
  # plot_DO_preds(DO_preds)

  # compare ODE methods (3 days, default day_start&end)
  dat3 <- data_metab('3', res='30')
  mmE <- metab_mle(specs('m_np_oi_eu_km.nlm'), data=dat3)
  mmP <- metab_mle(specs('m_np_oi_pm_km.nlm'), data=dat3)
  expect_lt(rmse_DO(predict_DO(mmP)), rmse_DO(predict_DO(mmE))) #, info="pairmeans should be more accurate than Euler")
  # plot_DO_preds(predict_DO(mmE))
  # plot_DO_preds(predict_DO(mmP))
  
  # fix K600 (uses mmP from above)
  K600 <- dplyr::transmute(predict_metab(mmP), date=date, K600=c(NA, 20, 21))
  mmK <- metab_mle(get_specs(mmP), data=dat3, data_daily=K600)
  expect_equal(predict_metab(mmK)[1,1:10], predict_metab(mmP)[1,1:10]) # whole first date should be identical
  expect_equal(predict_metab(mmK)[2:3,"K600"], K600$K600[2:3]) # K got fixed on days 2 & 3
  # plot_DO_preds(predict_DO(mmK), y_var="pctsat")
})

test_that("metab_mle outputs look like Bob's", {
  
  dat <- data_metab('1', day_start=-2, day_end=30)
  
  # PRK
  mms <- metab_mle(specs(mm_name('mle', ode_method='Euler'), day_start=-2, day_end=30), data=dat)
  mmb <- streamMetabolizer:::load_french_creek_std_mle(
    dat, estimate='PRK', 
    start=c(dates="09/17/12", times="22:00:00"),
    end=c(dates="09/19/12", times="06:00:00"))
  expect_equal(get_fit(mms)[1,"GPP"], mmb[1,"GPP"], tol=0.001) # we handle light slightly differently. i prefer the sM way
  expect_equal(get_fit(mms)[1,"ER"], mmb[1,"ER"], tol=0.00001)
  expect_equal(get_fit(mms)[1,"K600"], mmb[1,"K"], tol=0.00001)
  expect_equal(get_fit(mms)[1,"minimum"], mmb[1,"lik"], tol=0.00001)
  
  # PR
  mms <- metab_mle(specs(mm_name('mle', ode_method='Euler'), day_start=-2, day_end=30), 
                   data=dat, data_daily=data.frame(date=as.Date("2012-09-18"), K600=35))
  mmb <- streamMetabolizer:::load_french_creek_std_mle(
    dat, estimate='PR', K=35, 
    start=c(dates="09/17/12", times="22:00:00"),
    end=c(dates="09/19/12", times="06:00:00"))
  expect_equal(get_fit(mms)[1,"GPP"], mmb[1,"GPP"], tol=0.001) # we handle light slightly differently. i prefer the sM way
  expect_equal(get_fit(mms)[1,"ER"], mmb[1,"ER"], tol=0.00001)
  expect_equal(get_fit(mms)[1,"K600"], mmb[1,"K"], tol=0.00001)
  expect_equal(get_fit(mms)[1,"minimum"], mmb[1,"lik"], tol=0.00001)
  
})

test_that("metab_models can be saved & reloaded (see helper-save_load.R)", {
  
  # save and reload
  mm <- metab_mle(data=data_metab('1'))
  
  # see if saveRDS with gzfile, compression=9 works well
  rdstimes <- save_load_timing(mm, reps=1) # autoloaded b/c script begins with 'helper' and is in this directory
  expect_true('gz6' %in% rdstimes$typelevel[1:3], info="gz6 is reasonably efficient for saveRDS")
  # plot_save_load_timing(rdstimes)
  
  # save and load the mm, make sure it stays the same
  test_save_load_recovery(mm)
  
})
