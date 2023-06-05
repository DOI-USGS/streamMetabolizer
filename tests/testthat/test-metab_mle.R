context("metab_mle")

test_that("metab_mle models can be created", {

  mm <- metab_mle(data=data_metab('1', res='30'))

  # check basic structure
  expect_is(mm, "metab_mle")
  expect_is(slot(mm, "fit"), "data.frame")
  expect_is(slot(mm, "specs"), "list")
  expect_is(slot(mm, "data"), "data.frame")
  expect_is(slot(mm, "pkg_version"), "character")

  # should work when data is a tibble, too
  mmt <- metab_mle(data=tibble::as_tibble(data_metab('1', res='30')))
  expect_equal(get_params(mmt), get_params(mmt))
})

test_that("metab_mle works with fancy GPP, ER functions", {
  # setup for example fanciness
  satlight_q10temp_params <- c('Pmax','alpha','ER20','K600.daily')
  dat <- data_metab('1', '30')

  # specs should contain inits for the relevant parameters
  sp <- specs(mm_name('mle', GPP_fun='satlight', ER_fun='q10temp'))
  expect_true(all(paste0('init.',satlight_q10temp_params) %in% names(sp)))

  # model fitting should run without error
  mm <- metab_mle(sp, dat)

  # get_params should return values for the relevant parameters
  expect_true(all(satlight_q10temp_params %in% names(get_params(mm))))

  # predict_metab should return values
  mp <- predict_metab(mm)
  expect_true(!is.na(mp$GPP))
  expect_true(!is.na(mp$ER))
})

test_that("metab_mle treats data flaws correctly", {

  # missing end
  dat <- data_metab('3','15',flaws='missing end')
  mm <- metab_mle(data=dat)
  expect_true(is.na(get_params(mm)$GPP.daily[3]))
  expect_equal(get_params(mm)$errors[3], "data don't start when expected")

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
  mmE <- metab_mle(specs('m_np_oi_eu_plrckm.nlm'), data=dat3)
  mmP <- metab_mle(specs('m_np_oi_tr_plrckm.nlm'), data=dat3)
  expect_lt(rmse_DO(predict_DO(mmP)), rmse_DO(predict_DO(mmE))) #, info="trapezoid should be more accurate than euler")
  # plot_DO_preds(predict_DO(mmE))
  # plot_DO_preds(predict_DO(mmP))

  # fix K600 (uses mmP from above)
  K600 <- dplyr::transmute(predict_metab(mmP), date=date, K600.daily=c(NA, 20, 21))
  mmK <- metab_mle(get_specs(mmP), data=dat3, data_daily=K600)
  expect_equal(predict_metab(mmK)[1,1:6], predict_metab(mmP)[1,1:6]) # whole first date should be identical
  expect_equal(get_params(mmK)[2:3,"K600.daily"], K600[2:3,"K600.daily"]) # K got fixed on days 2 & 3
  # plot_DO_preds(predict_DO(mmK), y_var="pctsat")
})

test_that("metab_mle outputs look like Bob's", {

  dat <- data_metab('1', day_start=-2, day_end=30)

  # PRK
  mms <- metab_mle(specs(mm_name('mle', ode_method='euler'), day_start=-2, day_end=30), data=dat)
  mmb <- streamMetabolizer:::load_french_creek_std_mle(
    dat, estimate='PRK',
    start=c(dates="09/17/12", times="22:00:00"),
    end=c(dates="09/19/12", times="06:00:00"))
  expect_equal(predict_metab(mms)[1,"GPP"], mmb[1,"GPP"], tol=0.001) # we handle light slightly differently. i prefer the sM way
  expect_equal(predict_metab(mms)[1,"ER"], mmb[1,"ER"], tol=0.0001)
  expect_equal(get_params(mms)[1,"K600.daily"], mmb[1,"K"], tol=0.0001)
  expect_equal(get_fit(mms)[1,"minimum"], mmb[1,"lik"], tol=0.00001)

  # PR
  mms <- metab_mle(specs(mm_name('mle', ode_method='euler'), day_start=-2, day_end=30),
                   data=dat, data_daily=data.frame(date=as.Date("2012-09-18"), K600.daily=35))
  mmb <- streamMetabolizer:::load_french_creek_std_mle(
    dat, estimate='PR', K=35,
    start=c(dates="09/17/12", times="22:00:00"),
    end=c(dates="09/19/12", times="06:00:00"))
  expect_equal(predict_metab(mms)[1,"GPP"], mmb[1,"GPP"], tol=0.001) # we handle light slightly differently. i prefer the sM way
  expect_equal(predict_metab(mms)[1,"ER"], mmb[1,"ER"], tol=0.00001)
  expect_equal(get_params(mms)[1,"K600.daily"], mmb[1,"K"], tol=0.00001)
  expect_equal(get_fit(mms)[1,"minimum"], mmb[1,"lik"], tol=0.00001)

})

test_that("metab_models can be saved & reloaded efficiently (see helper-save_load.R)", {

  skip("skipping check of gz9 as a great file compression choice; just doesn't seem that important")
    
  # save and reload
  mm <- metab_mle(data=data_metab('1'))
  
  # see if saveRDS with gzfile, compression=9 works well
  rdstimes <- save_load_timing(mm, reps=1) # autoloaded b/c script begins with 'helper' and is in this directory
  expect_true('gz6' %in% rdstimes$typelevel[1:3], info="gz6 is reasonably efficient for saveRDS")
  # plot_save_load_timing(rdstimes)
  
})

test_that("metab_models can be saved & reloaded without loss of information (see helper-save_load.R)", {

  # save and reload
  mm <- metab_mle(data=data_metab('1'))

  # save and load the mm, make sure it stays the same
  test_save_load_recovery(mm)

})
