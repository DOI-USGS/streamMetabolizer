context("metab_sim")

test_that("metab_sim predictions (predict_metab, predict_DO) make sense", {
  
  # generate data
  dat <- data_metab('3')
  data_date <- mm_model_by_ply(mm_model_by_ply_prototype, dat, day_start=4, day_end=28)$date
  dd <- data.frame(date=data_date, DO.mod.1=7.5, GPP=4, ER=-c(NA,2,4), K600=30)
  
  # should be able to fit by specifying either data$DO.obs[d,1] or data_daily$DO.mod.1 
  mm <- metab_sim(data=select(dat, -DO.obs), data_daily=dd)
  mm2 <- metab_sim(data=dat, data_daily=select(dd, -DO.mod.1))
  expect_equal(select(get_fit(mm), -DO.mod.1), get_fit(mm2))
  
  # predict_metab
  expect_equal(get_fit(mm), select(predict_metab(mm), -warnings, -errors))
  expect_equal(get_fit(mm2), select(predict_metab(mm2), -warnings, -errors))
  
  # predict_DO - DO.mod (no error) and DO.obs (with any error) should still be pretty close
  expect_true(rmse_DO(predict_DO(mm)) < 0.15, "DO.mod tracks DO.obs with not too much error")
  # plot_DO_preds(predict_DO(mm))
  
  # predict_DO - DO.obs should be different each time unless seed is set. DO.mod should always be the same
  expect_true(!isTRUE(all.equal(predict_DO(mm)$DO.obs, predict_DO(mm)$DO.obs)))
  expect_true(isTRUE(all.equal(predict_DO(mm)$DO.mod, predict_DO(mm)$DO.mod)))
  mm <- metab_sim(data=dat, data_daily=select(dd, -DO.mod.1), 
                  specs=specs('s_np_oipcpi_eu_.rnorm', sim.seed=626, day_start=-1, day_end=23))
  expect_true(isTRUE(all.equal(predict_DO(mm)$DO.obs, predict_DO(mm)$DO.obs)))
  expect_true(isTRUE(all.equal(predict_DO(mm)$DO.mod, predict_DO(mm)$DO.mod)))
  
  # predict_DO - using default (just err.obs.sigma), should have basically no autocorrelation in errors
  mm <- metab_sim(data=select(dat, -DO.obs), data_daily=dd)
  DO_preds <- predict_DO(mm, date_start="2012-09-19")
  acf_out <- acf(DO_preds$DO.mod - DO_preds$DO.obs, plot=FALSE)
  expect_lt(acf_out$acf[acf_out$lag==1], 0.1)
  # plot_DO_preds(predict_DO(mm))
  
  # predict_DO - autocorrelation should be bigger when there's process error
  mm <- metab_sim(data=select(dat, -DO.obs), data_daily=dd,
                  specs=specs('s_np_oipcpi_eu_.rnorm', err.obs.sigma=0, err.proc.sigma=0.05))
  DO_preds <- predict_DO(mm, date_start="2012-09-19")
  acf_out <- acf(DO_preds$DO.mod - DO_preds$DO.obs, plot=FALSE)
  expect_gt(acf_out$acf[acf_out$lag==1], 0.6)
  # plot_DO_preds(predict_DO(mm))
  
  # should be able to switch ODE methods in fitting
  dat <- data_metab('3', res='30')
  mmE <- metab_sim(specs('s_np_oipcpi_eu_.rnorm', err.obs.sigma=0, err.proc.sigma=0.05, sim.seed=4),
                   data=select(dat, -DO.obs), data_daily=dd)
  mmP <- metab_sim(specs('s_np_oipcpi_pm_.rnorm', err.obs.sigma=0, err.proc.sigma=0.05, sim.seed=4),
                  data=select(dat, -DO.obs), data_daily=dd)
  # DO_preds <- bind_rows(
  #   data.frame(predict_DO(mmE), method="Euler", stringsAsFactors=FALSE),
  #   data.frame(predict_DO(mmP), method="pairmeans", stringsAsFactors=FALSE))
  # library(ggplot2)
  # ggplot(DO_preds, aes(x=solar.time, y=100*DO.mod/DO.sat, color=method)) + geom_line() + theme_bw()
  # ggplot(DO_preds, aes(x=solar.time, y=100*DO.obs/DO.sat, color=method)) + geom_line() + theme_bw()
  
})
