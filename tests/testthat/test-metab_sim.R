context("metab_sim")

test_that("metab_sim predictions (predict_metab, predict_DO) make sense", {
  
  library(dplyr)
  
  # generate data
  dat <- data_metab('3')
  data_date <- mm_model_by_ply(mm_model_by_ply_prototype, dat, day_start=4, day_end=28)$date
  dd <- data.frame(date=data_date, DO.mod.1=7.5, GPP.daily=4, ER.daily=-c(NA,2,4), K600.daily=30)
  
  # should be able to fit by specifying either data$DO.obs[d,1] or data_daily$DO.mod.1 
  mm <- metab_sim(data=select(dat, -DO.obs), data_daily=dd)
  mm2 <- metab_sim(data=dat, data_daily=select(dd, -DO.mod.1))
  
  # get_params
  expect_equal(select(get_params(mm), -DO.mod.1), get_params(mm2))
  
  # predict_metab
  expect_equal(select(get_params(mm), GPP=GPP.daily, ER=ER.daily), select(predict_metab(mm), GPP, ER))
  expect_equal(select(get_params(mm2), GPP=GPP.daily, ER=ER.daily), select(predict_metab(mm2), GPP, ER))
  
  # predict_DO - DO.mod.1 should follow specifications
  expect_equal(predict_DO(mm) %>% group_by(date) %>% summarize(first.DO.mod = DO.mod[1]) %>% .$first.DO.mod,
               dd %>% .$DO.mod.1 %>% {.*c(NA,1,1)} )
  expect_equal(predict_DO(mm2) %>% group_by(date) %>% summarize(first.DO.mod = DO.mod[1]) %>% .$first.DO.mod, 
               dat %>% filter(format(solar.time, '%H:%M') == '04:00') %>% .$DO.obs %>% {.*c(NA,1,1)} )
               
  # predict_DO - DO.mod (no error) and DO.obs (with any error) should still be pretty close
  expect_true(rmse_DO(predict_DO(mm)) < get_specs(mm)$err.obs.sigma*1.5, "DO.mod tracks DO.obs with not too much error")
  expect_true(rmse_DO(predict_DO(mm2)) < get_specs(mm2)$err.obs.sigma*1.5, "DO.mod tracks DO.obs with not too much error")
  # plot_DO_preds(predict_DO(mm))
  
  # predict_DO - DO.obs should be different each time unless seed is set. DO.mod should always be the same
  expect_true(!isTRUE(all.equal(predict_DO(mm)$DO.obs, predict_DO(mm)$DO.obs)))
  expect_true(isTRUE(all.equal(predict_DO(mm)$DO.mod, predict_DO(mm)$DO.mod)))
  mm <- metab_sim(data=dat, data_daily=select(dd, -DO.mod.1), 
                  specs=specs('s_np_oipcpi_eu_plrckm.rnorm', sim.seed=626, day_start=-1, day_end=23))
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
                  specs=specs('s_np_oipcpi_eu_plrckm.rnorm', err.obs.sigma=0, err.proc.sigma=0.5))
  DO_preds <- predict_DO(mm, date_start="2012-09-19")
  acf_out <- acf(DO_preds$DO.mod - DO_preds$DO.obs, plot=FALSE)
  expect_gt(acf_out$acf[acf_out$lag==1], 0.6)
  # plot_DO_preds(predict_DO(mm))
  
  # should be able to switch ODE methods in fitting
  dat <- select(data_metab('3', res='30'), -DO.obs)
  mmE <- metab_sim(specs('s_np_oipcpi_eu_plrckm.rnorm', err.obs.sigma=0, err.proc.sigma=0.05, sim.seed=4), data=dat, data_daily=dd)
  mmP <- metab_sim(specs('s_np_oipcpi_tr_plrckm.rnorm', err.obs.sigma=0, err.proc.sigma=0.05, sim.seed=4), data=dat, data_daily=dd)
  rmseEP <- sqrt(mean((predict_DO(mmE)$DO.obs - predict_DO(mmP)$DO.obs)^2, na.rm=TRUE))
  expect_gt(rmseEP, 0.001)
  expect_lt(rmseEP, 0.1)
  # DO_preds <- bind_rows(
  #   data.frame(predict_DO(mmE), method="euler", stringsAsFactors=FALSE),
  #   data.frame(predict_DO(mmP), method="trapezoid", stringsAsFactors=FALSE))
  # library(ggplot2)
  # ggplot(DO_preds, aes(x=solar.time, y=100*DO.mod/DO.sat, color=method)) + geom_line() + theme_bw()
  # ggplot(DO_preds, aes(x=solar.time, y=100*DO.obs/DO.sat, color=method)) + geom_line() + theme_bw()

})
