context("metab_Kmodel")

library(dplyr)

test_that("metab_Kmodel predictions (predict_metab, predict_DO) make sense", {
    
  # make up some data.frames of the right sizes
  set.seed(4438)
  # dat = data.frame of instantaneous solar.time and discharge
  dat <- data_metab('3') %>%
    select(solar.time) %>% 
    mutate(discharge=exp(rnorm(length(solar.time))))
  # ddat = data.frame of daily date, discharge.daily, and K600 (made-up data)
  ddat <- data.frame(date=seq(as.Date("2012-08-15"),as.Date("2012-09-15"),
    as.difftime(1,units='days')), discharge.daily=exp(rnorm(32,2,1)), K600=rnorm(32,30,4)) %>%
    mutate(K600.lower=K600-5, K600.upper=K600+6)
  # ddat1 = data.frame of just one day of data
  ddat1 <- data.frame(date=as.Date("2012-08-24"), K600=20, K600.lower=15, K600.upper=25, stringsAsFactors=FALSE)
  
  # 1-day tests: show that Kmodel(mean) won't break even if only 1 data point is entered
  expect_warning(mm <- metab_Kmodel(data=NULL, data_daily=ddat1, specs=specs(mm_name("Kmodel", engine='mean'))), "no SE available")
  expect_equal(predict_metab(mm)$K600, ddat1$K600)
  # show that Kmodel(lm) and Kmodel(loess) do break, on model-specific errors
  expect_error(metab_Kmodel(data=dat, data_daily=ddat1, specs=specs(mm_name("Kmodel", engine='lm'))), "0 (non-NA) cases", fixed=TRUE)
  expect_error(metab_Kmodel(data=dat, data_daily=ddat1, specs=specs(mm_name("Kmodel", engine='loess'))), "invalid 'x'", fixed=TRUE)
  
  # mean
  expect_warning(mm_mean <- metab_Kmodel(data_daily=ddat, specs=specs(mm_name('Kmodel', engine='mean'))), "no SE available")
  
  # lm
  ms <- specs(mm_name('Kmodel', engine='lm'), predictors=c())
  mm <- metab_Kmodel(data_daily=ddat, specs=ms)
  expect_equal(predict_metab(mm_mean)$K600, predict_metab(mm)$K600)
  mm <- metab_Kmodel(replace(ms, 'predictors', 'discharge.daily'), data_daily=ddat)
  mm <- metab_Kmodel(replace(ms, 'predictors', 'date'), data_daily=ddat)
  mm <- metab_Kmodel(specs(mm_name('Kmodel', engine='lm'), predictors=c('discharge.daily','date')), data_daily=ddat)
  
  # loess
  mm <- metab_Kmodel(specs(mm_name('Kmodel', engine='loess'), predictors='date'), data_daily=ddat)
  mm <- metab_Kmodel(specs(mm_name('Kmodel', engine='loess'), predictors='date', other_args=list(span=1.4)), data_daily=ddat)
  mm <- metab_Kmodel(specs(mm_name('Kmodel', engine='loess'), predictors=c('discharge.daily','date')), data_daily=ddat)
  mm <- metab_Kmodel(specs(mm_name('Kmodel', engine='loess'), predictors='discharge.daily', other_args=list(span=2)), data_daily=ddat)
  
  # use graphical check to inspect above mm results
  # plot_metab_preds(predict_metab(mm))
  expect_error(predict_DO(mm), "can only predict K, not DO, from metab_Kmodel")
    
})

test_that("try a complete PRK-K-PR workflow", {
    
  dat <- data_metab('10', res='10', day_start=-1, day_end=25)
  
  # fit a first-round MLE and extract the K estimates
  mm1 <- metab(specs(mm_name('mle'), day_start=-1, day_end=25), data=dat)
  K600_mm1 <- predict_metab(mm1) %>% select(date, K600, K600.lower, K600.upper)
  
  # smooth the K600s
  expect_warning({
    mm2 <- metab(specs(mm_name("Kmodel", engine='mean'), transforms=c(K600='log'), weights=c(), predictors=c()), data_daily=K600_mm1)
  }, "no SE available")
  K600_mm2 <- predict_metab(mm2) %>% select(date, K600)
  
  # refit the MLE with fixed K
  mm3 <- metab(specs(mm_name('mle'), day_start=-1, day_end=25), data=dat, data_daily=K600_mm2)
  expect_equal(length(unique(predict_metab(mm3)$K600)), 1)
  
  # that should have reduced variance in the GPP & ER predictions
  expect_lt(sd(predict_metab(mm3)$GPP), sd(predict_metab(mm1)$GPP))
  expect_lt(sd(predict_metab(mm3)$ER), sd(predict_metab(mm1)$ER))
  
})