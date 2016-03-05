context("metab_Kmodel")

library(dplyr)
vfrench <- streamMetabolizer:::load_french_creek(attach.units=FALSE)

test_that("metab_Kmodel predictions (predict_metab, predict_DO) make sense", {
    
  # fr = data.frame of instantaneous solar.time and discharge
  fr <- streamMetabolizer:::load_french_creek() %>% unitted::v() %>%
    filter(format(solar.time, "%Y-%m-%d") == "2012-08-24") %>%
    select(solar.time) %>% 
    mutate(discharge=exp(rnorm(length(solar.time))))
  # frk2 = data.frame of daily date, discharge.daily, and K600 (made-up data)
  frk2 <- data.frame(date=seq(as.Date("2012-08-15"),as.Date("2012-09-15"),
    as.difftime(1,units='days')), discharge.daily=exp(rnorm(32,2,1)), K600=rnorm(32,30,4)) %>%
    mutate(K600.lower=K600-5, K600.upper=K600+6)
  # frk = data.frame of just one day of data
  frk <- data.frame(date=as.Date("2012-08-24"), K600=20, K600.lower=15, K600.upper=25, stringsAsFactors=FALSE)
  
  # 1-day tests: show that Kmodel(mean) won't break even if only 1 data point is entered
  expect_warning(mm <- metab_Kmodel(data=NULL, data_daily=frk, specs=specs(mm_name("Kmodel", engine='mean'))), "no SE available")
  expect_equal(predict_metab(mm)$K600, frk$K600)
  # show that Kmodel(lm) and Kmodel(loess) do break, on model-specific errors
  expect_error(metab_Kmodel(data=fr, data_daily=frk, specs=specs(mm_name("Kmodel", engine='lm'))), "Error in lm.wfit", fixed=TRUE)
  expect_error(metab_Kmodel(data=fr, data_daily=frk, specs=specs(mm_name("Kmodel", engine='loess'))), "Error in simpleLoess", fixed=TRUE)
  
  # mean
  expect_warning(mm_mean <- metab_Kmodel(data_daily=frk2, specs=specs(mm_name('Kmodel', engine='mean'))), "no SE available")
  
  # lm
  ms <- specs(mm_name('Kmodel', engine='lm'), predictors=c())
  mm <- metab_Kmodel(data_daily=frk2, specs=ms)
  expect_equal(predict_metab(mm_mean)$K600, predict_metab(mm)$K600)
  mm <- metab_Kmodel(data_daily=frk2, specs=replace(ms, 'predictors', 'discharge.daily'))
  mm <- metab_Kmodel(data_daily=frk2, specs=replace(ms, 'predictors', 'date'))
  mm <- metab_Kmodel(data_daily=frk2, specs=specs(mm_name('Kmodel', engine='lm'), predictors=c('discharge.daily','date')))
  
  # loess
  mm <- metab_Kmodel(data_daily=frk2, specs=specs(mm_name('Kmodel', engine='loess'), predictors='date'))
  mm <- metab_Kmodel(data_daily=frk2, specs=specs(mm_name('Kmodel', engine='loess'), predictors='date', other_args=list(span=1.4)))
  mm <- metab_Kmodel(data_daily=frk2, specs=specs(mm_name('Kmodel', engine='loess'), predictors=c('discharge.daily','date')))
  mm <- metab_Kmodel(data_daily=frk2, specs=specs(mm_name('Kmodel', engine='loess'), predictors='discharge.daily', other_args=list(span=2)))
  
  # use graphical check to inspect above mm results
  # plot_metab_preds(predict_metab(mm))
  expect_error(predict_DO(mm), "can only predict K, not DO, from metab_Kmodel")
    
})

test_that("try a complete PRK-K-PR workflow", {
    
  # fit a first-round MLE and extract the K estimates
  expect_warning({mm1 <- metab_mle(data=vfrench, day_start=-1, day_end=23)}, "temperature out of range")
  K600_mm1 <- predict_metab(mm1) %>% select(date, K600, K600.lower, K600.upper)
  
  # smooth the K600s
  expect_warning({mm2 <- metab_Kmodel(data_daily=K600_mm1, specs=specs(
    mm_name("Kmodel", engine='mean'), transforms=c(K600='log'), weights=c(), predictors=c()))}, "no SE available")
  K600_mm2 <- predict_metab(mm2) %>% select(date, K600)
  
  # refit the MLE with fixed K
  expect_warning({mm3 <- metab_mle(data=vfrench, data_daily=K600_mm2, day_start=-1, day_end=23)}, "temperature out of range")
  predict_metab(mm3)
  
  # compare the two MLE results
  # library(ggplot2)
  # library(tidyr)
  # preds <- bind_rows(
  #   mutate(predict_metab(mm1), model='PRK'),
  #   mutate(predict_metab(mm3), model='PR')) %>%
  #   filter(!is.na(K600)) %>%
  #   gather(estimate, value, GPP, ER, K600)
  # ggplot(preds, aes(x=date, y=value, color=model)) + geom_point() + geom_line() + theme_bw() + facet_grid(estimate ~ ., scales = "free_y")

})