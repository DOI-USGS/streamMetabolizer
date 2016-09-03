context('predict')

test_that('predict_DO works as expected', {
  
  # empty model
  mm <- metab_model()
  expect_null(predict_DO(mm))

})

test_that('predict_metab works on allmodel types', {
  dat <- data_metab('3','15')
  expected_cols <- c('date','GPP','GPP.lower','GPP.upper','ER','ER.lower','ER.upper','msgs.fit','warnings','errors')
  
  # empty model
  mm <- metab_model()
  expect_null(predict_metab(mm))
  
  # metab_mle
  mm <- metab_mle(data=dat)
  mp <- predict_metab(mm)
  expect_equal(names(mp), expected_cols)
  expect_equal(nrow(mp), 3)
  
  # metab_night
  mm <- metab_night(data=dat)
  mp <- predict_metab(mm)
  expect_equal(names(mp), expected_cols)
  expect_equal(nrow(mp), 3)
  
  # metab_bayes
  mm <- metab_bayes(specs(mm_name('bayes'), burnin_steps=100, saved_steps=100), data=dat) # breaks
  mp <- predict_metab(mm)
  expect_equal(names(mp), expected_cols)
  expect_equal(nrow(mp), 3)
  
  # metab_sim
  dat_daily <- data.frame(date=as.Date(paste0("2012-09-", 18:20)), GPP.daily=2, ER.daily=-3, K600.daily=21)
  mm <- metab_sim(specs(mm_name('sim')), data=dat, data_daily=dat_daily)
  mp <- predict_metab(mm)
  expect_equal(names(mp), expected_cols[!expected_cols %in% c('warnings','errors')])
  expect_equal(nrow(mp), 3)
  
})

test_that('predict_metab works as expected', {
  
  # should break reasonably for missing data, easy case for metab_mle which won't produce params on those days
  dat <- data_metab('3','15',flaws='missing end')
  mm <- metab_mle(data=dat)
  mp <- predict_metab(mm)
  expect_equal(mp[3,'GPP'], NA_real_)
  expect_equal(mp[3,'ER'], NA_real_)
  expect_equal(mp[3,'messages.fit'], "data don't start when expected")
  expect_equal(mp[3,'errors'], "NAs in DO.mod")
  
  # shouldn't warn but not break for nighttime regression
  dat <- data_metab('3', res='30', day_start=18, day_end=30)
  mm <- metab_night(specs(mm_name('night'), day_start=18, day_end=30), data=dat)
  mp <- predict_metab(mm)
  expect_equal() # should see warning
  
  # should break reasonably even for metab_sim, which won't have broken on model fitting
  dat <- data_metab('3','15',flaws='missing end')
  dat_daily <- data.frame(date=as.Date(paste0("2012-09-", 18:20)), GPP.daily=2, ER.daily=-3, K600.daily=21)
  mm <- metab_sim(specs(mm_name('sim')), data=dat, data_daily=dat_daily)
  mp <- predict_metab(mm)
  expect_equal(mp[3,'GPP'], NA_real_) # breaks for now
  expect_equal(mp[3,'ER'], NA_real_)
  expect_equal(mp[3,'errors'], "data don't start when expected")
  
  # should warn if specified day length isn't 24 hours
  dat <- data_metab('3','30', day_start=2)
  mm <- metab_mle(specs(mm_name('mle'), day_start=2), data=dat)
  predict_metab(mm)
  
})