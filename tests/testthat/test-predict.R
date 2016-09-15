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
  mm <- metab_night(data=data_metab('3', day_start=12, day_end=36))
  mp <- predict_metab(mm)
  expect_equal(names(mp), expected_cols)
  expect_equal(nrow(mp), 3)
  
  # metab_bayes
  mm <- metab_bayes(specs(mm_name('bayes'), burnin_steps=100, saved_steps=100), data=dat)
  mp <- predict_metab(mm)
  expect_equal(names(mp), expected_cols)
  expect_equal(nrow(mp), 3)
  
  # metab_sim
  dat_daily <- data.frame(date=as.Date(paste0("2012-09-", 18:20)), GPP.daily=2, ER.daily=-3, K600.daily=21)
  mm <- metab_sim(specs(mm_name('sim')), data=dat, data_daily=dat_daily)
  mp <- predict_metab(mm)
  expect_equal(names(mp), expected_cols)
  expect_equal(nrow(mp), 3)
  
})

test_that('predict_metab works as expected for bad inputs', {
  
  # should stop on fitting for missing data and/or fitted parameters
  dat <- data_metab('3','15',flaws='missing end')
  # don't bother predicting on days where we didn't get a model fit
  mm <- metab_mle(data=dat)
  mp <- predict_metab(mm)
  expect_equal(mp[3,'GPP'], NA_real_)
  expect_equal(mp[3,'ER'], NA_real_)
  expect_equal(mp[3,'msgs.fit'], '      E')
  expect_equal(mp[3,'warnings'], NA_character_)
  expect_equal(get_params(mm)[3,'errors'], "data don't start when expected")
  # notice bad days for metab_sim, which won't have broken on model fitting
  dat_daily <- data.frame(date=as.Date(paste0("2012-09-", 18:20)), GPP.daily=2, ER.daily=-3, K600.daily=21)
  mm <- metab_sim(specs(mm_name('sim')), data=dat, data_daily=dat_daily)
  mp <- predict_metab(mm, use_saved=FALSE)
  expect_true(is.na(mp[3,'GPP']) && is.na(mp[3,'ER']))
  expect_true(is.na(mp[3,'msgs.fit']))
  expect_equal(mp[3,'errors'], "data don't start when expected")
  
  # should NOT stop on fitting if we said not to test
  mm <- metab_mle(specs(mm_name('mle'), day_tests=c()), data=dat)
  mp <- predict_metab(mm)
  expect_true(all(mp$msgs.fit == '       '))
  expect_true(all(mp$warnings == ''))
  expect_true(all(mp$errors == ''))
  mm <- metab_sim(specs(mm_name('sim'), day_tests=c()), data=dat, data_daily=dat_daily)
  mp <- predict_metab(mm)
  expect_equal(mp[3,'GPP'], get_params(mm)[3,'GPP.daily'])
  expect_equal(mp[3,'ER'], get_params(mm)[3,'ER.daily'])
  expect_true(all(is.na(mp$msgs.fit)))
  expect_true(all(mp$warnings == ''))
  expect_true(all(mp$errors == ''))
  
  # should give message and force day length to 24 hours for prediction
  dat <- data_metab('3','30', day_start=2)
  expect_message(mm <- metab_mle(specs(mm_name('mle'), day_start=2), data=dat), "differs from the model-fitting range")
  mp <- predict_metab(mm)
  expect_true(all(!is.na(mp$GPP)))
  expect_true(all(mp$msgs.fit == '       '))
  expect_true(all(mp$warnings == ''))
  expect_true(all(mp$warnings == ''))
  expect_message(expect_error(predict_metab(mm, day_start=20, day_end=28, use_saved=FALSE), 'day_end - day_start < 24 hours'), "differs from the model-fitting range")
  expect_error(predict_metab(mm, day_start=2, day_end=28, use_saved=FALSE), 'day_end - day_start must not exceed 24 hours')
  # same for metab_night when requested day is too long
  dat <- data_metab('3','30', day_start=7, day_end=36)
  expect_message(mm <- metab(specs(mm_name('night'), day_start=7), data=dat), "differs from the model-fitting range")
  expect_error(mp <- predict_metab(mm, day_start=7, day_end=36, use_saved=FALSE), 'day_end - day_start must not exceed 24 hours')
  # but should allow <24 hours if needed for nighttime regression
  expect_message(mm <- metab(specs(mm_name('night'), day_start=12, day_end=24), data=dat), "GPP estimates are 0 because they're for nighttime only")
  mp <- predict_metab(mm)
  expect_true(all(!is.na(mp$ER)))
  expect_true(all(mp$msgs.fit == '       '))
  expect_true(all(mp$warnings == ''))
  expect_true(all(mp$warnings == ''))
  # and should notice if we're ignoring the day_start and day_end values
  expect_warning(mp <- predict_metab(mm, day_start=7, day_end=36), 'using saved daily metabolism values')
})
