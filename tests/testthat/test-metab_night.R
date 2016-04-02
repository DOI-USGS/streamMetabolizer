context("metab_night")

test_that("metab_night models can be created", {
  
  mm <- metab_night(data=data_metab('1', res='30', day_start=12, day_end=36))
  
  # check basic metab_night
  expect_is(mm, "metab_night")
  expect_is(slot(mm, "fit"), "data.frame")
  expect_is(slot(mm, "specs"), "list")
  expect_is(slot(mm, "data"), "data.frame")
  expect_is(slot(mm, "pkg_version"), "character")
  
})

test_that("metab_night predictions (predict_metab, predict_DO) make sense", {
  
  # 1 day
  mm <- metab_night(data=data_metab('1', res='30', day_start=12, day_end=36))
  metab <- predict_metab(mm)
  DO_preds <- predict_DO(mm)
  expect_true(rmse_DO(DO_preds) < 0.1, info="DO.mod tracks DO.obs with not too much error")
  # plot_DO_preds(DO_preds)
  
  # 10 days
  mm <- metab_night(data=data_metab('10', res='5', day_start=12, day_end=36))
  metab <- predict_metab(mm)
  DO_preds <- predict_DO(mm)
  expect_true(rmse_DO(DO_preds) < 0.1, info="DO.mod tracks DO.obs with not too much error")
  # plot_DO_preds(DO_preds)
  
})

test_that("day_tests=c('full_day','include_sunset') get handled appropriately", {
  
  # specs should include both tests by default
  expect_true(all(c('full_day','include_sunset') %in% specs(mm_name('night'))$day_tests))
  
  # when the date bounds match & span the full night, definitely no errors
  mm <- metab_night(specs(mm_name('night'), day_start=12, day_end=35), data=data_metab('1', day_start=12, day_end=35))
  expect_equal(predict_metab(mm)$errors, "")
  
  # error if day starts after day_start & dusk
  sp <- specs(mm_name('night'), day_start=12, day_end=35, day_tests=c('include_sunset','full_day'))
  dat <- data_metab('1',        day_start=22, day_end=35)
  expect_equal(predict_metab(metab_night(replace(sp, 'day_tests', 'full_day'), data=dat))$errors, "data don't start when expected")
  expect_equal(predict_metab(metab_night(replace(sp, 'day_tests', 'include_sunset'), data=dat))$errors, "data don't include day-night transition")
  expect_equal(predict_metab(metab_night(sp, data=dat))$errors, "data don't include day-night transition; data don't start when expected")
  # plot_DO_preds(predict_DO(metab_night(replace(sp, 'day_tests', c()), data=dat)))
  
  # full_day is OK, include_sunset isn't if day starts after dusk but before/on day_start
  sp <- specs(mm_name('night'), day_start=20, day_end=35, day_tests=c('include_sunset','full_day'))
  dat <- data_metab('1',        day_start=19, day_end=35)
  expect_equal(predict_metab(metab_night(replace(sp, 'day_tests', 'full_day'), data=dat))$errors, "")
  expect_equal(predict_metab(metab_night(replace(sp, 'day_tests', 'include_sunset'), data=dat))$errors, "data don't include day-night transition")
  # plot_DO_preds(predict_DO(metab_night(replace(sp, 'day_tests', 'full_day'), data=dat)))
  
  # full_day & include_sunset are OK if day starts after day_start but before/on dusk
  sp <- specs(mm_name('night'), day_start=12, day_end=35, day_tests=c('include_sunset','full_day'))
  dat <- data_metab('1',        day_start=14, day_end=35)
  expect_equal(predict_metab(metab_night(sp, data=dat))$errors, "")
  # plot_DO_preds(predict_DO(metab_night(sp, data=dat)))
  
  # full_day is OK if day ends before dawn but exactly on day_end
  sp <- specs(mm_name('night'), day_start=12, day_end=23)
  dat <- data_metab('1',        day_start=12, day_end=23)
  expect_equal(predict_metab(metab_night(sp, data=dat))$errors, "")
  
  # full_day is not OK if day ends before dawn & before day_end
  sp <- specs(mm_name('night'), day_start=12, day_end=36)
  dat <- data_metab('1',        day_start=12, day_end=24)
  expect_equal(predict_metab(metab_night(sp, data=dat))$errors, "data don't end when expected")
  
  # full_day is great if day ends after dawn, no matter whether it ends before or on day_end
  sp <- specs(mm_name('night'), day_start=12, day_end=36)
  dat <- data_metab('1',        day_start=12, day_end=32)
  expect_equal(predict_metab(metab_night(sp, data=dat))$errors, "")
  
})

test_that("metab_night predictions can be passed back into metab_mle", {
  
  # metab_night
  kdat <- data_metab('10', day_start=12, day_end=36)
  mmk <- metab_night(data=kdat)
  # plot_metab_preds(predict_metab(mmk))
  
  # metab_mle
  mledat <- data_metab('10')
  mm <- metab_mle(data=mledat, data_daily=predict_metab(mmk)[c('date', 'K600')])
  expect_equal(predict_metab(mm)$K600, predict_metab(mmk)$K600)
  expect_true(rmse_DO(predict_DO(mm)) < 0.2, info="DO.mod tracks DO.obs with not too much error")
  # plot_DO_preds(predict_DO(mm))
  
})

test_that("metab_night predictions match Bob's", {
  
  dat <- data_metab('1', day_start=12, day_end=36)
  mms <- metab_night(
    specs('n_np_pi_eu_.lm', day_start=18, day_end=29, day_tests=c('full_day','even_timesteps','complete_data')), data=dat)
  mmb <- streamMetabolizer:::load_french_creek_std_mle(
    dat, estimate='K', start=c(dates="09/18/12", times="18:00:00"), end=c(dates="09/19/12", times="05:00:00"))
  expect_equal(predict_metab(mms)$K600, mmb$K)
})
