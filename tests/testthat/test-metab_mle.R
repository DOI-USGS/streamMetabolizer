context("metab_mle")

# comparing old method (manual euler/pairmeans) to new (ode() solutions + dDOdt analytical pairmeans)
manual_temporary_test <- function() {
    
  # temporary b/c after checking the transition, i want to remove the 'old'
  # mle_method from the package
  
  # get a big enough dataset to take time for mle. this is why the test is manual
  dat <- mda.streams::get_metab_data('nwis_08062500', start_date='2011-04-01', end_date='2011-07-01')
  expect_equal(dim(dat), c(8700, 6))
  
  # debug
  dat <- data_metab('10','30')[(48*6+1):(48*7),]
  sp <- specs(mm_name('mle', ode_method='rk4'))
  mm.debug.old <- metab(replace(sp, 'mle_method', 'old'), dat)
  mm.debug.new <- metab(replace(sp, 'mle_method', 'new'), dat)
  
  # compare Eulers
  dat <- data_metab('10','30')
  sp <- specs(mm_name('mle', ode_method='Euler'))
  mm.debug.old <- metab(replace(sp, 'mle_method', 'old'), dat)
  mm.debug.new <- metab(replace(sp, 'mle_method', 'new'), dat)
  plot_metab_preds(mm.debug.old)
  plot_metab_preds(mm.debug.new)
  plot_DO_preds(mm.debug.old)
  plot_DO_preds(mm.debug.new)
  predict_metab(mm.debug.new)[c('GPP','ER','K600')]-predict_metab(mm.debug.old)[c('GPP','ER','K600')]
  predict_metab(mm.debug.new)[c('GPP','ER','K600')]/predict_metab(mm.debug.old)[c('GPP','ER','K600')]
  
  # test with old & new implementations
  dat <- data_metab('3','30')
  mm <- list(rk2=list(new=NA), rk4=list(new=NA), lsoda=list(new=NA),
             Euler=list(old=NA, new=NA), pairmeans=list(old=NA, new=NA))
  for(ode_method in names(mm)) {
    message(ode_method)
    sp <- specs(mm_name('mle', ode_method=ode_method))
    for(mle_method in names(mm[[ode_method]])) {
      mm[[ode_method]][[mle_method]] <- metab(replace(sp, 'mle_method', mle_method), dat)
    }
  }
  
  # check that we haven't lost much runtime efficiency with the new 
  # implementation relative to the old; find out how much longer the more
  # complex methods take
  fitting.times <- lapply(mm, function(mm_list) {
    lapply(mm_list, function(mm_obj) {
      get_fitting_time(mm_obj)[['elapsed']]
    })
  })
  
  # check that intercepts are 0 and slopes are 1 and correlations are perfect
  # between old and new implementations of Euler and pairmeans
  old.new.lms <- lapply(mm, function(mm_list) {
    if(all(c('old','new') %in% names(mm_list))) {
      preds.old <- predict_metab(mm_list$old)
      preds.new <- predict_metab(mm_list$new)
      sapply(c("GPP","GPP.lower","GPP.upper","ER","ER.lower","ER.upper","K600","K600.lower","K600.upper"), function(var) {
        onlm <- lm(preds.new[[var]] ~ preds.old[[var]])
        c(coef(onlm), r.squared=suppressWarnings(summary(onlm))$r.squared, 
          cor=cor(preds.new[[var]], preds.old[[var]], use="complete.obs"))
      })
    }
  })
  
  # compare predictions between ODE methods
  preds <- lapply(mm, function(mm_list) {
    predict_metab(mm_list$new)
  })
  pred.cors <- lapply(setNames(nm=c("GPP","GPP.lower","GPP.upper","ER","ER.lower","ER.upper","K600","K600.lower","K600.upper")), function(var) {
    cor(sapply(preds, function(p) p[[var]]))
  })
  rrmse <- function(mat) {
    rmses <- sapply(1:ncol(mat), function(c1) {
      sapply(1:ncol(mat), function(c2) {
        sqrt(mean(((mat[,c1] - mat[,c2])/mat[,c1])^2, na.rm=TRUE))
      })
    })
    colnames(rmses) <- rownames(rmses) <- colnames(mat)
    rmses
  }
  pred.rmse <- lapply(setNames(nm=c("GPP","GPP.lower","GPP.upper","ER","ER.lower","ER.upper","K600","K600.lower","K600.upper")), function(var) {
    rrmse(sapply(preds, function(p) p[[var]]))
  })
  
}


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
  mmE <- metab_mle(specs('m_np_oi_eu_plrckm.nlm'), data=dat3)
  mmP <- metab_mle(specs('m_np_oi_pm_plrckm.nlm'), data=dat3)
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
