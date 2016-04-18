context("metab_model class and inheriting classes")

test_that("metab_model objects can be created and accessed", {
  
  # basic structure of metab_model parent class
  mm <- metab_model()
  expect_output(print(slot(mm, "fit")), "generic metab_model class; no actual fit")
  expect_true(all(names(formals(metab)) %in% names(getSlots('metab_model'))), info="slots should match args to metab()")
  expect_is(slot(mm, "fitting_time"), "proc_time", info="time should be recorded")
  expect_is(slot(mm, "data"), "data.frame", info="default should populate with example data")
  expect_is(slot(mm, "pkg_version"), "character", info="pkg version should be autopopulated")

  # display
  expect_output(print(mm), "metab_model", info="model type should be shown")
  
  # accessors
  expect_equal(slot(mm, "fit"), get_fit(mm))
  expect_equal(slot(mm, "specs"), get_specs(mm))
  expect_equal(slot(mm, "data"), get_data(mm))
  expect_equal(slot(mm, "pkg_version"), get_version(mm))
  expect_equal(slot(mm, "info"), get_info(mm))
})

test_that("metab_models have default predict_metab and predict_DO methods", {
  
  mm <- metab_model()
  
  # predict_metab
  expect_warning(metab <- predict_metab(mm), "model is missing columns")
  expect_equal(names(metab), c("date","GPP","GPP.lower","GPP.upper","ER","ER.lower","ER.upper","K600","K600.lower","K600.upper"))
  
  # can't predict_DO
  expect_warning(expect_error(DO_preds <- predict_DO(mm), "day_start must be specified"), "missing columns for estimates")

})

test_that("metab_models can be saved & reloaded (see helper-save_load.R)", {
  
  # save and reload
  mm <- metab_model()
  
  # see if saveRDS with gzfile, compression=9 works well
  rdstimes <- save_load_timing(mm, reps=1) # autoloaded b/c script begins with 'helper' and is in this directory
  expect_true('gz6' %in% rdstimes$typelevel[1:3], info="gz6 is reasonably efficient for saveRDS")
  # plot_save_load_timing(rdstimes)
  
  # save and load the mm, make sure it stays the same
  test_save_load_recovery(mm)
  
})
