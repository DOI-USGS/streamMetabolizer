context("metab_model class and inheriting classes")


test_that("metab_model objects can be created and accessed", {
  
  # basic structure of metab_model parent class
  mm <- metab_model()
  expect_output(slot(mm, "fit"), "generic")
  expect_is(slot(mm, "args"), "list")
  expect_is(slot(mm, "data"), "data.frame")
  expect_is(slot(mm, "pkg_version"), "character")

  # display
  expect_output(mm, "metab_model")
  
  # accessors
  expect_equal(slot(mm, "fit"), get_fit(mm))
  expect_equal(slot(mm, "args"), get_args(mm))
  expect_equal(slot(mm, "data"), get_data(mm))
  expect_equal(slot(mm, "pkg_version"), get_version(mm))
  expect_equal(slot(mm, "info"), get_info(mm))
})


test_that("metab_models can be saved", {
  
  # save and reload
  mm <- metab_model()
  mm@fit
  save(mm, file="test_temp.RData")
  mm_orig <- mm
  load("test_temp.RData")
  file.remove("test_temp.RData")
  expect_equal(mm_orig, mm)
  expect_equal(get_fit(mm_orig), get_fit(mm))
  expect_equal(get_data(mm_orig), get_data(mm))
  
})

test_that("metab_models have default predict_metab and predict_DO methods", {
  
  mm <- metab_model()
  
  # predict_metab
  expect_warning(metab <- predict_metab(mm), "model is missing columns")
  expect_equal(names(metab), c("local.date","GPP","GPP.lower","GPP.upper","ER","ER.lower","ER.upper","K600","K600.lower","K600.upper"))
  
  # predict_DO
  expect_warning(DO_preds <- predict_DO(mm), "model is missing columns")
  expect_equal(names(DO_preds), c("local.date", names(get_data(mm)), "DO.mod"))

})