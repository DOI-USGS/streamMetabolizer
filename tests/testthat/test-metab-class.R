context("metab_model class and inheriting classes")

test_that("metab_model objects can be created and accessed", {
  
  # metab_model parent class
  mm <- metab_model()
  expect_output(slot(mm, "fit"), "generic")
  expect_is(slot(mm, "args"), "list")
  expect_is(slot(mm, "data"), "data.frame")
  expect_is(slot(mm, "pkg_version"), "character")

  # display
  expect_output(mm, "metab_model")
  
  # access
  get_fit(mm)
  get_args(mm)
  get_data(mm)
  get_version(mm)
})

test_that("metab_model objects generate predictions", {
  # metab_model parent class
  mm <- metab_model()
  
  # predict
  #predict_metab(mm)
  #predict_DO(mm)
})