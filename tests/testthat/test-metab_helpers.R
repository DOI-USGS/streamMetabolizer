context("metab_model_helpers")

test_that("mm_data()", {
  expect_equal(mm_data(), get_data(metab_model()), info="mm_data() and metab_model()@data are the same")
  expect_equal(mm_data(date.time), mm_data()["date.time"], info="mm_data can be used for to select columns")
  expect_equal(mm_data(depth, temp.water, date.time), mm_data()[c("depth","temp.water","date.time")], info="mm_data can be used for to select columns")
})

test_that("mm_validate_data works", {
  expect_error(mm_validate_data(dplyr::select(mm_data(),-temp.water), "metab_mle"), "missing these columns: temp.water")
  expect_error(mm_validate_data(dplyr::mutate(mm_data(),temp.air=9), "metab_mle"), "should omit these extra columns: temp.air")
  expect_error(mm_validate_data(dplyr::mutate(mm_data(),temp.water=u(temp.water, "degF"), DO.obs=u(DO.obs, "mgO L^-1")), "metab_mle"), "unexpected units")
  expect_equal(mm_validate_data(mm_data(), "metab_mle"), mm_data())
})