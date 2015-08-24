context("metab_night")

french <- streamMetabolizer:::load_french_creek()

test_that("metab_night models can be created", {
  
  mm <- metab_night(data=v(french))
  
  # check basic structure
  expect_is(mm, "metab_night")
  expect_is(slot(mm, "fit"), "data.frame")
  expect_is(slot(mm, "args"), "list")
  expect_is(slot(mm, "data"), "data.frame")
  expect_is(slot(mm, "pkg_version"), "character")
  
})