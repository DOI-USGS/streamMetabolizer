context("output of metabolism models")

test_that("metabolism models run & produce reasonable output", {
  
   # use a subset of data from Bob
   french <- streamMetabolizer:::load_french_broad()
   
   # metab_simple
   mm <- metab_simple(data=v(french))
   expect_is(mm, "metab_model") # should be a metab_simple soon
   expect_is(slot(mm, "fit"), "list") # specific to this model
   expect_equal(names(slot(mm, "fit")), c("minimum","estimate","gradient","code","iterations")) # specific to this model
   expect_is(slot(mm, "args"), "list")
   expect_is(slot(mm, "data"), "data.frame")
   expect_is(slot(mm, "pkg_version"), "character")
})