context('mm_name')

test_that("mm_name can generate names and mm_parse_name can parse them", {
  # missing args OK
  expect_equal(mm_name(), "b_np_oipi_pm_km.stan")
  expect_equal(mm_name('bayes'), "b_np_oipi_pm_km.stan")
  expect_equal(mm_name('m'), "m_np_oi_pm_km.nlm") # even abbreviations work! lazy, though
  expect_equal(mm_name('sim'), "s_np_oipcpi_eu_.rnorm")
  expect_equal(mm_name(pooling='none'), "b_np_oipi_pm_km.stan")
  expect_equal(mm_name('b', pooling='none', err_proc_acor=TRUE), "b_np_oipcpi_pm_km.stan")
  
  # catches bad arg combos
  expect_error(mm_name('b', pooling='none', err_proc_acor=TRUE, bayes_software='nlm'), 'mismatch')
  expect_error(mm_name('m', err_proc_iid=TRUE), 'not among valid')
  expect_error(mm_name('s', err_proc_iid=FALSE), 'not among valid')
  expect_error(mm_name('n', ode_method='pairmeans'), 'not among valid')
})

test_that("mm_valid_names and mm_validate_names check model names", {
  # all the model names we know about are returned
  expect_equal(length(mm_valid_names()), 51)
  
  # the models given by mm_valid_names() are all valid by mm_validate_name()
  validations <- sapply(mm_valid_names(), mm_validate_name)
  expect_true(all(sapply(validations, is.null)))
})

test_that("specs uses any valid mm_name", {
  specs_list <- lapply(mm_valid_names(), specs)
  # specs list lengths differ by model type. the exact range could change
  expect_equal(range(sapply(tryme, length)), c(1,27))
})

