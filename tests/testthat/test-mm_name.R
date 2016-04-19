context('mm_name')

test_that("mm_name can generate names", {
  # missing args OK
  expect_equal(mm_name(), "m_np_oi_pm_km.nlm")
  expect_equal(mm_name('bayes'), "b_np_oipi_pm_km.stan")
  expect_equal(mm_name('m'), "m_np_oi_pm_km.nlm") # even abbreviations work! lazy, though
  expect_equal(mm_name('sim'), "s_np_oipcpi_pm_.rnorm")
  expect_equal(mm_name('Kmodel'), "K_np___.lm")
  expect_equal(mm_name('b', pool_K600='none'), "b_np_oipi_pm_km.stan")
  expect_equal(mm_name('b', pool_K600='none', err_proc_acor=TRUE), "b_np_oipcpi_pm_km.stan")
  
  # catches bad arg combos
  expect_error(mm_name('b', pool_K600='none', err_proc_acor=TRUE, engine='nlm'), 'mismatch')
  expect_error(mm_name('m', err_proc_iid=TRUE), 'not among valid')
  expect_error(mm_name('s', err_proc_iid=FALSE), 'not among valid')
  expect_error(mm_name('n', ode_method='pairmeans'), 'not among valid')
})

test_that("mm_parse_name can parse names", {
  # parse a name
  expect_is(mm_parse_name("m_np_oi_pm_km.nlm"), "data.frame")
  expect_equal(dim(mm_parse_name("m_np_oi_pm_km.nlm")), c(1,8))
  expect_equal(mm_parse_name("n_np_pi_eu_.lm")$ode_method, "Euler")
  expect_equal(mm_parse_name("s_np_oipcpi_eu_.rnorm")$pool_K600, "none")
  expect_equal(mm_parse_name("b_Kl_oipcpi_eu_.rnorm")$pool_K600, "linear")
  expect_equal(mm_parse_name(mm_valid_names("Kmodel"))$engine, c('lm','mean','loess'))
  
  # parse and then rebuild a name
  expect_equal(do.call(mm_name, mm_parse_name("b_np_oipi_pm_km.stan")), "b_np_oipi_pm_km.stan")
  expect_equal(do.call(mm_name, mm_parse_name("K_np___.lm")[c('type','engine')]), "K_np___.lm")
})

test_that("mm_valid_names and mm_validate_names check model names", {
  # all the model names we know about are returned
  expect_lt(53, length(mm_valid_names()))
  
  # the models given by mm_valid_names() are all valid by mm_validate_name()
  validations <- sapply(mm_valid_names(), mm_validate_name)
  expect_true(all(sapply(validations, is.null)))
})

test_that("specs uses any valid mm_name", {
  specs_list <- lapply(mm_valid_names(), specs)
  # specs list lengths differ by model type. the exact range could change
  expect_equal(range(sapply(specs_list, length)), c(4,34))
})

