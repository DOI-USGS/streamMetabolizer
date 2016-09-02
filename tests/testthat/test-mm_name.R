context('mm_name')

test_that("mm_name can generate names", {
  # missing args OK
  expect_equal(mm_name(), "m_np_oi_tr_plrckm.nlm")
  expect_equal(mm_name('bayes'), "b_np_oipi_tr_plrckm.stan")
  expect_equal(mm_name('m'), "m_np_oi_tr_plrckm.nlm") # even abbreviations work! lazy, though
  expect_equal(mm_name('sim'), "s_np_oipcpi_tr_plrckm.rnorm")
  expect_equal(mm_name('Kmodel'), "K_Kc___.lm")
  expect_equal(mm_name('b', pool_K600='none'), "b_np_oipi_tr_plrckm.stan")
  expect_equal(mm_name('b', pool_K600='none', err_proc_acor=TRUE), "b_np_oipcpi_tr_plrckm.stan")
  
  # catches bad arg combos
  expect_error(mm_name('b', pool_K600='none', err_proc_acor=TRUE, engine='nlm'), 'mismatch')
  expect_error(mm_name('m', err_proc_iid=TRUE), 'not among valid')
  expect_error(mm_name('s', err_proc_iid=FALSE), 'not among valid')
  expect_error(mm_name('n', ode_method='trapezoid'), 'not among valid')
})

test_that("mm_parse_name can parse names", {
  # parse a name
  expect_is(mm_parse_name("m_np_oi_tr_km.nlm"), "data.frame")
  expect_equal(dim(mm_parse_name("m_np_oi_tr_km.nlm")), c(1,10))
  expect_equal(mm_parse_name("n_np_pi_eu_rckf.lm")$ode_method, "euler")
  expect_equal(mm_parse_name("s_np_oipcpi_eu_plrckm.rnorm")$pool_K600, "none")
  expect_equal(mm_parse_name("b_Kl_oipcpi_eu_plrcko.rnorm")$pool_K600, "linear")
  expect_equal(mm_parse_name(mm_valid_names("Kmodel"))$engine, c('lm','mean','loess'))
  
  # parse and then rebuild a name
  expect_equal(do.call(mm_name, mm_parse_name("b_np_oipi_tr_plrckm.stan")), "b_np_oipi_tr_plrckm.stan")
  expect_equal(do.call(mm_name, mm_parse_name("K_Kc___.lm")[c('type','engine')]), "K_Kc___.lm")
})

test_that("mm_valid_names and mm_validate_names check model names", {
  # all the model names we know about are returned
  expect_lt(53, length(mm_valid_names()))
  
  # the models given by mm_valid_names() are all valid by mm_validate_name().
  # this is too slow to check completely now, so just pick a random sample
  nms <- sample(mm_valid_names(), size=10)
  validations <- lapply(nms, mm_validate_name)
  expect_true(all(sapply(validations, is.null)), info = paste0('validating ', paste(nms, collapse=', ')))
})

test_that("specs uses any valid mm_name", {
  # subsample because there are >1000 valid model names
  mnames <- mm_valid_names()
  expect_gte(length(mnames), 1000)
  specs_list <- lapply(sample(mnames, 50), specs)
  # specs list lengths differ by model type. the exact range could change
  sp_len_range <- range(sapply(specs_list, length))
  expect_gte(sp_len_range[1], 3)
  expect_lte(sp_len_range[2], 35)
})

