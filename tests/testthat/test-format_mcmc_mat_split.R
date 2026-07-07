test_that("metab_bayes(split_dates=TRUE) runs end to end and produces the expected predict_metab() shape", {
  skip_on_cran()

  # split_dates=TRUE requires a pooled K600 model, which in turn requires a
  # discharge column in the input data.
  dat <- data_metab('3', '30')
  dat$discharge <- 3
  sp <- suppressWarnings(specs(
    'b_Kl_oi_eu_plrckm.stan', split_dates=TRUE, n_chains=2, n_cores=2,
    burnin_steps=200, saved_steps=100))

  # 1. metab_bayes() completes without error
  mm <- metab_bayes(specs=sp, data=dat)

  # 2. predict_metab(mm) returns a data frame with more than 0 rows
  mp <- predict_metab(mm)
  expect_true(is.data.frame(mp))
  expect_gt(nrow(mp), 0)

  # 3. column names for split_dates=TRUE output: predict_metab() produces the
  # same standardized column set regardless of split_dates, since it's built
  # from get_fit(mm)$daily via the same predict_metab.metab_bayes code path
  expect_true(all(c('date','GPP','GPP.lower','GPP.upper','ER','ER.lower','ER.upper') %in% names(mp)))

  # 4. date ordering/structure: 3 consecutive, strictly increasing days with
  # 1-day gaps, same as the nosplit case
  expect_equal(nrow(mp), 3)
  expect_equal(mp$date, as.Date(c('2012-09-18','2012-09-19','2012-09-20')))
  expect_true(all(diff(mp$date) > 0))
  expect_equal(as.numeric(diff(mp$date)), c(1, 1))
})
