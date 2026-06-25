
test_that("plot_metab_preds retains all groups when fit data are present", {
  mm <- metab_night(specs(mm_name('night')), data=data_metab('3', day_start=12, day_end=36))
  g <- plot_metab_preds(mm)

  expect_s3_class(g, "ggplot")

  preds <- g$data
  # Column presence checked with %in% rather than exact order — group_modify()
  # moves the grouping key to first position, unlike the old do() code.
  expect_true(all(c('date', 'as', 'fit', 'lwr', 'upr', 'var') %in% names(preds)))

  # 3 dates × 2 groups (GPP, ER) = 6 rows; both groups retained
  expect_equal(nrow(preds), 6)
  expect_setequal(unique(preds$as), c('GPP', 'ER'))
})

test_that("plot_metab_preds drops groups where all fit values are NA", {
  # Verifies the group_modify() replacement for do() in plot_metab_preds.R.
  # plot_metab_preds() accepts a predict_metab() data frame directly.
  mm <- metab_night(specs(mm_name('night')), data=data_metab('3', day_start=12, day_end=36))
  mp <- predict_metab(mm)

  # Force GPP to all-NA so the GPP group should be dropped
  mp_no_gpp <- mp
  mp_no_gpp$GPP <- NA_real_

  g <- plot_metab_preds(mp_no_gpp)
  preds <- g$data

  # GPP group dropped; only ER rows remain (3 dates)
  expect_false('GPP' %in% unique(preds$as))
  expect_true('ER' %in% unique(preds$as))
  expect_equal(nrow(preds), 3)
})
