
test_that("plot_DO_preds runs with style='dygraphs' without error", {
  skip_if_not_installed("dygraphs")
  skip_if_not_installed("xts")

  mm <- metab_night(specs(mm_name('night')), data=data_metab('3', day_start=12, day_end=36))
  DO_preds <- predict_DO(mm)
  result <- plot_DO_preds(DO_preds, style='dygraphs', y_var='conc')

  expect_s3_class(result, "dygraphs")
})

test_that("dygraphs gap-row logic inserts one NA sentinel row per date group", {
  # Verifies the group_modify() replacement for do() in plot_DO_preds.R.
  # The transformation appends a NA-padded copy of each group's last row to
  # create a visual break between dates in dygraphs output.
  mm <- metab_night(specs(mm_name('night')), data=data_metab('3', day_start=12, day_end=36))
  DO_preds <- predict_DO(mm)

  # Build the conc subset of preds_xts_input — same structure as inside plot_DO_preds()
  input <- unitted::v(DO_preds) %>%
    mutate(as='conc', pure=NA, mod=DO.mod, obs=DO.obs) %>%
    arrange(solar.time)

  n_dates <- length(unique(input$date))  # 3 with data_metab('3')

  result <- input %>%
    group_by(date) %>%
    group_modify(~ {
      out <- .x[c(seq_len(nrow(.x)), nrow(.x)), ]
      out[nrow(.x)+1, c('pure','mod','obs')] <- NA
      out
    }) %>%
    ungroup()

  # Exactly one sentinel appended per date group
  expect_equal(nrow(result), nrow(input) + n_dates)

  # All key columns survive the transformation
  expect_true(all(c('date', 'solar.time', 'pure', 'mod', 'obs', 'as') %in% names(result)))

  # The last row of each date group is the NA sentinel
  sentinels <- result %>% group_by(date) %>% slice_tail(n=1) %>% ungroup()
  expect_equal(nrow(sentinels), n_dates)
  expect_true(all(is.na(sentinels$mod)))
  expect_true(all(is.na(sentinels$obs)))
})
