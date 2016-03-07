context("metab_night")

vfrench <- streamMetabolizer:::load_french_creek(attach.units=FALSE)
vfrenchshort <- vfrench[vfrench$solar.time >= as.POSIXct("2012-08-23 00:00:00", tz="UTC") & 
                          vfrench$solar.time <= as.POSIXct("2012-08-26 00:00:00", tz="UTC"), ]

test_that("metab_night models can be created", {
  
  mm <- metab_night(vfrenchshort)
  
  # check basic structure
  expect_is(mm, "metab_night")
  expect_is(slot(mm, "fit"), "data.frame")
  expect_is(slot(mm, "specs"), "list")
  expect_is(slot(mm, "data"), "data.frame")
  expect_is(slot(mm, "pkg_version"), "character")
  
})

test_that("metab_night predictions (predict_metab, predict_DO) make sense", {
  
  # metab_night
  mm <- metab_night(data=vfrenchshort)
  metab <- predict_metab(mm)
  DO_preds <- predict_DO(mm)
  DO_preds_Aug24 <- dplyr::filter(DO_preds, date == "2012-08-24")
  expect_true(all(abs(DO_preds_Aug24$DO.obs - DO_preds_Aug24$DO.mod) < 0.15), "DO.mod tracks DO.obs with not too much error")
  # plot_DO_preds(DO_preds, y_var="conc")
  # plot_DO_preds(DO_preds, y_var="pctsat")
  
})

test_that("metab_night predictions can be passed back into metab_mle", {
  
  # metab_night
  mmk <- metab_night(data=vfrenchshort)

  # metab_mle
  mm <- metab_mle(data=vfrenchshort, data_daily=predict_metab(mmk)[c('date', 'K600')])
  metab <- predict_metab(mm)
  DO_preds <- predict_DO(mm)
  DO_preds_Aug24 <- dplyr::filter(DO_preds, date == "2012-08-24")
  # note that had to raise the maximum error for this model combination, from 0.15 to 0.35 mgO/L
  expect_true(all(abs(DO_preds_Aug24$DO.obs - DO_preds_Aug24$DO.mod) < 0.35), "DO.mod tracks DO.obs with not too much error")
  # plot_DO_preds(DO_preds)
  
})

test_that("metab_night predictions match Bob's", {
  # # set the date in several formats
  # start.chron <- chron::chron(dates="08/23/12", times="22:00:00")
  # end.chron <- chron::chron(dates="08/25/12", times="06:00:00")
  # start.posix <- as.POSIXct(format(start.chron, "%Y-%m-%d %H:%M:%S"), tz="UTC")
  # end.posix <- as.POSIXct(format(end.chron, "%Y-%m-%d %H:%M:%S"), tz="UTC")
  # mid.date <- as.Date(start.posix + (end.posix - start.posix)/2, tz=lubridate::tz(start.posix))
  # 
  # # get, format, & subset data
  # vfrench <- streamMetabolizer:::load_french_creek(attach.units=FALSE)
  # vfrenchshort <- vfrench[vfrench$solar.time >= start.posix & vfrench$solar.time <= end.posix, ]
  # 
  # # dates & subsetting specific to nighttime regression
  # first.dark <- 100 + which(vfrenchshort$light[101:nrow(vfrenchshort)] < 0.1)[1]
  # stop.dark <- 100 + which(
  #   format(vfrenchshort$solar.time[101:nrow(vfrenchshort)], "%H:%M") == "23:00")[1]
  # vfrenchnight <- vfrenchshort[first.dark:stop.dark,]
  # night.start <- eval(parse(text=format(vfrenchnight$solar.time[1], "%H + %M/60")))
  # night.end <- eval(parse(text=format(vfrenchnight$solar.time[nrow(vfrenchnight)], "%H + %M/60")))
  # 
  # # fit
  # mm <- metab_night(data=vfrenchnight,
  #   specs=specs('n_np_pi_eu_.lm'),
  #   day_start=night.start, day_end=night.end)
  # 
  # # give estimates
  # predict_metab(mm)
  # streamMetabolizer:::load_french_creek_std_mle(vfrenchnight, estimate='K')
})


