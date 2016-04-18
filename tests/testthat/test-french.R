context("french creek data")

test_that("French Creek data are similar for streamMetabolizer & Bob Hall's code", {
  
  # load both datasets
  fx <- streamMetabolizer:::load_french_creek(attach.units=FALSE)
  fy <- streamMetabolizer:::load_french_creek_std()
  expect_equal(dim(fx), dim(fy))
  expect_equal(names(fx), names(fy))
  expect_equal(sapply(unitted::v(fx), class), sapply(fy, class))
  expect_equal(unitted::v(fx$DO.obs), fy$DO.obs)
  expect_equal(unitted::v(fx$solar.time), fy$solar.time)
  
  # combine for further comparison
  fxy <- dplyr::full_join(unitted::v(fx), fy, by="solar.time")
  expect_equal(nrow(fx), nrow(fxy))
  expect_equal(fxy$solar.time.x, fxy$solar.time.y)
  
  # check values that should be completely equal
  expect_equal(fxy$DO.obs.x, fxy$DO.obs.y)
  expect_equal(fxy$depth.x, fxy$depth.y)
  expect_equal(fxy$temp.water.x, fxy$temp.water.y)
  expect_equal(fxy$DO.sat.x, fxy$DO.sat.y)
  
  # check values that should be pretty much equal - light
  expect_gt(cor(fxy$light.x, fxy$light.y), 0.9999)
  expect_lt(abs(coef(lm(fxy$light.y ~ fxy$light.x))['(Intercept)']), 0.4)
  expect_lt(abs(coef(lm(fxy$light.y ~ fxy$light.x))['fxy$light.x'] - 1), 0.01)
  # library(ggplot2); ggplot(fxy, aes(x=light.x, y=light.y)) + geom_abline() + geom_point(alpha=0.5) + theme_bw()
  # library(ggplot2); ggplot(fxy, aes(x=solar.time)) + geom_line(aes(y=light.x), color='red') + geom_line(aes(y=light.y), color='blue') + theme_bw()
  # library(ggplot2); ggplot(fxy[7100:7350,], aes(x=solar.time)) + geom_line(aes(y=light.x), color='red') + geom_line(aes(y=light.y), color='blue') + theme_bw()
  
})


test_that("French Creek predictions are similar for streamMetabolizer & Bob Hall's code", {
  
  # set the date in several formats
  start.chron <- chron::chron(dates="08/23/12", times="22:00:00")
  end.chron <- chron::chron(dates="08/25/12", times="06:00:00")
  start.posix <- as.POSIXct(format(start.chron, "%Y-%m-%d %H:%M:%S"), tz="UTC")
  end.posix <- as.POSIXct(format(end.chron, "%Y-%m-%d %H:%M:%S"), tz="UTC")
  mid.date <- as.Date(start.posix + (end.posix - start.posix)/2)
  start.numeric <- as.numeric(start.posix - as.POSIXct(format(mid.date, "%Y-%m-%d 00:00:00"), tz="UTC"), units='hours')
  end.numeric <- as.numeric(end.posix - as.POSIXct(format(mid.date, "%Y-%m-%d 00:00:00"), tz="UTC"), units='hours')
  
  # get, format, & subset data
  vfrench <- streamMetabolizer:::load_french_creek(attach.units=FALSE)
  vfrenchshort <- vfrench[vfrench$solar.time >= start.posix & vfrench$solar.time <= end.posix, ]
  
  # dates & subsetting specific to nighttime regression
  first.dark <- 100 + which(vfrenchshort$light[101:nrow(vfrenchshort)] < 0.1)[1]
  stop.dark <- 100 + which(format(vfrenchshort$solar.time[101:nrow(vfrenchshort)], "%H:%M") == "23:00")[1]
  vfrenchnight <- vfrenchshort[first.dark:stop.dark,]
  night.start <- eval(parse(text=format(vfrenchnight$solar.time[1], "%H + %M/60")))
  night.end <- eval(parse(text=format(vfrenchnight$solar.time[nrow(vfrenchnight)], "%H + %M/60")))
  
  # PRK (metab_mle)
  smest <- get_fit(metab(
    specs=specs('m_np_oi_eu_km.nlm', day_start=start.numeric, day_end=end.numeric),
    data=vfrenchshort))[1,c("GPP","ER","K600","minimum")]
  bobest <- streamMetabolizer:::load_french_creek_std_mle(vfrenchshort, estimate='PRK')
  expect_lt(abs(smest$GPP - bobest$GPP), 0.01) #, info=paste0("GPP by SM: ", smest$GPP, "; by Bob: ", bobest$GPP))
  expect_lt(abs(smest$ER - bobest$ER), 0.01) #, info=paste0("ER by SM: ", smest$ER, "; by Bob: ", bobest$ER))
  expect_lt(abs(smest$K600 - bobest$K), 0.01) #, info=paste0("K600 by SM: ", smest$K600, "; by Bob: ", bobest$K))
  expect_lt(abs(smest$minimum - bobest$lik), 0.000001)
  
  # K (metab_night)
  smest <- predict_metab(metab(
    specs=specs(mm_name('night'), day_start=night.start, day_end=night.end, day_tests=c('full_day', 'even_timesteps', 'complete_data')),
    data=vfrenchnight))[c("GPP","ER","K600")]
  bobest <- streamMetabolizer:::load_french_creek_std_mle(
    vfrenchnight, estimate='K', start=chron::chron(dates="08/24/12", times="18:45:00"), end=chron::chron(dates="08/24/12", times="22:56:00"))
  expect_lt(abs(smest$K600 - bobest$K), 0.0001) #, info=paste0("K600 by SM: ", smest$K600, "; by Bob: ", bobest$K))
  
  # PR (metab_mle)
  smest <- get_fit(metab_mle(
    specs=specs('m_np_oi_eu_km.nlm', day_start=start.numeric, day_end=end.numeric),
    data=vfrenchshort, data_daily=data.frame(date=mid.date, K600=35)))[,c("GPP","ER","K600","minimum")]
  bobest <- streamMetabolizer:::load_french_creek_std_mle(vfrenchshort, estimate='PR', K=35)
  expect_lt(abs(smest$GPP - bobest$GPP), 0.02) #, info=paste0("GPP by SM: ", smest$GPP, "; by Bob: ", bobest$GPP))
  expect_lt(abs(smest$ER - bobest$ER), 0.01) #, info=paste0("ER by SM: ", smest$ER, "; by Bob: ", bobest$ER))
  expect_lt(abs(smest$minimum - bobest$lik), 0.000001)
  
  # Bayes w/ Bob's MLE-PRK for comparison - really loose criteria for prediction agreement
  prkest <- get_fit(metab_mle(
    specs=specs('m_np_oi_eu_km.nlm', day_start=start.numeric, day_end=end.numeric),
    data=vfrenchshort))[,c("GPP","ER","K600","minimum")]
  bobest <- streamMetabolizer:::load_french_creek_std_mle(vfrenchshort, estimate='PRK')
  mb <- metab_bayes(
    specs=specs('b_np_oi_eu_km.jags', saved_steps=4000, day_start=start.numeric, day_end=end.numeric), 
    data=vfrenchshort)
  smest <- predict_metab(mb)[,c("GPP","ER","K600")]
  expect_lt(abs(smest$GPP - bobest$GPP), 0.2) #, info=paste0("GPP by SM: ", smest$GPP, "; by Bob: ", bobest$GPP))
  expect_lt(abs(smest$ER - bobest$ER), 0.2) #, info=paste0("ER by SM: ", smest$ER, "; by Bob: ", bobest$ER))
  expect_lt(abs(smest$K600 - bobest$K), 3) #, info=paste0("K600 by SM: ", smest$K600, "; by Bob: ", bobest$K))
  
})
