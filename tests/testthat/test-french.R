context("french creek data")

test_that("French Creek data are similar for streamMetabolizer & Bob Hall's code", {
  
  # load both datasets
  fx <- streamMetabolizer:::load_french_creek(attach.units=FALSE)
  fy <- streamMetabolizer:::load_french_creek_std()
  expect_equal(dim(fx), dim(fy))
  expect_equal(names(fx), names(fy))
  expect_equal(sapply(unitted::v(fx), class), sapply(fy, class))
  expect_equal(unitted::v(fx$DO.obs), fy$DO.obs)
  expect_equal(unitted::v(fx$local.time), fy$local.time)
  
  # combine for further comparison
  fxy <- dplyr::full_join(unitted::v(fx), fy, by="local.time")
  expect_equal(nrow(fx), nrow(fxy))
  expect_equal(fxy$local.time.x, fxy$local.time.y)
  
  # check values that should be completely equal
  expect_equal(fxy$DO.obs.x, fxy$DO.obs.y)
  expect_equal(fxy$depth.x, fxy$depth.y)
  expect_equal(fxy$temp.water.x, fxy$temp.water.y)
  expect_equal(fxy$DO.sat.x, fxy$DO.sat.y)
  
  # check values that should be pretty much equal - light
  expect_more_than(cor(fxy$light.x, fxy$light.y), 0.9999)
  expect_less_than(abs(coef(lm(fxy$light.y ~ fxy$light.x))['(Intercept)']), 0.4)
  expect_less_than(abs(coef(lm(fxy$light.y ~ fxy$light.x))['fxy$light.x'] - 1), 0.01)
  # library(ggplot2); ggplot(fxy, aes(x=light.x, y=light.y)) + geom_abline() + geom_point(alpha=0.5) + theme_bw()
  # library(ggplot2); ggplot(fxy, aes(x=local.time)) + geom_line(aes(y=light.x), color='red') + geom_line(aes(y=light.y), color='blue') + theme_bw()
  # library(ggplot2); ggplot(fxy[7100:7350,], aes(x=local.time)) + geom_line(aes(y=light.x), color='red') + geom_line(aes(y=light.y), color='blue') + theme_bw()
  
})


test_that("French Creek predictions are similar for streamMetabolizer & Bob Hall's code", {
  
  # set the date in several formats
  start.chron <- chron::chron(dates="08/23/12", times="22:00:00")
  end.chron <- chron::chron(dates="08/25/12", times="06:00:00")
  start.posix <- as.POSIXct(format(start.chron, "%Y-%m-%d %H:%M:%S"), tz="Etc/GMT+7")
  end.posix <- as.POSIXct(format(end.chron, "%Y-%m-%d %H:%M:%S"), tz="Etc/GMT+7")
  mid.date <- as.Date(start.posix + (end.posix - start.posix)/2)
  start.numeric <- as.numeric(start.posix - as.POSIXct(format(mid.date, "%Y-%m-%d 00:00:00"), tz="Etc/GMT+7"), units='hours')
  end.numeric <- as.numeric(end.posix - as.POSIXct(format(mid.date, "%Y-%m-%d 00:00:00"), tz="Etc/GMT+7"), units='hours')
  
  # get, format, & subset data
  vfrench <- streamMetabolizer:::load_french_creek(attach.units=FALSE)
  vfrenchshort <- vfrench[vfrench$local.time >= start.posix & vfrench$local.time <= end.posix, ]
  
  # dates & subsetting specific to nighttime regression
  first.dark <- 100 + which(vfrenchshort$light[101:nrow(vfrenchshort)] < 0.1)[1]
  stop.dark <- 100 + which(format(vfrenchshort$local.time[101:nrow(vfrenchshort)], "%H:%M") == "23:00")[1]
  vfrenchnight <- vfrenchshort[first.dark:stop.dark,]
  night.start <- eval(parse(text=format(vfrenchnight$local.time[1], "%H + %M/60")))
  night.end <- eval(parse(text=format(vfrenchnight$local.time[nrow(vfrenchnight)], "%H + %M/60")))
  
  # PRK (metab_mle)
  smest <- get_fit(metab_mle(data=vfrenchshort, day_start=start.numeric, day_end=end.numeric, 
                             model_specs=specs_mle_nopool_oi(ODE_method="Euler")))[2,c("GPP","ER","K600","minimum")]
  bobest <- streamMetabolizer:::load_french_creek_std_mle(vfrenchshort, estimate='PRK')
  expect_less_than(abs(smest$GPP - bobest$GPP), 0.01, info=paste0("GPP by SM: ", smest$GPP, "; by Bob: ", bobest$GPP))
  expect_less_than(abs(smest$ER - bobest$ER), 0.01, info=paste0("ER by SM: ", smest$ER, "; by Bob: ", bobest$ER))
  expect_less_than(abs(smest$K600 - bobest$K), 0.01, info=paste0("K600 by SM: ", smest$K600, "; by Bob: ", bobest$K))
  expect_less_than(abs(smest$minimum - bobest$lik), 0.000001)
  
  # K (metab_night)
  smest <- predict_metab(metab_night(data=vfrenchnight, day_start=night.start, day_end=night.end))[c("GPP","ER","K600")]
  bobest <- streamMetabolizer:::load_french_creek_std_mle(vfrenchnight, estimate='K')
  expect_less_than(abs(smest$K600 - bobest$K), 0.0001, info=paste0("K600 by SM: ", smest$K600, "; by Bob: ", bobest$K))
  
  # PR (metab_mle)
  smest <- get_fit(metab_mle(data=vfrenchshort, data_daily=data.frame(local.date=mid.date, K600=35), day_start=start.numeric, day_end=end.numeric, 
                             model_specs=specs_mle_nopool_oi(ODE_method="Euler")))[2,c("GPP","ER","K600","minimum")]
  bobest <- streamMetabolizer:::load_french_creek_std_mle(vfrenchshort, estimate='PR', K=35)
  expect_less_than(abs(smest$GPP - bobest$GPP), 0.02, info=paste0("GPP by SM: ", smest$GPP, "; by Bob: ", bobest$GPP))
  expect_less_than(abs(smest$ER - bobest$ER), 0.01, info=paste0("ER by SM: ", smest$ER, "; by Bob: ", bobest$ER))
  expect_less_than(abs(smest$minimum - bobest$lik), 0.000001)
  
  # Bayes w/ Bob's MLE-PRK for comparison
  prkest <- get_fit(metab_mle(data=vfrenchshort, day_start=start.numeric, day_end=end.numeric, 
                              model_specs=specs_mle_nopool_oi(ODE_method="Euler")))[2,c("GPP","ER","K600","minimum")]
  bobest <- streamMetabolizer:::load_french_creek_std_mle(vfrenchshort, estimate='PRK')
  mb <- metab_bayes(data=vfrenchshort, model_specs=specs_bayes_jags_nopool_oi(saved_steps = 4000, model_file="nopool_oi_Euler.jags"),
                    day_start=start.numeric, day_end=end.numeric)
  smest <- predict_metab(mb)[2,c("GPP","ER","K600")]
  expect_less_than(abs(smest$GPP - bobest$GPP), 0.05, info=paste0("GPP by SM: ", smest$GPP, "; by Bob: ", bobest$GPP))
  expect_less_than(abs(smest$ER - bobest$ER), 0.05, info=paste0("ER by SM: ", smest$ER, "; by Bob: ", bobest$ER))
  expect_less_than(abs(smest$K600 - bobest$K), 0.5, info=paste0("K600 by SM: ", smest$K600, "; by Bob: ", bobest$K))
  
})
