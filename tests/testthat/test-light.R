
context("basic light model")
test_that("can generate light predictions from basic light model", {
  # radian-degree conversions
  expect_equal(streamMetabolizer:::to_radians(180), pi)
  expect_equal(streamMetabolizer:::to_degrees(pi*3/2), 270)

  # declination angle
  decdf <- data.frame(
    jday=1:366, 
    dec=streamMetabolizer:::calc_declination_angle(1:366, format="degrees"),
    dec_rad=streamMetabolizer:::calc_declination_angle(1:366, format="radia"))
  expect_true(all(!is.na(decdf$dec)), "valid return for days between 1 and 366")
  equinox <- which(decdf$jday == as.numeric(format(as.Date("2015-03-20"), "%j")))
  solstice <- which(decdf$jday == as.numeric(format(as.Date("2015-06-21"), "%j")))
  expect_more_than(decdf[solstice,"dec"], decdf[equinox,"dec"], "solstice more declined than equinox")
  expect_equal(max(decdf$dec), -1*(min(decdf$dec)), info="equal & opposite declination on summer/winter solstices")
  expect_equal(streamMetabolizer:::to_radians(decdf$dec), decdf$dec_rad, info="degree & radian formats agree")
  
  # hour angle
  hourdf <- data.frame(
    hour=c(0:24), 
    hragl=streamMetabolizer:::calc_hour_angle(0:24))
  expect_equal(hourdf[hourdf$hour==12, "hragl"], 0, info="hour angle is 0 at solar noon")
  expect_equal(-1*hourdf[0:13, "hragl"], hourdf[25:13, "hragl"], info="hour angles reverse at noon")
  expect_equal(streamMetabolizer:::to_radians(hourdf$hragl), streamMetabolizer:::calc_hour_angle(0:24, format="radians"))
  
  # zenith angle
  library(dplyr)
  zendf <- 
    data.frame(
      lat=rep(c(0,20,40,60, 80), each=365),
      jday=rep(1:365, times=5), 
      hour=12) %>%
    transform(
      dec=streamMetabolizer:::calc_declination_angle(jday),
      hragl=streamMetabolizer:::calc_hour_angle(hour)) %>%
    transform(
      zen=streamMetabolizer:::calc_zenith_angle(lat, dec, hragl))
  expect_true(all(zendf[zendf$lat==80, "zen"] > zendf[zendf$lat==40, "zen"]), "higher latitude, higher zenith angle (lower sun)")
  
  # insolation
  insdf <- data.frame(
    lat=rep(c(0,20,40,60,80), each=24*4),
    jday=rep(rep(c(1,101,201,301), each=24), times=5), 
    hour=rep(c(0:12,13.5:23.5), times=4*5))
  insdf <- transform(insdf, ins=calc_solar_insolation(jday, hour, lat))
  expect_true(all(insdf$ins >= 0), "non-negative insolation, always")
  expect_true(all(
    (insdf %>% group_by(lat, jday) %>% 
      summarize(daily_peak=max(ins)) %>% 
      summarize(summer_more_than_winter=all(daily_peak[jday==201] >= daily_peak[jday==1])))$summer_more_than_winter), 
    info="summertime noon insolation exceeds wintertime noon insolation at all latitudes")
  ins_u <- transform(insdf, ins=calc_solar_insolation(jday, hour, lat, attach.units=TRUE))$ins
  expect_true(verify_units(ins_u, "J s^-1 m^-2", list(TRUE,FALSE)), "if requested, returns expected units")
  
})


test_that("calc_solar_insolation has consistent output with that of calc_sun_rise_set and calc_is_daytime", {
  library(dplyr)
  insdf <- data.frame(
    #   lat=rep(c(0,20,40,60,80), each=366),
    #   jday=rep(1:366, times=5)) %>%
    lat=rep(c(0,20,40,60,80), each=18),
    jday=rep(seq(1, 360, by=20), times=5)) %>%
    transform(
      date.time=strptime(sprintf("2000-%d 00",jday), format="%Y-%j %H")) %>% 
    transform(
      lm_sunrise=calc_sun_rise_set(date.time, lat)[,1],
      lm_sunset=calc_sun_rise_set(date.time, lat)[,2]) %>%
    group_by(jday, lat) %>%
    do(with(., {
      # compare to streamMetabolizer method, which determines light at any given time
      hours <- seq(0,23.95,by=0.05)
      insol <- calc_solar_insolation(jday, hour=hours, lat=lat)
      whichdaytime <- which(insol > 0.00001)
      if(any(!is.na(whichdaytime)))
        sm_daytime <- date.time + hours[range(whichdaytime)]*60*60
      else 
        sm_daytime <- c(NA,NA)
      
      # compare to calc_is_daytime (id) method (from LakeMetabolizer), which determines whether it is light at any given time
      isday <- calc_is_daytime(date.time+hours*60*60, lat=lat)
      whichdaytime <- which(isday)
      if(any(!is.na(whichdaytime)))
        id_daytime <- date.time + hours[range(whichdaytime)]*60*60
      else 
        id_daytime <- c(NA,NA)
      
      # put together
      data.frame(., sm_sunrise=sm_daytime[1], sm_sunset=sm_daytime[2], id_sunrise=id_daytime[1], id_sunset=id_daytime[2])
    }))
  
  # Visual inspection
  # library(tidyr)
  # instidy <- suppressWarnings(insdf %>% 
  #   gather(pkg, sunrise, lm_sunrise, sm_sunrise) %>%
  #   gather(pkg, sunset, lm_sunset, sm_sunset)) %>%
  #   transform(
  #     sunrise=as.numeric(as.POSIXct(sunrise, origin="1970-01-01") - trunc(as.POSIXct(sunrise, origin="1970-01-01"), "day"), units="hours"), 
  #     sunset=as.numeric(as.POSIXct(sunset, origin="1970-01-01") - trunc(as.POSIXct(sunset, origin="1970-01-01"), "day"), units="hours"))
  # library(ggplot2)
  # ggplot(instidy, aes(x=date.time, y=sunrise, color=pkg, group=lat)) + geom_point() + theme_bw()
  
  # expect_that inspection
  diffs <- insdf %>%
    transform(
      diff_sunrise = as.numeric(sm_sunrise - lm_sunrise, units="mins"),
      diff_sunset = as.numeric(sm_sunset - lm_sunset, units="mins")) %>%
    group_by(lat) %>%
    summarize(
      max_diff_sunrise = max(abs(diff_sunrise), na.rm=TRUE),
      max_diff_sunset = max(abs(diff_sunset), na.rm=TRUE))
  expect_less_than(max(filter(diffs, lat==0)[,c("max_diff_sunrise","max_diff_sunset")]), 3.5)
  expect_less_than(max(filter(diffs, lat==20)[,c("max_diff_sunrise","max_diff_sunset")]), 5)
  expect_less_than(max(filter(diffs, lat==40)[,c("max_diff_sunrise","max_diff_sunset")]), 7.5)
  expect_less_than(max(filter(diffs, lat==60)[,c("max_diff_sunrise","max_diff_sunset")]), 12)
  expect_less_than(max(filter(diffs, lat==80)[,c("max_diff_sunrise","max_diff_sunset")]), 100) # difference of up to 97 minutes at high lat when run for all days
})