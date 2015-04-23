
context("basic light model")

test_that("can generate light predictions from basic light model", {
  # radian-degree conversions
  expect_equal(streamMetabolizer:::rad(180), pi)
  expect_equal(streamMetabolizer:::deg(pi*3/2), 270)

  # declination angle
  decdf <- data.frame(
    jday=1:366, 
    dec=streamMetabolizer:::model_declination_angle(1:366, format="degrees"),
    dec_rad=streamMetabolizer:::model_declination_angle(1:366, format="radia"))
  expect_true(all(!is.na(decdf$dec)), "valid return for days between 1 and 366")
  equinox <- which(decdf$jday == as.numeric(format(as.Date("2015-03-20"), "%j")))
  solstice <- which(decdf$jday == as.numeric(format(as.Date("2015-06-21"), "%j")))
  expect_more_than(decdf[solstice,"dec"], decdf[equinox,"dec"], "solstice more declined than equinox")
  expect_equal(max(decdf$dec), -1*(min(decdf$dec)), info="equal & opposite declination on summer/winter solstices")
  expect_equal(rad(decdf$dec), decdf$dec_rad, info="degree & radian formats agree")
  
  # hour angle
  hourdf <- data.frame(
    hour=c(0:24), 
    hragl=streamMetabolizer:::model_hour_angle(0:24))
  expect_equal(hourdf[hourdf$hour==12, "hragl"], 0, info="hour angle is 0 at solar noon")
  expect_equal(-1*hourdf[0:13, "hragl"], hourdf[25:13, "hragl"], info="hour angles revers at noon")
  expect_equal(rad(hourdf$hragl), streamMetabolizer:::model_hour_angle(0:24, format="radians"))
  
  # zenith angle
  library(dplyr)
  zendf <- 
    data.frame(
      lat=rep(c(0,20,40,60, 80), each=365),
      jday=rep(1:365, times=5), 
      hour=12) %>%
    transform(
      dec=streamMetabolizer:::model_declination_angle(jday),
      hragl=streamMetabolizer:::model_hour_angle(hour)) %>%
    transform(
      zen=streamMetabolizer:::model_zenith_angle(lat, dec, hragl))
  expect_true(all(filter(zendf, lat==80)$zen > filter(zendf, lat==40)$zen), "higher latitude, higher zenith angle (lower sun)")
  
  # insolation
  insdf <- data.frame(
    lat=rep(c(0,20,40,60,80), each=24*4),
    jday=rep(rep(c(1,101,201,301), each=24), times=5), 
    hour=rep(c(0:12,13.5:23.5), times=4*5))
  insdf <- transform(insdf, ins=model_solar_insolation(jday, hour, lat))
  expect_true(all(insdf$ins >= 0), "non-negative insolation, always")
  expect_true(all(
    (insdf %>% group_by(lat, jday) %>% 
      summarize(daily_peak=max(ins)) %>% 
      summarize(summer_more_than_winter=all(daily_peak[jday==201] >= daily_peak[jday==1])))$summer_more_than_winter), 
    info="summertime noon insolation exceeds wintertime noon insolation at all latitudes")
  ins_u <- transform(insdf, ins=model_solar_insolation(jday, hour, lat, attach.units=TRUE))$ins
  expect_true(verify_units(ins_u, "J s^-1 m^-2", list(TRUE,FALSE)), "if requested, returns expected units")
  
})