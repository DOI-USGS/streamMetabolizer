context("conversion functions")

test_that("converting between k600 and kgas works", {
  # k600 to kgas
  expect_equal(convert_k600_to_kGAS(k600=20, temperature=15), 18.698, tol=0.001)
  expect_equal(convert_k600_to_kGAS(k600=20, temperature=15, gas="O2"), 18.698, tol=0.001)
  expect_equal(convert_k600_to_kGAS(k600=20, temperature=15, gas="CO2"), 17.361, tol=0.001)
  # kgas to k600
  expect_equal(convert_kGAS_to_k600(kGAS=18.698, temperature=15), 20, tol=0.001)
  # there and back
  expect_equal(convert_k600_to_kGAS(convert_kGAS_to_k600(kGAS=15, temperature=15), temperature=15), 15)
  expect_equal(convert_kGAS_to_k600(convert_k600_to_kGAS(k600=10, temperature=15), temperature=15), 10)
})

test_that("converting between SW and PAR works", {
  # sw to par
  expect_equal(convert_SW_to_PAR(sw=800), 1691.2)
  expect_equal(convert_SW_to_PAR(sw=800, coeff=2), 1600)
  # par to sw
  expect_equal(convert_PAR_to_SW(par=400), 189.2)
  expect_equal(convert_PAR_to_SW(par=400, coeff=0.5), 200)
  # there and back
  expect_equal(convert_PAR_to_SW(convert_SW_to_PAR(sw=800)), 800, tol=0.0001)
  expect_equal(convert_SW_to_PAR(convert_PAR_to_SW(par=800)), 800, tol=0.0001)
})

test_that("converting between date and DOY works", {
  # date to doyhr
  expect_equal(convert_date_to_doyhr(as.POSIXct("2020-01-01 00:00:00", tz="UTC")), 1, tol=0.000001, info="Jan 1 should be 0")
  expect_equal(convert_date_to_doyhr(as.POSIXct("2020-01-01 00:00:00", tz="CST6CDT")), 1, tol=0.000001, info="use the same timezone as the arg date")
  expect_equal(convert_date_to_doyhr(as.POSIXct("2020-01-01 00:00:00")), 1, tol=0.000001, info="should work in any tester's default tz")
  expect_equal(convert_date_to_doyhr(as.POSIXct("2020-01-01 01:00:00")), 1+1/24, tol=0.000001, info="decimal should include hours")
  expect_equal(convert_date_to_doyhr(as.POSIXct("2020-01-01 01:03:58")), 1+(60+3+58/60)/(24*60), tol=0.000001, info="decimal should include minutes and seconds")
  expect_equal(convert_date_to_doyhr(as.POSIXct("2004-12-01 00:00:00")), 1+convert_date_to_doyhr(as.POSIXct("2019-12-01 00:00:00")), tol=0.000001, info="should catch leap days")
  expect_equal(convert_date_to_doyhr(as.POSIXct("2016-05-29 01:00:00", tz="America/Chicago")), 150, info="treat numbers as true time since jan 1, ignoring daylight time")
  # doyhr to date
  expect_equal(convert_doyhr_to_date(1, year=1920), as.POSIXct("1920-01-01 00:00:00", tz="UTC"), info="1 should be Jan 1, default UTC")
  expect_equal(convert_doyhr_to_date(3, year=2016, tz="CST6CDT"), as.POSIXct("2016-01-03 00:00:00", tz="CST6CDT"), info="tz should stay as indicated")
  expect_equal(convert_doyhr_to_date(150, year=2016, tz="CST6CDT"), as.POSIXct("2016-05-29 01:00:00", tz="CST6CDT"), info="treat numbers as true time since jan 1, ignoring daylight time")
  # there and back
  expect_equal(convert_date_to_doyhr(convert_doyhr_to_date(12, year=2024, tz="CST6CDT")), 12, info="should preserve DOY in there&back")
  expect_equal(convert_date_to_doyhr(convert_doyhr_to_date(120, year=1998, tz="America/Chicago")), 120, info="should preserve DOY in there&back even during daylight savings")
  expect_equal(convert_doyhr_to_date(convert_date_to_doyhr(as.POSIXct("2007-05-29 01:00:00", tz="America/Denver")), 2007, tz="America/Denver"), 
               as.POSIXct("2007-05-29 01:00:00", tz="America/Denver"), info="preserve date in there&back") 
})

test_that("converting between UTC and solar time works", {
  library(unitted)
  adate <- as.POSIXct("2014-04-01 00:00:00", tz="UTC")
  somedates <- seq(adate, adate+as.difftime(365*2, units="days"), by=as.difftime(10.35, units="days"))
  # UTC to solar
  expect_equal(convert_UTC_to_solartime(adate, longitude=u(0, "degW"), time.type="mean solar"), adate)
  expect_equal(convert_UTC_to_solartime(adate, longitude=u(0, "degW"), time.type="apparent solar"), adate+as.difftime(-4.661701, units="mins"), tol=0.0001)
  expect_error(convert_UTC_to_solartime(adate, longitude=u(0, "degW"), time.type="not a type"), "should be one of", info="only accept valid time.types")
  expect_lt(as.numeric(convert_UTC_to_solartime(adate, longitude=u(105.3, "degE"), time.type="mean solar") - (adate + as.difftime(7, units="hours")), units='secs'), 10) #go east, be later
  expect_lt(as.numeric(convert_UTC_to_solartime(adate, longitude=u(105.3, "degW"), time.type="mean solar") - (adate - as.difftime(7, units="hours")), units='secs'), 10) #go west, be earlier
  expect_equal(convert_UTC_to_solartime(adate, longitude=u(89, "degW"), time.type="apparent solar"), 
               convert_UTC_to_solartime(adate, longitude=-89, time.type="apparent solar"), info="negative degrees are degW")
  expect_equal(convert_UTC_to_solartime(u(adate), longitude=71, time.type="apparent solar"), 
               convert_UTC_to_solartime(adate, longitude=u(-71, "degW"), time.type="apparent solar"), info="mix&match units is OK for this fun")
  expect_equal(convert_UTC_to_solartime(somedates, longitude=u(0, "degW"), time.type="mean solar"), somedates, info="handle multiple dates")
  # solar to UTC
  expect_equal(convert_solartime_to_UTC(adate, longitude=u(0, "degW"), time.type="mean solar"), adate)
  expect_equal(convert_solartime_to_UTC(adate, longitude=u(0, "degW"), time.type="apparent solar"), adate+as.difftime(+4.661701, units="mins"), tol=0.0001)
  expect_error(convert_solartime_to_UTC(adate, longitude=u(0, "degW"), time.type="not a type"), "should be one of", info="only accept valid time.types")
  expect_lt(as.numeric(convert_solartime_to_UTC(adate, longitude=u(105.3, "degE"), time.type="mean solar") - (adate - as.difftime(7, units="hours")), units='secs'), 10) #, info="go east, be later")
  expect_lt(as.numeric(convert_solartime_to_UTC(adate, longitude=u(105.3, "degW"), time.type="mean solar") - (adate + as.difftime(7, units="hours")), units='secs'), 10) #, info="go west, be earlier")
  expect_equal(convert_solartime_to_UTC(somedates, longitude=u(0, "degW"), time.type="mean solar"), somedates, info="handle multiple dates")
  expect_equal(as.numeric(convert_solartime_to_UTC(somedates, longitude=u(0, "degW"), time.type="apparent solar")), as.numeric(somedates), tol=1000, info="handle multiple dates")
  # there and back
  expect_equal(convert_UTC_to_solartime(convert_solartime_to_UTC(adate, longitude=-103.8, time.type="mean solar"), longitude=-103.8, time.type="mean solar"), adate)
  expect_equal(convert_solartime_to_UTC(convert_UTC_to_solartime(adate, longitude=-103.8, time.type="mean solar"), longitude=-103.8, time.type="mean solar"), adate)
  expect_equal(as.numeric(convert_UTC_to_solartime(convert_solartime_to_UTC(adate, longitude=-103.8, time.type="app"), longitude=-103.8, time.type="apparent solar")), as.numeric(adate), tol=6)
  expect_equal(as.numeric(convert_solartime_to_UTC(convert_UTC_to_solartime(adate, longitude=-103.8, time.type="app"), longitude=-103.8, time.type="appar")), as.numeric(adate), tol=6)
  expect_equal(convert_UTC_to_solartime(convert_solartime_to_UTC(somedates, longitude=-103.8, time.type="mean solar"), longitude=-103.8, time.type="mean solar"), somedates)
  expect_equal(convert_solartime_to_UTC(convert_UTC_to_solartime(somedates, longitude=-103.8, time.type="mean solar"), longitude=-103.8, time.type="mean solar"), somedates)
  expect_equal(as.numeric(convert_UTC_to_solartime(convert_solartime_to_UTC(somedates, longitude=-103.8, time.type="app"), longitude=-103.8, time.type="apparent solar")), as.numeric(somedates), tol=10)
  expect_equal(as.numeric(convert_solartime_to_UTC(convert_UTC_to_solartime(somedates, longitude=-103.8, time.type="app"), longitude=-103.8, time.type="appar")), as.numeric(somedates), tol=10)  
})


test_that("converting between UTC and local time works", {
  library(unitted)
  adate <- as.POSIXct("2014-02-01 00:00:00", tz="UTC")
  asummerdate <- as.POSIXct("2014-07-04 12:14:16", tz="UTC")
  somedates <- seq(adate, adate+as.difftime(365*2, units="days"), by=as.difftime(10.35, units="days"))
  # UTC to local
  #   what happens in london stays in london:
  expect_equal(convert_UTC_to_localtime(adate, latitude=u(51.48, "degN"), longitude=u(0, "degW"), time.type="standard local"), adate)
  expect_equal(as.numeric(convert_UTC_to_localtime(adate, latitude=u(51.48, "degN"), longitude=u(0, "degW"), time.type="daylight local")), as.numeric(adate), info="still the same #s")
  expect_equal(lubridate::tz(convert_UTC_to_localtime(adate, latitude=u(51.48, "degN"), longitude=u(0, "degW"), time.type="daylight local")), "Europe/London", info="different tz name")
  #   error checking
  expect_error(convert_UTC_to_localtime(adate, latitude=u(51.48, "degN"), longitude=0, time.type="standard"), "unitted")
  expect_error(convert_UTC_to_localtime(adate, latitude=u(51.48, "degN"), longitude=u(0, "degW"), time.type="not a type"), 'should be one of "standard local", "daylight local"', info="only accept valid time.types")
  #   real time changes
  # "POSIX has positive signs west of Greenwich" - http://opensource.apple.com/source/system_cmds/system_cmds-230/zic.tproj/datfiles/etcetera
  expect_equal(lubridate::tz(convert_UTC_to_localtime(adate, latitude=u(40, "degN"), longitude=u(105.3, "degE"), time.type="standard")), "Etc/GMT-8", info="go east, be POSIX-negative")
  expect_equal(lubridate::tz(convert_UTC_to_localtime(adate, latitude=u(37, "degN"), longitude=u(105.3, "degW"), time.type="standard")), "Etc/GMT+7", info="go west, be POSIX-positive")
  expect_equal(lubridate::tz(convert_UTC_to_localtime(adate, latitude=u(37, "degN"), longitude=u(-105.3, "degE"), time.type="standard")), "Etc/GMT+7", info="go west, be POSIX-positive")
  expect_equal(convert_UTC_to_localtime(somedates, latitude=u(34, "degN"), longitude=u(80, "degW"), time.type="daylight"), lubridate::with_tz(somedates, "America/New_York"), info="handle multiple dates")
  expect_equal(convert_UTC_to_localtime(somedates, latitude=u(34, "degN"), longitude=u(80, "degW"), time.type="standard"), lubridate::with_tz(somedates, "Etc/GMT+5"), info="handle multiple dates")

  # local to UTC
  #   std/daylight, winter/summer
  expect_equal(convert_localtime_to_UTC(adate), adate)
  expect_equal(convert_localtime_to_UTC(asummerdate), asummerdate)
  expect_equal(convert_localtime_to_UTC(lubridate::with_tz(adate, "Etc/GMT+8")), adate)
  expect_equal(convert_localtime_to_UTC(lubridate::with_tz(adate, "Etc/GMT-5")), adate)
  
  # there and back
  expect_equal(convert_UTC_to_localtime(convert_localtime_to_UTC(lubridate::with_tz(adate, "America/Denver")), latitude=40, longitude=-105.3, time.type="daylight"), lubridate::with_tz(adate, "America/Denver"))
  expect_equal(convert_localtime_to_UTC(convert_UTC_to_localtime(adate, latitude=40, longitude=-105.3, time.type="daylight")), adate)
  # not sure why only this next line would fail on Travis-CI, but it does. It works on my machine.
  expect_equal(convert_localtime_to_UTC(convert_UTC_to_localtime(adate, latitude=40, longitude=-103.8, time.type="standard")), adate)
  
})
