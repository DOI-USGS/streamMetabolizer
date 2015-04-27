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
  expect_equal(convert_date_to_doyhr(as.POSIXct("2020-01-01 00:00:00", tz="GMT")), 1, tol=0.000001, info="Jan 1 should be 1")
  expect_equal(convert_date_to_doyhr(as.POSIXct("2020-01-01 00:00:00", tz="CST6CDT")), 1, tol=0.000001, info="use the same timezone as the arg date")
  expect_equal(convert_date_to_doyhr(as.POSIXct("2020-01-01 00:00:00")), 1, tol=0.000001, info="should work in any tester's default tz")
  expect_equal(convert_date_to_doyhr(as.POSIXct("2020-01-01 01:00:00")), 1+1/24, tol=0.000001, info="decimal should include hours")
  expect_equal(convert_date_to_doyhr(as.POSIXct("2020-01-01 01:03:58")), 1+(60+3+58/60)/(24*60), tol=0.000001, info="decimal should include minutes and seconds")
  expect_equal(convert_date_to_doyhr(as.POSIXct("2004-12-01 00:00:00")), 1+convert_date_to_doyhr(as.POSIXct("2019-12-01 00:00:00")), tol=0.000001, info="should catch leap days")
  expect_equal(convert_date_to_doyhr(as.POSIXct("2016-05-29 01:00:00", tz="CST6CDT")), 150, info="treat numbers as true time since jan 1, ignoring daylight time")
  # doyhr to date
  expect_equal(convert_doyhr_to_date(1, year=1920), as.POSIXct("1920-01-01 00:00:00", tz="GMT"), info="1 should be Jan 1, default GMT")
  expect_equal(convert_doyhr_to_date(3, year=2016, tz="CST6CDT"), as.POSIXct("2016-01-03 00:00:00", tz="CST6CDT"), info="tz should stay as indicated")
  expect_equal(convert_doyhr_to_date(150, year=2016, tz="CST6CDT"), as.POSIXct("2016-05-29 01:00:00", tz="CST6CDT"), info="treat numbers as true time since jan 1, ignoring daylight time")
  # there and back
  expect_equal(convert_date_to_doyhr(convert_doyhr_to_date(12, year=2024, tz="CST6CDT")), 12, info="should preserve DOY in there&back")
  expect_equal(convert_date_to_doyhr(convert_doyhr_to_date(120, year=1998, tz="CST6CDT")), 120, info="should preserve DOY in there&back even during daylight savings")
  expect_equal(convert_doyhr_to_date(convert_date_to_doyhr(as.POSIXct("2007-06-15 12:00:00", tz="MST7MDT")), 2007, tz="MST7MDT"), 
               as.POSIXct("2007-06-15 12:00:00", tz="MST7MDT"), info="preserve date in there&back") 
})