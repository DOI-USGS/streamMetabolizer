## ----setup, include=FALSE------------------------------------------------
knitr::opts_chunk$set(echo = TRUE)

## ------------------------------------------------------------------------
num.time <- 1471867200
(posix.time.localtz <- as.POSIXct(num.time, origin='1970-01-01', tz='UTC'))

## ------------------------------------------------------------------------
num.time <- 1471867200
(posix.time.nominalUTC <- as.POSIXct(num.time, origin='1970-01-01', tz='UTC')) # the numbers get treated as UTC no matter what tz you request
(posix.time.localtz <- lubridate::force_tz(posix.time.nominalUTC, 'Etc/GMT+7')) # +7 = mountain standard time

## ------------------------------------------------------------------------
text.time <- '2016-08-22 12:00:00'
(posix.time.utc <- as.POSIXct(text.time, tz='UTC'))

## ------------------------------------------------------------------------
text.time <- '2016-08-22 12:00:00'
(posix.time.localtz <- as.POSIXct(text.time, format="%Y-%m-%d %H:%M:%S", tz='America/New_York'))

## ------------------------------------------------------------------------
text.time <- '2016-08-22 12:00:00'
(posix.time.localtz <- as.POSIXct(text.time, format="%Y-%m-%d %H:%M:%S", tz='Etc/GMT+5'))

## ------------------------------------------------------------------------
chron.time <- chron::chron('08/22/16', '12:00:00')
time.format <- "%Y-%m-%d %H:%M:%S"
text.time <- format(chron.time, time.format) # direct as.POSIXct time works poorly
(posix.time.localtz <- as.POSIXct(text.time, format=time.format, tz='Etc/GMT+5'))

## ------------------------------------------------------------------------
lubridate::tz(posix.time.localtz) # yep, we want and have the code for EST
(posix.time.solar <- streamMetabolizer::calc_solar_time(posix.time.localtz, longitude=-106.3))

