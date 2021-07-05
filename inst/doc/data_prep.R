## ----setup, include=FALSE-----------------------------------------------------
knitr::opts_chunk$set(echo = TRUE)

## ----data, message=FALSE------------------------------------------------------
library(streamMetabolizer)
dat <- data_metab(num_days='3', res='15')

## ----data_check---------------------------------------------------------------
dim(dat)
dat[c(1,48,96,240,288),] # some example rows

## ---- message=FALSE-----------------------------------------------------------
library(dplyr)
library(tidyr)
library(ggplot2)

## ----viz_inputs_DO, fig.width=7, fig.height=3---------------------------------
dat %>% 
  mutate(DO.pctsat = 100 * (DO.obs / DO.sat)) %>%
  select(solar.time, starts_with('DO')) %>%
  gather(type, DO.value, starts_with('DO')) %>%
  mutate(units=ifelse(type == 'DO.pctsat', 'DO\n(% sat)', 'DO\n(mg/L)')) %>%
  ggplot(aes(x=solar.time, y=DO.value, color=type)) + geom_line() + 
  facet_grid(units ~ ., scale='free_y') + theme_bw() +
  scale_color_discrete('variable')

## ----viz_inputs_other, fig.width=7, fig.height=4------------------------------
labels <- c(depth='depth\n(m)', temp.water='water temp\n(deg C)', light='PAR\n(umol m^-2 s^-1)')
dat %>% 
  select(solar.time, depth, temp.water, light) %>%
  gather(type, value, depth, temp.water, light) %>%
  mutate(
    type=ordered(type, levels=c('depth','temp.water','light')),
    units=ordered(labels[type], unname(labels))) %>%
  ggplot(aes(x=solar.time, y=value, color=type)) + geom_line() + 
  facet_grid(units ~ ., scale='free_y') + theme_bw() +
  scale_color_discrete('variable')

## ----data_needs---------------------------------------------------------------
metab_inputs('mle', 'data')

## -----------------------------------------------------------------------------
num.time <- 1471867200
(posix.time.localtz <- as.POSIXct(num.time, origin='1970-01-01', tz='UTC'))

## -----------------------------------------------------------------------------
num.time <- 1471867200
(posix.time.nominalUTC <- as.POSIXct(num.time, origin='1970-01-01', tz='UTC')) # the numbers get treated as UTC no matter what tz you request
(posix.time.localtz <- lubridate::force_tz(posix.time.nominalUTC, 'Etc/GMT+7')) # +7 = mountain standard time

## -----------------------------------------------------------------------------
text.time <- '2016-08-22 12:00:00'
(posix.time.utc <- as.POSIXct(text.time, tz='UTC'))

## -----------------------------------------------------------------------------
text.time <- '2016-08-22 12:00:00'
(posix.time.localtz <- as.POSIXct(text.time, format="%Y-%m-%d %H:%M:%S", tz='America/New_York'))

## -----------------------------------------------------------------------------
text.time <- '2016-08-22 12:00:00'
(posix.time.localtz <- as.POSIXct(text.time, format="%Y-%m-%d %H:%M:%S", tz='Etc/GMT+5'))

## -----------------------------------------------------------------------------
chron.time <- chron::chron('08/22/16', '12:00:00')
time.format <- "%Y-%m-%d %H:%M:%S"
text.time <- format(chron.time, time.format) # direct as.POSIXct time works poorly
(posix.time.localtz <- as.POSIXct(text.time, format=time.format, tz='Etc/GMT+5'))

## -----------------------------------------------------------------------------
lubridate::tz(posix.time.localtz) # yep, we want and have the code for EST
(posix.time.solar <- streamMetabolizer::calc_solar_time(posix.time.localtz, longitude=-106.3))

