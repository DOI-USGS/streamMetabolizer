## ----knitr_init, echo=FALSE, cache=FALSE---------------------------------
library(knitr)
library(rmdformats)

## Global options
options(max.print="75")
opts_chunk$set(echo=TRUE,
	             prompt=FALSE,
               tidy=TRUE)
opts_knit$set(width=75)

## ---- messages=TRUE, warnings=TRUE, errors=TRUE--------------------------
library(streamMetabolizer)

## ------------------------------------------------------------------------
suppressPackageStartupMessages(library(dplyr))

## ------------------------------------------------------------------------
dat <- data_metab('3', '15')

## ------------------------------------------------------------------------
# the Classic: linear GPP, constant ER (also the default)
mm_classic <- 
  mm_name('mle', GPP_fun='linlight', ER_fun='constant') %>% 
  specs() %>%
  metab(dat)
mm_classic

## ------------------------------------------------------------------------
# the Saturator: GPP saturating with light, constant ER
mm_saturator <- 
  mm_name('mle', GPP_fun='satlight', ER_fun='constant') %>% 
  specs() %>%
  metab(dat)
mm_saturator

## ------------------------------------------------------------------------
get_params(mm_saturator) %>% select(date, warnings, errors)

## ------------------------------------------------------------------------
predict_metab(mm_saturator) %>% select(date, warnings, errors)

## ------------------------------------------------------------------------
predict_DO(mm_saturator) %>% head
plot_DO_preds(mm_saturator)

## ------------------------------------------------------------------------
mm_saturator2 <- 
  mm_name('mle', GPP_fun='satlight', ER_fun='constant') %>% 
  specs() %>%
  metab(dat, data_daily=select(get_params(mm_saturator), date, init.Pmax=Pmax, init.alpha=alpha))
get_params(mm_saturator2)
mm_saturator3 <- 
  mm_name('mle', GPP_fun='satlight', ER_fun='constant') %>% 
  specs(init.Pmax=6.2, init.alpha=0.008) %>%
  metab(dat)
get_params(mm_saturator3)
mm_saturator4 <- 
  mm_name('mle', GPP_fun='satlight', ER_fun='constant') %>% 
  specs(init.Pmax=6.2, init.alpha=0.008) %>%
  metab(dat, transmute(get_params(mm_saturator), date, init.Pmax=Pmax[1], init.alpha=alpha[1])[2,])
get_params(mm_saturator4)

## ------------------------------------------------------------------------
plot_DO_preds(mm_classic)
plot_DO_preds(mm_saturator4)

