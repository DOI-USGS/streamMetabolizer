## ----knitr_init, echo=FALSE, cache=FALSE----------------------------------------------------------
knitr::opts_chunk$set(echo = TRUE)
options(width=100)

## ---- message=FALSE-------------------------------------------------------------------------------
library(streamMetabolizer)
library(dplyr)
library(ggplot2)

## -------------------------------------------------------------------------------------------------
dat <- data_metab('3','30')

## -------------------------------------------------------------------------------------------------
mm_euler <- metab(specs(mm_name('mle', ode_method='euler')), dat)

## -------------------------------------------------------------------------------------------------
mm_trapezoid <- metab(specs(mm_name('mle', ode_method='trapezoid')), dat)

## -------------------------------------------------------------------------------------------------
mm_rk4 <- metab(specs(mm_name('mle', ode_method='rk4')), dat) 

## -------------------------------------------------------------------------------------------------
mm_lsoda <- metab(specs(mm_name('mle', ode_method='lsoda')), dat) 

## -------------------------------------------------------------------------------------------------
DO.standard <- rep(predict_DO(mm_rk4)$'DO.mod', times=4)
ode_preds <- bind_rows(
  mutate(predict_DO(mm_euler), method='euler'),
  mutate(predict_DO(mm_trapezoid), method='trapezoid'),
  mutate(predict_DO(mm_rk4), method='rk4'),
  mutate(predict_DO(mm_lsoda), method='lsoda')) %>%
  mutate(DO.mod.diffeuler = DO.mod - DO.standard)

## ---- fig.width=6---------------------------------------------------------------------------------
ggplot(ode_preds, aes(x=solar.time)) +
  geom_point(aes(y=DO.obs), color='grey', alpha=0.3) +
  geom_line(aes(y=DO.mod, color=method), size=1) +
  theme_bw()

## ---- fig.width=6---------------------------------------------------------------------------------
ggplot(ode_preds, aes(x=solar.time)) +
  geom_point(aes(y=pmax(-0.2, pmin(0.2, DO.mod.diffeuler)), color=method), size=1, alpha=0.8) +
  scale_y_continuous(limits=c(-0.2,0.2)) +
  ylab("Deviations from rk4 (capped at +/- 0.2)") +
  theme_bw()

