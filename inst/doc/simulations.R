## ----knitr_init, echo=FALSE, cache=FALSE----------------------------------------------------------
knitr::opts_chunk$set(echo = TRUE)
options(width=100)

## ---- messages=TRUE, warnings=TRUE, errors=TRUE---------------------------------------------------
library(streamMetabolizer)

## -------------------------------------------------------------------------------------------------
suppressPackageStartupMessages(library(dplyr))

## -------------------------------------------------------------------------------------------------
dat <- data_metab('3', '15')

## -------------------------------------------------------------------------------------------------
name_sim_q10 <- mm_name('sim', ER_fun='q10temp')

## ---- error=TRUE----------------------------------------------------------------------------------
mm_sim_q10_trial <- metab(specs(name_sim_q10), dat)
get_params(mm_sim_q10_trial)

## -------------------------------------------------------------------------------------------------
params_sim_q10a <- data.frame(date=as.Date(paste0('2012-09-',18:20)), GPP.daily=2.1, ER20=-5:-3, K600.daily=16)
params_sim_q10a

## -------------------------------------------------------------------------------------------------
mm_mle_q10 <- metab(specs(mm_name('mle', ER_fun='q10temp')), data=dat)

## -------------------------------------------------------------------------------------------------
params_sim_q10b <- get_params(mm_mle_q10, uncertainty='none', messages=FALSE)
params_sim_q10b

## ---- eval=FALSE----------------------------------------------------------------------------------
#  ?mm_name

## -------------------------------------------------------------------------------------------------
specs_sim_q10 <- specs(name_sim_q10, err.obs.sigma=0, err.proc.sigma=1)
specs_sim_q10

## -------------------------------------------------------------------------------------------------
mm_sim_q10a <- metab(specs_sim_q10, dat, data_daily=params_sim_q10a)
mm_sim_q10b <- metab(specs_sim_q10, dat, data_daily=params_sim_q10b)

## ---- fig.height=2--------------------------------------------------------------------------------
plot_DO_preds(mm_sim_q10a, y_var='conc')
plot_DO_preds(mm_sim_q10b, y_var='conc')

## -------------------------------------------------------------------------------------------------
specs_sim_sat <- specs(mm_name('sim', GPP_fun='satlight'), err.obs.sigma=0, err.proc.sigma=1)
params_sim_sat <- get_params(metab(specs(mm_name('mle', GPP_fun='satlight')), data=dat), uncertainty='none', messages=FALSE)

## ---- fig.height=2--------------------------------------------------------------------------------
mm_sim_sat_i <- metab(specs_sim_sat, dat, data_daily=params_sim_sat)
plot_DO_preds(mm_sim_sat_i, y_var='conc')
plot_DO_preds(mm_sim_sat_i, y_var='conc')

## ---- fig.height=2--------------------------------------------------------------------------------
mm_sim_sat_f <- metab(revise(specs_sim_sat, sim.seed=47), dat, data_daily=params_sim_sat)
plot_DO_preds(mm_sim_sat_f, y_var='conc')
plot_DO_preds(mm_sim_sat_f, y_var='conc')

## -------------------------------------------------------------------------------------------------
dat <- data_metab('10', '30')
params <- data.frame(date=as.Date(paste0('2012-09-',18:27)), Pmax=rnorm(10, 6, 2), alpha=rnorm(10, 0.01, 0.001), ER20=rnorm(10, -4, 2), K600.daily=16)
specs <- specs(mm_name('sim', GPP_fun='satlight', ER_fun='q10temp'), err.obs.sigma=0.2, err.proc.sigma=1)
mm <- metab(specs, data=dat, data_daily=params)

## -------------------------------------------------------------------------------------------------
mm

## -------------------------------------------------------------------------------------------------
plot_metab_preds(mm)

