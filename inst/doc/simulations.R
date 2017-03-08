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
specs_sim_q10 <- specs(name_sim_q10, err_obs_sigma=0, err_proc_sigma=1, K600_daily=NULL, GPP_daily=NULL, ER20=NULL)
specs_sim_q10

## -------------------------------------------------------------------------------------------------
mm_sim_q10a <- metab(specs_sim_q10, dat, data_daily=params_sim_q10a)
mm_sim_q10b <- metab(specs_sim_q10, dat, data_daily=params_sim_q10b)

## ---- fig.height=2--------------------------------------------------------------------------------
head(predict_DO(mm_sim_q10a))
head(predict_DO(mm_sim_q10b))
plot_DO_preds(mm_sim_q10a, y_var='conc')
plot_DO_preds(mm_sim_q10b, y_var='conc')

## -------------------------------------------------------------------------------------------------
specs_sim_sat <- specs(mm_name('sim', GPP_fun='satlight'), err_obs_sigma=0, err_proc_sigma=1, K600_daily=NULL, Pmax=NULL, alpha=NULL, ER_daily=NULL)
params_sim_sat <- get_params(metab(specs(mm_name('mle', GPP_fun='satlight')), data=dat), uncertainty='none', messages=FALSE)

## ---- fig.height=2--------------------------------------------------------------------------------
mm_sim_sat_i <- metab(specs_sim_sat, dat, data_daily=params_sim_sat)
plot_DO_preds(mm_sim_sat_i, y_var='conc')
plot_DO_preds(mm_sim_sat_i, y_var='conc')

## ---- fig.height=2--------------------------------------------------------------------------------
mm_sim_sat_f <- metab(revise(specs_sim_sat, sim_seed=47), dat, data_daily=params_sim_sat)
plot_DO_preds(mm_sim_sat_f, y_var='conc')
plot_DO_preds(mm_sim_sat_f, y_var='conc')

## -------------------------------------------------------------------------------------------------
dat <- data_metab('10', '30')
params <- data.frame(date=as.Date(paste0('2012-09-',18:27)), Pmax=rnorm(10, 6, 2), alpha=rnorm(10, 0.01, 0.001), ER20=rnorm(10, -4, 2), K600.daily=16)
specs <- specs(mm_name('sim', GPP_fun='satlight', ER_fun='q10temp'), err_obs_sigma=0.2, err_proc_sigma=1, K600_daily=NULL, Pmax=NULL, alpha=NULL, ER20=NULL)
mm <- metab(specs, data=dat, data_daily=params)

## -------------------------------------------------------------------------------------------------
mm

## -------------------------------------------------------------------------------------------------
plot_metab_preds(mm)

## -------------------------------------------------------------------------------------------------
dat <- data_metab('10','15')
datlen <- as.numeric(diff(range(dat$solar.time)) + as.difftime(15, units='mins'), units='days')
dat20 <- bind_rows(lapply((0:1)*10, function(add) {
  mutate(dat, solar.time = solar.time + as.difftime(add, units='days'))
}))

## -------------------------------------------------------------------------------------------------
sp <- specs(mm_name('sim'))
lapply(unclass(sp)[c('K600_daily','GPP_daily','ER_daily')], function(fun) {
  list(code=attr(fun, 'srcref'), example_vals=fun(n=10))
})

## -------------------------------------------------------------------------------------------------
mm <- metab(sp, dat20, data_daily=NULL)
get_params(mm)[c('date','K600.daily','GPP.daily','ER.daily')]

## -------------------------------------------------------------------------------------------------
sp <- specs('sim', err_obs_sigma=function(n, ...) -0.01*((1:n) - (n/2))^2 + 1, err_proc_sigma=function(n, ...) rnorm(n, 0.1, 0.005), err_proc_phi=seq(0, 1, length.out=20), GPP_daily=3, ER_daily=-4, K600_daily=16)
mm <- metab(sp, dat20)
get_params(mm)
plot_DO_preds(mm)

## -------------------------------------------------------------------------------------------------
sp <- specs('sim', err_obs_sigma=0.01, err_proc_sigma=function(n, ...) rep(c(0.5, 4), each=10), err_proc_phi=rep(seq(0, 0.8, length.out=10), times=2), GPP_daily=3, ER_daily=-4, K600_daily=16)
mm <- metab(sp, dat20)
get_params(mm)
plot_DO_preds(mm)

## -------------------------------------------------------------------------------------------------
sp <- specs('sim', err_obs_sigma=0.01, err_proc_sigma=0.4, K600_daily=16, GPP_daily=function(n, ...) round(rnorm(n, 4, 1), 1), ER_daily=function(GPP.daily, ...) GPP.daily*-2)
mm <- metab(sp, dat20)
get_params(mm)

## -------------------------------------------------------------------------------------------------
sp <- specs(mm_name('sim', pool_K600='binned', ER_fun='q10temp'), sim_seed=6332)

## -------------------------------------------------------------------------------------------------
mm <- metab(revise(sp, K600_daily=function(n, K600_daily_predlog, ...) pmax(0, rnorm(n, exp(K600_daily_predlog), 0.4))), dat20)
pars <- get_params(mm)
pars

## -------------------------------------------------------------------------------------------------
attr(pars, 'K600_eqn')

## -------------------------------------------------------------------------------------------------
library(ggplot2)
KQ <- as.data.frame(attr(pars, 'K600_eqn')[c('K600_lnQ_nodes_centers', 'lnK600_lnQ_nodes')])
Kpred <- mutate(select(pars, date, discharge.daily, K600.daily), K600_daily_predlog=attr(pars, 'K600_eqn')$K600_daily_predlog)
ggplot(KQ, aes(x=K600_lnQ_nodes_centers, y=lnK600_lnQ_nodes)) + geom_line(color='blue') + geom_point(color='blue') + 
  geom_point(data=Kpred, aes(x=log(discharge.daily), y=K600_daily_predlog), color='purple') +
  geom_point(data=Kpred, aes(x=log(discharge.daily), y=log(K600.daily)), color='red')

