## ----knitr_init, echo=FALSE, cache=FALSE--------------------------------------
knitr::opts_chunk$set(echo = TRUE)
options(width=80)

## ----libs, warning=FALSE, message=FALSE---------------------------------------
library(streamMetabolizer)
library(dplyr)

## ----data---------------------------------------------------------------------
dat <- data_metab(num_days='3', res='15', day_start=4, day_end=28)

## ----bayes_name---------------------------------------------------------------
bayes_name <- mm_name(type='bayes', pool_K600='none', err_obs_iid=TRUE, err_proc_iid=TRUE)
bayes_name

## ----bayes_specs--------------------------------------------------------------
bayes_specs <- specs(bayes_name)
bayes_specs

## ----bayes_specs2-------------------------------------------------------------
# one way to alter specifications: call specs() again
bayes_specs <- specs(bayes_name, burnin_steps=100, saved_steps=200, n_cores=1, GPP_daily_mu=3, GPP_daily_sigma=2)
# another way: use revise()
bayes_specs <- revise(bayes_specs, burnin_steps=100, saved_steps=200, n_cores=1, GPP_daily_mu=3, GPP_daily_sigma=2)

## ----bayes_fit----------------------------------------------------------------
mm <- metab(bayes_specs, data=dat)

## -----------------------------------------------------------------------------
mm

## ----bayes_pred_tbl-----------------------------------------------------------
predict_metab(mm)

## ----bayes_pred_fig, fig.width=5, fig.height=5--------------------------------
plot_metab_preds(mm)

## -----------------------------------------------------------------------------
get_params(mm)

## ----bayes_pdo_tbl, results='asis'--------------------------------------------
predict_DO(mm) %>% head()

## ----bayes_pdo_fig, fig.width=5, fig.height=5---------------------------------
plot_DO_preds(mm)

## ---- fig.width=5, fig.height=5-----------------------------------------------
mcmc <- get_mcmc(mm)
rstan::traceplot(mcmc, pars='K600_daily', nrow=3)

## -----------------------------------------------------------------------------
get_fit(mm)$overall %>%
  select(ends_with('Rhat'))

## -----------------------------------------------------------------------------
get_fit(mm) %>%
  lapply(names)

