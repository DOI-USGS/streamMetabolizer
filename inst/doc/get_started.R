## ----knitr_init, echo=FALSE, cache=FALSE----------------------------------------------------------
knitr::opts_chunk$set(echo = TRUE)
options(width=100)

## ---- eval=FALSE----------------------------------------------------------------------------------
#  install.packages("streamMetabolizer", dependencies=TRUE,
#    repos=c("https://owi.usgs.gov/R","https://cran.rstudio.com"))

## ---- eval=FALSE----------------------------------------------------------------------------------
#  update.packages(oldPkgs=c("streamMetabolizer","unitted"),
#    dependencies=TRUE, repos=c("https://owi.usgs.gov/R", "https://cran.rstudio.com"))

## ---- eval=FALSE----------------------------------------------------------------------------------
#  install.packages(
#    c("LakeMetabolizer","unitted","dplyr","lazyeval","lubridate","magrittr",
#      "tidyr","chron","dygraphs","ggplot2","RCurl","rstan","XML","xts"),
#    repos=c("https://owi.usgs.gov/R","https://cran.rstudio.com"))

## ---- eval=FALSE----------------------------------------------------------------------------------
#  devtools::install_github("USGS-R/streamMetabolizer", ref="develop")

## ----libs, warning=FALSE, message=FALSE-----------------------------------------------------------
library(dplyr)
library(tidyr)
library(ggplot2)
library(streamMetabolizer)

## ----data-----------------------------------------------------------------------------------------
dat <- data_metab(num_days='3', res='15', day_start=4, day_end=28, attach.units=TRUE)

## ----data_check-----------------------------------------------------------------------------------
dim(dat)
dat[c(1,48,96,240,288),] # some example rows

## ----viz_inputs_DO, fig.width=7, fig.height=3-----------------------------------------------------
dat %>% unitted::v() %>%
  mutate(DO.pctsat = 100 * (DO.obs / DO.sat)) %>%
  select(solar.time, starts_with('DO')) %>%
  gather(type, DO.value, starts_with('DO')) %>%
  mutate(units=ifelse(type == 'DO.pctsat', 'DO\n(% sat)', 'DO\n(mg/L)')) %>%
  ggplot(aes(x=solar.time, y=DO.value, color=type)) + geom_line() + 
  facet_grid(units ~ ., scale='free_y') + theme_bw() +
  scale_color_discrete('variable')

## ----viz_inputs_other, fig.width=7, fig.height=4--------------------------------------------------
labels <- c(depth='depth\n(m)', temp.water='water temp\n(deg C)', light='PAR\n(umol m^-2 s^-1)')
dat %>% unitted::v() %>%
  select(solar.time, depth, temp.water, light) %>%
  gather(type, value, depth, temp.water, light) %>%
  mutate(
    type=ordered(type, levels=c('depth','temp.water','light')),
    units=ordered(labels[type], unname(labels))) %>%
  ggplot(aes(x=solar.time, y=value, color=type)) + geom_line() + 
  facet_grid(units ~ ., scale='free_y') + theme_bw() +
  scale_color_discrete('variable')

## ----data_needs-----------------------------------------------------------------------------------
metab_inputs('mle', 'data')

## ----data_str-------------------------------------------------------------------------------------
data.frame(
  colname=names(dat), 
  class=unname(sapply(unitted::v(dat), function(col) paste(class(col), collapse=','))),
  units=unname(unitted::get_units(dat)))

## ----names----------------------------------------------------------------------------------------
three_names <- c(
  mm_name(type='mle'), # the default MLE model
  mm_name(type='mle', ode_method='euler'), # override the default ode_method
  mm_name(type='bayes')) # the default Bayesian model
three_names

# parse the above model names
mm_parse_name(three_names)

## ----mle_name-------------------------------------------------------------------------------------
mle_name <- mm_name(type='mle')
mle_name

## ----mle_specs------------------------------------------------------------------------------------
mle_specs <- specs(mle_name)
mle_specs

## ----specs_details--------------------------------------------------------------------------------
mle_specs <- specs(mle_name, init.GPP.daily=2, init.ER.daily=-1, init.K600.daily=3)

## ----mle_fit, warning=FALSE-----------------------------------------------------------------------
mle_fit <- metab(mle_specs, data=dat, info=c(site='French Creek, WY', source='Bob Hall'))

## ----show-----------------------------------------------------------------------------------------
mle_fit

## ----info-----------------------------------------------------------------------------------------
get_info(mle_fit)
head(get_data(mle_fit))

## ----info2----------------------------------------------------------------------------------------
get_fitting_time(mle_fit) # the time it took to fit the model
get_version(mle_fit) # the streamMetabolizer version used to fit the model
get_specs(mle_fit) # the specifications we passed in

## ----plot_metab1, fig.width=7, fig.height=4.5-----------------------------------------------------
plot_metab_preds(mle_fit)

## ----plot_metab2, fig.width=7, fig.height=6-------------------------------------------------------
plot_DO_preds(mle_fit)

## ----pred_dfs-------------------------------------------------------------------------------------
met_preds <- predict_metab(mle_fit)
met_preds
oxy_preds <- predict_DO(mle_fit)
head(oxy_preds)

## ----bayes_name-----------------------------------------------------------------------------------
bayes_name <- mm_name(
  type='bayes', pool_K600='none', 
  err_obs_iid=TRUE, err_proc_acor=FALSE, err_proc_iid=TRUE, 
  ode_method='trapezoid')
bayes_name

## ----bayes_specs----------------------------------------------------------------------------------
bayes_specs <- specs(bayes_name)
bayes_specs

## ----bayes_specs2---------------------------------------------------------------------------------
# one way to alter specifications: call specs() again
bayes_specs <- specs(bayes_name, burnin_steps=100, saved_steps=200, n_cores=1, GPP_daily_mu=3, GPP_daily_sigma=2)
# another way: use revise()
bayes_specs <- revise(bayes_specs, burnin_steps=100, saved_steps=200, n_cores=1, GPP_daily_mu=3, GPP_daily_sigma=2)

## ----bayes_fit------------------------------------------------------------------------------------
bayes_fit <- metab(bayes_specs, data=dat)

## -------------------------------------------------------------------------------------------------
bayes_fit

## ----bayes_warning--------------------------------------------------------------------------------
# here's where you'd find fitting messages:
select(get_params(bayes_fit), warnings, errors)
get_fit(bayes_fit)$warnings
get_fit(bayes_fit)$errors
# and prediction messages
select(predict_metab(bayes_fit), warnings, errors)

## ----bayes_pred_tbl, results='asis'---------------------------------------------------------------
predict_metab(bayes_fit) %>% 
  lapply(function(col) if(is.numeric(col)) round(col, 2) else col ) %>%
  as.data.frame() %>%
  knitr::kable()

## ---- fig.width=7, fig.height=6-------------------------------------------------------------------
plot_DO_preds(bayes_fit)

## -------------------------------------------------------------------------------------------------
valid_names <- mm_valid_names(type=c('bayes','mle','night'))
length(valid_names)
c(valid_names[seq(1,length(valid_names),length.out=20)], '...')

