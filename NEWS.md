# 0.12.1

* updated README.md, code.json, and other package metadata files

# 0.12.0

* fixed compatibility issues with R 4.1.0 and recent (OK, not super recent)
changes to dplyr

* beginning to deprecate the use of `attach.units` arguments

# 0.11.4

* fixed issue with `data_metab` that arises with newer versions of `tidyverse`

# 0.11.3

* new options to compute and report GPP_pseudo_R2 and DO_R2 for Bayesian models

* smarter default for err_mult_GPP_sdlog_sigma (changed from 0.05 to 1, which
seems to be better for more streams)

# 0.11.2

* changed `err_proc_dayiid` to `err_proc_GPP` and switched to a more GPP-oriented
process error for those models

# 0.11.1

* added `err_proc_dayiid` option for Bayesian models

# 0.11.0

* added `GPP_fun='satlight'` option for Bayesian models

# 0.10.10

* added a test to mm_is_valid_day that checks for depths <= 0, which would break
or seriously confuse a model

* fixed a bug where a single day of data passed to a b_Kl or b_Kb model would
cause Stan to fail because the log(discharge) values were being formatted
improperly for Stan.

# 0.10.7

* non-standard parameters from custom Bayesian models are now sorted in a less
bewildering order in the model fit (as from `get_fit()`)

# 0.10.4

* `plot_DO_preds` with dygraphs works again, and now you can subset dates using
args `date_start` and `date_end` directly in `plot_DO_preds`

# 0.10.3

* adding a tiny bit of support toward bringing in new light-saturation models in
the future

# 0.10.2

* test for day lengths < 0 in `mm_model_by_ply`, thanks to @weisoon

* bugfix for `get_params` for sim models with fixed vector of Q and/or K values

# 0.10.1

* upgrade `unitted` dependency to v0.2.8 to accommodate recent change in 
function exporting requirements in R 3.3

* increase robustness of `mm_is_valid_day` for dates with 1 observation

# 0.10.0

* more comprehensive assignment of parameters from Stan output into model fit 
data.frames in Bayesian models

* `daily` and `inst` data.frames in Bayesian model fit now get date/timestamps

* structural changes in Bayesian models: fitted `DO_mod[1]` for state space 
models, probability constraints on `DO_mod[1]` and `DO_mod_partial[1]`

* structural change in Bayesian models: reindexing `err_obs_iid` and 
`err_proc_iid` to match other `inst` variables

# 0.9.48

* Extend filtering with `required_timestep` to metab_bayes

# 0.9.47

* Bug fix in use of `required_timestep` in `mm_is_valid_day`

# 0.9.46

* Removed the warning discouraging setting params_out within specs()

# 0.9.45

* Bayesian specs defaults now reflect a little more of the literature and our 
experience modeling metabolism

# 0.9.44

* New specs element `required_timestep` allows you to require that each date has
the specified numeric timestep in days

# 0.9.42-0.9.43

* Bugfixes for new test that excludes days with non-positive discharge

# 0.9.40

* Update to roxygen2 6.0.1

# 0.9.39

* Bayesian, MLE, and nighttime regression models can now all test for and 
exclude days with non-positive discharge

# 0.9.36

* Models should now be able to accept `tbl_df`s (dplyr/tibble format) for the 
`data` and `data_daily` arguments to `metab()`

# 0.9.35

* Bayesian models now distinguish between compilation time and fitting time

* Updates to `plot_distribs` for recent changes to Bayesian models

# 0.9.34

* Bayesian models with `pool_K600 != 'none'` can how have their 
`K600_daily_sigma` (or `K600_daily_sdlog`) be a fitted value, a fixed value, or 
a value fixed at 0

# 0.9.33

* new function: `calc_light_merged`, which merges modeled and observed light 
into a smooth curve

# 0.9.29

* functions in the `specs` of `sim` models can now refer to their own metabolism
model (and therefore also its `info` or `data_daily` slots)

# 0.9.28

* switched from `rlnorm` to `rnorm` for distribution of `K600_daily` around 
`K600_daily_pred` in linear and binned models

# 0.9.27

* `sim` models can now generate daily parameters from functions and can even 
generate binned K~Q relationships with random variation

# 0.9.25

* various improvements to flexibility and speed of `sim` models

# 0.9.22

* renamed `metab_sim` specs for consistency with other model arguments, e.g., 
`err.obs.sigma` is now `err_obs_sigma`

# 0.9.21

* simplified Bayesian K pooling models to make the pooling more effective (and 
also faster). Hierarchical bayesian models now fix rather than fit 
`K600_daily_sdlog`, the standard deviation of `K600_daily` relative to 
`K600_daily_predlog`

# 0.9.18

* solidified the model feature options, including error types, DO deficit 
source, ODE integration methods, and hierarchy

* revised hierarchical methods to be consistent with recent email exchanges.

# 0.9.15

* new function: `plot_distribs` to explore priors

* bayesian models can now track `err_obs_iid` and other parameters fit at the 
resolution of `data`

* new functions: `calc_solar_time` and `calc_light` for faster data preparation

* refined treatment of depth in trapezoid method - effect will mainly be 
noticeable for sites/times with rapidly changing depth

# 0.9.14

* all tests are passing and all examples are running [again]. tests are 
acceptably comprehensive [for now].

* known issue: K binning in bayesian models is producing unlikely estimates. all
bayesian models need further testing.

# 0.9.13

* ongoing bug fixes and re-integration of all model types into new ODE 
production/integration scheme

# 0.9.12

* removed JAGS from the package. Stan is better for us.

# 0.9.11

* this is the last version with JAGS in it.

# 0.9.10

* myriad bug fixes to accommodate changes from version 0.9.9.

# 0.9.9

* `metab_mle` and `metab_sim` now have several options for relationships among 
GPP, light, ER, and temperature.

* `metab_mle` and `metab_sim` now have several options for the method by which 
the differential equation for dDO/dt is numerically integrated to produce a time
series of DO predictions.

* `metab_mle` has taken a performance hit to become more flexible in the 
GPP-light relationship, the ODE method, and so on. Optimization is likely in the
future.

* `predict_metab` and `predict_DO` optionally attach units to their output.

# 0.9.8

* better error and warning handling in metab_bayes, including a new function 
`get_log()` that retrieves log file[s] from MCMC model compilation & run[s]. log
files are now retrieved for both JAGS and Stan models.

* updates to accommodate changes in dependency packages (tibble and dplyr)

* more efficient specifications of JAGS and Stan models

* incorporated feedback on vignette

# 0.9.6

* model names (from `mm_name()`) and bayesian model file names (in models 
folder) now include info on the GPP and ER functions - default is still pl = GPP
is a linear function of light, and rc = ER is constant over every 24-hour period

* more informative error messages for timesteps <= 0 in `mm_model_by_ply`

# 0.9.5.1

* in `metab_Kmodel`, now avoiding negative weights

# 0.9.5

* Bug fixes and error prevention

# 0.9.4

* Now automatically checks for available updates when you attach the package

* Improved units handling in `convert_k600`/`kGAS`

# 0.9.2

* Hierarchical constraints on K600 are now available! Options are 'normal', 
'linear', and 'binned'; see the details section on `pool_K600` in `?mm_name` and
the description of parameters starting with `K600_daily` in `?specs`.

* Interface change: `specs` lists now print more prettily and have class 'specs'
(though they're still fundamentally just lists)

* Vignette: see `vignette('getstarted')`

# 0.9.0

* New function: metab() serves as a gateway to all model types. You can now pass
specs to metab() and expect the appropriate model to be chosen and called based 
on model_name in the specs list.

* New function: data_metab() produces a dataset for testing/demonstration, with 
options for the resolution & flaws to introduce.

* Newly public function: mm_model_by_ply is now public. Its interface has also 
changed somewhat: tests has been renamed to day_tests, and validity tests are 
conducted within mm_model_by_ply if day_tests is not empty, and validity and 
timestep information are now passed to model_fun.

* Changed functionality: mm_model_by_ply_prototype() now produces a 1-row 
data.frame as well as a message, which helps this function demonstrate the 
workings of mm_model_by_ply(). mm_model_by_ply_prototype() is a lightweight 
example of a function that can be passed to mm_model_by_ply(), and its help file
describes the parameters such a function should expect.

* New function: mm_get_timestep() computes the mean and/or unique timestep[s] 
and optionally requires that there be just one unique timestep within a vector 
of times or dates.

* Interface change: the argument tests is now called day_tests in the metab(), 
metab_night(), etc., mm_model_by_ply(), and mm_is_valid_day().

* Interface change: day_start, day_end, and tests are now containined within 
specs rather than defined separately in the call to metab, metab_bayes, etc.

* Interface change: in metab(), metab_mle(), etc., the model_specs argument is 
now called specs.

* Interface change: metab functions now accept specs first, then data, 
data_daily, and info. (specs was renamed from model_specs; see above.) This 
permits chaining from mm_name to specs to metab.

* Interface change: get_args is now get_specs, and the result is a list of specs
as named in specs() rather than a list with an element called model_specs that 
is itself a list.

* Hierarchical bayesian models are now possible and include hierarchical 
parameters for distributions on error and K600 (normal, linear, and binned). 
Some models are known to work; complete testing for all models is forthcoming.

# 0.8.0

* Major interface change (renamed variable) to clarify types of time: solar.time
(mean solar time), app.solar.time (apparent solar time), local.time (time in 
local time zone). Metabolism models now accept solar.time rather than 
local.time, though it's still possible to pass in local time but just call it 
solar.time (as long as you don't have daylight savings time).

# 0.7.3

* Remove calc\_schmidt because it is never used

# 0.7.2

This package is not ready for use by many, but it does currently have:

* support for a wide range of non-hierarchical models, both Bayesian and 
MLE-based

* support for regressions of daily K versus discharge and/or velocity

* default specifications for every model

* a maturing user interface for fitting models (probably not quite fixed yet)

* convenience functions for calculating DO saturation concentrations, air 
pressure, depth, solar time, PAR, etc.

* functions for simulating data and error, for testing models with data having 
known underlying parameters

* two small datasets, courtesy of Bob Hall, for testing models with real data
