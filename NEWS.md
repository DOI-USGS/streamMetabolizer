# streamMetabolizer 0.9.5.1

* in metab_Kmodel, now avoiding negative weights

# streamMetabolizer 0.9.5

* Bug fixes and error prevention

# streamMetabolizer 0.9.4

* Now automatically checks for available updates when you attach the package

* Improved units handling in convert_k600/kGAS

# streamMetabolizer 0.9.2

* Hierarchical constraints on K600 are now available! Options are 'normal', 
'linear', and 'binned'; see the details section on pool_K600 in ?mm_name and the
description of parameters starting with K600_daily in ?specs.

* Interface change: specs lists now print more prettily and have class 'specs' 
(though they're still fundamentally just lists)

* Vignette: see vignette('getstarted')

# streamMetabolizer 0.9.0

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

# streamMetabolizer 0.8.0

* Major interface change (renamed variable) to clarify types of time: solar.time
(mean solar time), app.solar.time (apparent solar time), local.time (time in 
local time zone). Metabolism models now accept solar.time rather than 
local.time, though it's still possible to pass in local time but just call it 
solar.time (as long as you don't have daylight savings time).

# streamMetabolizer 0.7.3

* Remove calc\_schmidt because it is never used

# streamMetabolizer 0.7.2

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