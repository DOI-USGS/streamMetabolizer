# streamMetabolizer

Uses inverse modeling to estimate aquatic metabolism (photosynthesis and respiration) from time series data on dissolved oxygen, water temperature, depth, and light.

## Status

| Branch | Linux | Windows | Test Coverage | USGS Status | DOI |
|--------|-------|---------|---------------|-------------|-----|
| master | [![master Build Status](https://travis-ci.org/USGS-R/streamMetabolizer.svg?branch=master)](https://travis-ci.org/USGS-R/streamMetabolizer/branches) | [![master Build status](https://ci.appveyor.com/api/projects/status/605tgcru05jdgb22/branch/master?svg=true)](https://ci.appveyor.com/project/aappling-usgs/streammetabolizer/branch/master) | [![master Coverage Status](https://coveralls.io/repos/github/USGS-R/streamMetabolizer/badge.svg?branch=master)](https://coveralls.io/github/USGS-R/streamMetabolizer?branch=master) | [![USGS Status](https://img.shields.io/badge/USGS-Research-blue.svg)](https://owi.usgs.gov/R/packages.html#research) | [![DOI](https://zenodo.org/badge/34403148.svg)](https://zenodo.org/badge/latestdoi/34403148) |
| develop| [![develop Build Status](https://travis-ci.org/USGS-R/streamMetabolizer.svg?branch=develop)](https://travis-ci.org/USGS-R/streamMetabolizer/branches) | [![develop Build status](https://ci.appveyor.com/api/projects/status/605tgcru05jdgb22/branch/develop?svg=true)](https://ci.appveyor.com/project/aappling-usgs/streammetabolizer/branch/develop) | [![develop Coverage Status](https://coveralls.io/repos/github/USGS-R/streamMetabolizer/badge.svg?branch=develop)](https://coveralls.io/github/USGS-R/streamMetabolizer?branch=develop) | | |

## Questions and bug reports

Please report bugs and ask questions on the Issues page:
[https://github.com/USGS-R/streamMetabolizer/issues](https://github.com/USGS-R/streamMetabolizer/issues)

We can address your issues fastest if they are:

* New - search past issues first to see if your question has already been answered. Reopen an old issue if your question was almost but not quite answered.

* Complete - share all relevant code and console output/errors/warnings, in the order they were run and produced on your computer. Also include the output from a call to `devtools::session_info()` to tell us about your computer's configuration.

* Reproducible - include all the data and code necessary for us/others to recreate the problem locally. It's fine to make up data if you can't share yours, as long as the problem still comes through.

* Minimal - what is the smallest amount of data and code you can use to demonstrate the problem? This is less essential than the others but improves communication and our response time.


### Code of Conduct

We want to encourage a warm, welcoming, and safe environment for contributing to this project. See the [code of conduct](https://github.com/USGS-R/streamMetabolizer/blob/develop/CONDUCT.md) for more information.

### Package Support

`streamMetabolizer` was developed 2015-2017 with support from the USGS Powell Center (through a working group on Continental Patterns of Stream Metabolism), the USGS NAWQA program, and the USGS Office of Water Information. Ongoing package work is unfunded and therefore limited, though still enthusiastic.

![USGS](http://usgs-r.github.io/images/usgs.png)

Follow `@USGS_R` on Twitter for updates on USGS R packages:

[![Twitter Follow](https://img.shields.io/twitter/follow/USGS_R.svg?style=social&label=Follow%20USGS_R)](https://twitter.com/USGS_R)

## Installation

The most stable and current version of this package can be installed with this R command:
```r
install.packages("streamMetabolizer", dependencies=TRUE, 
  repos=c("https://owi.usgs.gov/R","https://cran.rstudio.com"))
```

### For the adventurous

The in-development version of the package can be installed with the `devtools` package. 
We can make no guarantees about the stability of this version, 
but it might have new features that you'll like.
If you go this route, first install as above to install all the package dependencies.
You can then install the most cutting edge version of streamMetabolizer with this command:
```r
devtools::install_github("USGS-R/streamMetabolizer", ref="develop")
```

### Software dependencies

If you plan to use Bayesian models, you will need an up-to-date installation of
[Rtools](https://cran.r-project.org/bin/windows/Rtools/). Run
`devtools::find_rtools()` to make sure Rtools is ready to go. (Rtools is broadly
useful for R packages and might become a stronger dependency of
`streamMetabolizer` in the future.) Having Rtools installed will allow you to
install `rstan`, the package that `streamMetabolizer` relies on to run MCMC
models.

Bayesian models require the [`rstan`](http://mc-stan.org/interfaces/rstan.html) 
interface to [Stan](http://mc-stan.org/). Sometimes this is as simple as
installing Rtools and calling the above `install.packages` command, but other
times everything seems fine until you try to run a Bayesian model in
`streamMetabolizer`. Symptoms of an imperfect `rstan` installation are probably
diverse. Here's one we've seen:
```r
> bayes_fit <- metab(specs('bayes'), data=mydat)
Warning message:
In metab_fun(specs = specs, data = data, data_daily = data_daily,  :
  Modeling failed: argument is of length zero

> get_fit(bayes_fit)
...
$warnings
[1] "running command ''/Library/Frameworks/R.framework/Resources/bin/R' CMD config CXX 2>/dev/null' had status 1"

$errors
[1] "argument is of length zero"
```
In such cases you should refer to the detailed instructions on the `rstan` website for 
[Mac and Linux](https://github.com/stan-dev/rstan/wiki/Installing-RStan-on-Mac-or-Linux) 
or [Windows](https://github.com/stan-dev/rstan/wiki/Installing-RStan-on-Windows).

## Disclaimer

This software is in the public domain because it contains materials that originally came from the U.S. Geological Survey  (USGS), an agency of the United States Department of Interior. For more information, see the official USGS copyright policy at [https://www.usgs.gov/visual-id/credit_usgs.html#copyright](https://www.usgs.gov/visual-id/credit_usgs.html#copyright)

Although this software program has been used by the USGS, no warranty, expressed or implied, is made by the USGS or the U.S. Government as to the accuracy and functioning of the program and related program material nor shall the fact of distribution constitute any such warranty, and no responsibility is assumed by the USGS in connection therewith.

This software is provided "AS IS."

 [
    ![CC0](https://i.creativecommons.org/p/zero/1.0/88x31.png)
  ](https://creativecommons.org/publicdomain/zero/1.0/)
