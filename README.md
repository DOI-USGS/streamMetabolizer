# streamMetabolizer

## Models for Estimating Aquatic Photosynthesis and Respiration

Alison P. Appling, Robert O. Hall, Jr., Maite Arroita, and Charles B. Yackulic

The streamMetabolizer R package uses inverse modeling to estimate aquatic
photosynthesis and respiration (collectively, metabolism) from time series
data on dissolved oxygen, water temperature, depth, and light. The package
assists with data preparation, handles data gaps during modeling, and
provides tabular and graphical reports of model outputs. Several
time-honored methods are implemented along with many promising new variants
that produce more accurate and precise metabolism estimates.

This package has been described, with special focus on the Bayesian model options, by
[Appling et al. 2018](https://doi.org/10.1002/2017JG004140).

Recommended citation: To see the recommended citation for this package, please run `citation('streamMetabolizer')` at the R prompt.

## Getting started

See http://usgs-r.github.io/streamMetabolizer for tutorials on package installation, getting started, and customizing your metabolism models.


## Package status

| Branch | Linux | Windows | Test Coverage | USGS Status | DOI |
|--------|-------|---------|---------------|-------------|-----|
| master | [![master Build Status](https://travis-ci.org/USGS-R/streamMetabolizer.svg?branch=master)](https://travis-ci.org/USGS-R/streamMetabolizer/branches) | [![master Build status](https://ci.appveyor.com/api/projects/status/605tgcru05jdgb22/branch/master?svg=true)](https://ci.appveyor.com/project/aappling-usgs/streammetabolizer/branch/master) | [![master Coverage Status](https://coveralls.io/repos/github/USGS-R/streamMetabolizer/badge.svg?branch=master)](https://coveralls.io/github/USGS-R/streamMetabolizer?branch=master) | [![USGS Status](https://img.shields.io/badge/USGS-Research-blue.svg)](https://owi.usgs.gov/R/packages.html#research) | |


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

`streamMetabolizer` was developed 2015-2017 with support from the USGS Powell Center (through a working group on Continental Patterns of Stream Metabolism), the USGS NAWQA program, and the USGS Office of Water Information. Ongoing package work is unfunded and therefore limited.

[![USGS](http://usgs-r.github.io/images/usgs.png)](https://www.usgs.gov/)

Follow [`@USGS_R` on Twitter](https://twitter.com/USGS_R) for updates on USGS R packages.

## Disclaimer

This software is in the public domain because it contains materials that originally came from the U.S. Geological Survey  (USGS), an agency of the United States Department of Interior. For more information, see the official USGS copyright policy at [https://www.usgs.gov/visual-id/credit_usgs.html#copyright](https://www.usgs.gov/visual-id/credit_usgs.html#copyright)

Although this software program has been used by the USGS, no warranty, expressed or implied, is made by the USGS or the U.S. Government as to the accuracy and functioning of the program and related program material nor shall the fact of distribution constitute any such warranty, and no responsibility is assumed by the USGS in connection therewith.

This software is provided "AS IS."

 [
    ![CC0](https://i.creativecommons.org/p/zero/1.0/88x31.png)
  ](https://creativecommons.org/publicdomain/zero/1.0/)
