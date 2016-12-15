# streamMetabolizer
stream metabolism R package


## Status

This package is in development. We are using it for our own early applications and welcome bold, flexible, resilient new users who can help us make the package better. Please contact us if you might be one of those.

| Name       | Status (develop branch)   |  Status (master branch) |
| :------------ |:-------------|:-------------| 
| Linux Build: | [![develop Build Status](https://travis-ci.org/USGS-R/streamMetabolizer.svg?branch=develop)](https://travis-ci.org/USGS-R/streamMetabolizer/branches)  | [![master Build Status](https://travis-ci.org/USGS-R/streamMetabolizer.svg?branch=master)](https://travis-ci.org/USGS-R/streamMetabolizer/branches) |
| Windows Build: | [![develop Build status](https://ci.appveyor.com/api/projects/status/605tgcru05jdgb22/branch/develop?svg=true)](https://ci.appveyor.com/project/aappling-usgs/streammetabolizer/branch/develop) | [![master Build status](https://ci.appveyor.com/api/projects/status/605tgcru05jdgb22/branch/master?svg=true)](https://ci.appveyor.com/project/aappling-usgs/streammetabolizer/branch/master) |  
| Package Tests: | [![develop Coverage Status](https://coveralls.io/repos/github/USGS-R/streamMetabolizer/badge.svg?branch=develop)](https://coveralls.io/github/USGS-R/streamMetabolizer?branch=develop) | [![master Coverage Status](https://coveralls.io/repos/github/USGS-R/streamMetabolizer/badge.svg?branch=master)](https://coveralls.io/github/USGS-R/streamMetabolizer?branch=master) |  


## Installation

### Recommended

The most stable+current version of this package can be installed with this R command:
```r
install.packages("streamMetabolizer", dependencies=TRUE, 
  repos=c("https://owi.usgs.gov/R","https://cran.rstudio.com"))
```
and updated with this command:
```r
update.packages(oldPkgs=c("streamMetabolizer","unitted"),
  dependencies=TRUE, repos=c("https://owi.usgs.gov/R", "https://cran.rstudio.com"))
```

### For the adventurous

The in-development version of the package can be installed with the `devtools` package. 
We can make no guarantees about the stability of this version, 
but it might have new features that you'll like.
If you go this route, you will need to install the package dependencies separately, like this:
```r
install.packages(
  c("LakeMetabolizer","unitted","dplyr","lazyeval","lubridate","magrittr",
    "tidyr","chron","dygraphs","ggplot2","RCurl","rstan","XML","xts"),
  repos=c("https://owi.usgs.gov/R","https://cran.rstudio.com"))
```
You can then install the most cutting edge version of streamMetabolizer with this command:
```r
devtools::install_github("USGS-R/streamMetabolizer", ref="develop")
```

### Software dependencies

If you plan to use Bayesian models, you will need an up-to-date installation of [Rtools](http://cran.r-project.org/bin/windows/Rtools/). Run `devtools::find_rtools()` to make sure Rtools is ready to go. (Rtools is broadly useful for R packages and might become a stronger dependency of `streamMetabolizer` in the future.) Having Rtools installed will allow you to install rstan, the package that `streamMetabolizer` relies on to run MCMC models.


## Disclaimer

This software is in the public domain because it contains materials that originally came from the U.S. Geological Survey, an agency of the United States Department of Interior. For more information, see the official USGS copyright policy at [http://www.usgs.gov/visual-id/credit_usgs.html#copyright](http://www.usgs.gov/visual-id/credit_usgs.html#copyright)

This information is preliminary or provisional and is subject to revision. It is being provided to meet the need for timely best science. The information has not received final approval by the U.S. Geological Survey (USGS) and is provided on the condition that neither the USGS nor the U.S. Government shall be held liable for any damages resulting from the authorized or unauthorized use of the information. Although this software program has been used by the USGS, no warranty, expressed or implied, is made by the USGS or the U.S. Government as to the accuracy and functioning of the program and related program material nor shall the fact of distribution constitute any such warranty, and no responsibility is assumed by the USGS in connection therewith.

This software is provided "AS IS."


 [
    ![CC0](http://i.creativecommons.org/p/zero/1.0/88x31.png)
  ](http://creativecommons.org/publicdomain/zero/1.0/)
