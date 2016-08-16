# streamMetabolizer
stream metabolism R package


## Status

This package is in development. We are using it for our own early applications and welcome bold, flexible, resilient new users who can help us make the package better. Please contact us if you might be one of those.

| Name       | Status (develop branch)   |  Status (master branch) |
| :------------ |:-------------|:-------------| 
| Linux Build: | [![develop Build Status](https://travis-ci.org/USGS-R/streamMetabolizer.svg?branch=develop)](https://travis-ci.org/USGS-R/streamMetabolizer/branches)  | [![master Build Status](https://travis-ci.org/USGS-R/streamMetabolizer.svg?branch=master)](https://travis-ci.org/USGS-R/streamMetabolizer/branches) |
| Windows Build: | [![develop Build status](https://ci.appveyor.com/api/projects/status/n2u0tpmkaetj7kjp/branch/develop?svg=true)](https://ci.appveyor.com/project/jread-usgs/streammetabolizer/branch/develop) | [![master Build status](https://ci.appveyor.com/api/projects/status/n2u0tpmkaetj7kjp/branch/master?svg=true)](https://ci.appveyor.com/project/jread-usgs/streammetabolizer/branch/master) |  
| Package Tests: | [![develop Coverage Status](https://coveralls.io/repos/github/USGS-R/streamMetabolizer/badge.svg?branch=develop)](https://coveralls.io/github/USGS-R/streamMetabolizer?branch=develop) | [![master Coverage Status](https://coveralls.io/repos/github/USGS-R/streamMetabolizer/badge.svg?branch=master)](https://coveralls.io/github/USGS-R/streamMetabolizer?branch=master) |  
| Priorities: | [![Issues Ready to Address](https://badge.waffle.io/USGS-R/streamMetabolizer.png?label=ready&title=Ready)](https://waffle.io/USGS-R/streamMetabolizer) [![Issues in Progress](https://badge.waffle.io/USGS-R/streamMetabolizer.png?label=In%20Progress&title=In%20Progress)](https://waffle.io/USGS-R/streamMetabolizer)| | 


## Installation

### Recommended

The most stable+current version of this package can be installed with this R command:
```r
install.packages("streamMetabolizer", dependencies=TRUE, 
  repos=c("http://owi.usgs.gov/R","https://cran.rstudio.com"))
```
and updated with this command:
```r
update.packages(oldPkgs=c("streamMetabolizer","unitted"),
  dependencies=TRUE, repos=c("http://owi.usgs.gov/R", "https://cran.rstudio.com"))
```

### For the adventurous

The in-development version of the package can be installed with the `devtools` package. 
We can make no guarantees about the stability of this version, 
but it might have new features that you'll like.
If you go this route, you will need to install the package dependencies separately, like this:
```r
install.packages(
  c("LakeMetabolizer","unitted","dplyr","lazyeval","lubridate","magrittr",
    "tidyr","chron","dygraphs","ggplot2","RCurl","runjags","rstan","XML","xts"),
  repos=c("http://owi.usgs.gov/R","https://cran.rstudio.com"))
```
You can then install the most cutting edge version of streamMetabolizer with this command:
```r
devtools::install_github("USGS-R/streamMetabolizer", ref="develop")
```

### Software dependencies

If you plan to use Bayesian models with 
`engine='jags'`, you will also need to install [JAGS](http://mcmc-jags.sourceforge.net/). (JAGS may be deprecated in the future because so far we've seen better performance with Stan. You should be able to switch seamlessly between the two MCMC engines.)

If you plan to use Bayesian models with 
`engine='stan'`, you will need an up-to-date installation of [Rtools](http://cran.r-project.org/bin/windows/Rtools/). Run `devtools::find_rtools()` to make sure Rtools is ready to go. (Rtools is broadly useful for R packages and might become a stronger dependency of `streamMetabolizer` in the future.)


## Disclaimer

This software is in the public domain because it contains materials that originally came from the U.S. Geological Survey, an agency of the United States Department of Interior. For more information, see the official USGS copyright policy at [http://www.usgs.gov/visual-id/credit_usgs.html#copyright](http://www.usgs.gov/visual-id/credit_usgs.html#copyright)

This information is preliminary or provisional and is subject to revision. It is being provided to meet the need for timely best science. The information has not received final approval by the U.S. Geological Survey (USGS) and is provided on the condition that neither the USGS nor the U.S. Government shall be held liable for any damages resulting from the authorized or unauthorized use of the information. Although this software program has been used by the USGS, no warranty, expressed or implied, is made by the USGS or the U.S. Government as to the accuracy and functioning of the program and related program material nor shall the fact of distribution constitute any such warranty, and no responsibility is assumed by the USGS in connection therewith.

This software is provided "AS IS."


 [
    ![CC0](http://i.creativecommons.org/p/zero/1.0/88x31.png)
  ](http://creativecommons.org/publicdomain/zero/1.0/)
