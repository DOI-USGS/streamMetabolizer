# streamMetabolizer: Models for Estimating Aquatic Photosynthesis and Respiration

``` diff
! In summer or fall 2023, this package will move from
! https://github.com/USGS-R/streamMetabolizer to
! https://github.com/DOI-USGS/streamMetabolizer.
! Please update your links accordingly.
```

The `streamMetabolizer` R package uses inverse modeling to estimate
aquatic photosynthesis and respiration (collectively, metabolism) from
time series data on dissolved oxygen, water temperature, depth, and
light. The package assists with data preparation, handles data gaps
during modeling, and provides tabular and graphical reports of model
outputs. Several time-honored methods are implemented along with many
promising new variants that produce more accurate and precise metabolism
estimates.

This package has been described, with special focus on the Bayesian
model options, by [Appling et
al. 2018a](https://doi.org/10.1002/2017JG004140). An application to 356
streams across the U.S. is described in [Appling et
al. 2018b](https://doi.org/10.1038/sdata.2018.292).

> Appling, A. P., Hall, R. O., Yackulic, C. B., & Arroita, M. (2018a).
> Overcoming equifinality: Leveraging long time series for stream
> metabolism estimation. Journal of Geophysical Research:
> Biogeosciences, 123(2), 624–645.
> <https://doi.org/10.1002/2017JG004140>

> Appling, A. P., Read, J. S., Winslow, L. A., Arroita, M., Bernhardt,
> E. S., Griffiths, N. A., Hall, R. O., Harvey, J. W., Heffernan, J. B.,
> Stanley, E. H., Stets, E. G., & Yackulic, C. B. (2018b). The metabolic
> regimes of 356 rivers in the United States. Scientific Data, 5(1),
> 180292. <https://doi.org/10.1038/sdata.2018.292>

To see the recommended citation for this package, please run
`citation('streamMetabolizer')` at the R prompt.

``` r
citation('streamMetabolizer')
## 
## To cite streamMetabolizer in publications, please use:
## 
##   Appling, Alison P., Robert O. Hall, Charles B. Yackulic, and Maite
##   Arroita. “Overcoming Equifinality: Leveraging Long Time Series for
##   Stream Metabolism Estimation.” Journal of Geophysical Research:
##   Biogeosciences 123, no. 2 (February 2018): 624–45.
##   https://doi.org/10.1002/2017JG004140.
## 
## A BibTeX entry for LaTeX users is
## 
##   @Article{,
##     author = {Alison P. Appling and Robert O. {Hall Jr.} and Charles B. Yackulic and Maite Arroita},
##     title = {Overcoming Equifinality: Leveraging Long Time Series for Stream Metabolism Estimation},
##     journal = {Journal of Geophysical Research: Biogeosciences},
##     year = {2018},
##     volume = {123},
##     number = {2},
##     doi = {10.1002/2017JG004140},
##     url = {https://github.com/USGS-R/streamMetabolizer},
##   }
```

## Installation

To install the `streamMetabolizer` package, use the `remotes` package
(running `install.packages('remotes')` first if needed). To use
`remotes::install_github()` it is convenient to set a [GitHub Personal
Access Token
(PAT)](https://docs.github.com/en/authentication/keeping-your-account-and-data-secure/managing-your-personal-access-tokens).
There are [several
methods](https://usethis.r-lib.org/articles/git-credentials.html) for
setting your PATs within R; the simplest is to call
\`Sys.setenv(GITHUB_PAT=“yyyy”), replacing yyyy with the PAT you
established on the GitHub website.

You may first need to install the `unitted` dependency:

``` r
remotes::install_github('appling/unitted')
```

You can then install the most cutting edge version of
`streamMetabolizer` with this command:

``` r
remotes::install_github(
  "USGS-R/streamMetabolizer", # soon to be "DOI-USGS/streamMetabolizer"
  build_vignettes = TRUE)
```

### Software dependencies for Bayesian models

The major dependency for Bayesian models is the `rstan` package, and
installation of that package is rarely as simple as a call to
`install.packages()`. Start at the [rstan wiki
page](https://github.com/stan-dev/rstan/wiki) for the most up-to-date
installation instructions, which differ by operating system.

## Getting started

After installing and loading `streamMetabolizer`, run `vignette()` in R
to see tutorials on getting started and customizing your metabolism
models.

``` r
vignette(package='streamMetabolizer')
## displays a list of available vignettes

vignette('get_started', package='streamMetabolizer')
## displays an html or pdf rendering of the 'get_started' vignette
```

You can also view pre-built html versions of these vignettes in the
“inst/doc” folder in the source code, e.g.,
[inst/doc/get_started.html](https://github.com/USGS-R/streamMetabolizer/blob/main/inst/doc/get_started.html),
which you can download and then open in a browser.

## Development and Maintenance Status

`streamMetabolizer` is a USGS Archive Research Package: [![USGS
Status](https://img.shields.io/badge/USGS-Research-blue.svg)](https://owi.usgs.gov/R/packages.html#research)

Project funding has ended and our maintenance time is limited, but we do
attempt to provide bug fixes and lightweight support as we are able.
Submit questions or suggestions to
<https://github.com/USGS-R/streamMetabolizer/issues>.

## Contributing

We want to encourage a warm, welcoming, and safe environment for
contributing to this project. See
[CODE_OF_CONDUCT.md](https://github.com/USGS-R/streamMetabolizer/blob/main/CODE_OF_CONDUCT.md)
for more information.

For technical details on how to contribute, see
[CONTRIBUTING.md](https://github.com/USGS-R/streamMetabolizer/blob/main/CONTRIBUTING.md)

### Development History

`streamMetabolizer` was developed 2015-2018 with support from the USGS
Powell Center (through a working group on Continental Patterns of Stream
Metabolism), the USGS National Water Quality Program, and the USGS
Office of Water Information.

## Model Archive

The following version of R and package dependencies were used most
recently to pass the embedded tests within this package. There is no
guarantee of reproducible results using future versions of R or updated
versions of package dependencies; however, we aim to test and update
future modeling environments.

<!-- Run and paste manually after edits, only when tests pass locally -->

``` r
sessioninfo::session_info()

## ─ Session info ───────────────────────────────────────────────────────────────────────────────────
##  setting  value
##  version  R version 4.2.3 (2023-03-15)
##  os       macOS Ventura 13.4.1
##  system   x86_64, darwin17.0
##  ui       RStudio
##  language (EN)
##  collate  en_US.UTF-8
##  ctype    en_US.UTF-8
##  tz       America/New_York
##  date     2023-07-02
##  rstudio  2023.06.0+421 Mountain Hydrangea (desktop)
##  pandoc   3.1.1 @ /Applications/RStudio.app/Contents/Resources/app/quarto/bin/tools/ (via rmarkdown)
## 
## ─ Packages ───────────────────────────────────────────────────────────────────────────────────────
##  package           * version  date (UTC) lib source
##  cli                 3.6.1    2023-03-23 [1] CRAN (R 4.2.0)
##  deSolve             1.35     2023-03-12 [1] CRAN (R 4.2.0)
##  digest              0.6.32   2023-06-26 [1] CRAN (R 4.2.0)
##  dplyr               1.1.2    2023-04-20 [1] CRAN (R 4.2.0)
##  evaluate            0.21     2023-05-05 [1] CRAN (R 4.2.0)
##  fansi               1.0.4    2023-01-22 [1] CRAN (R 4.2.0)
##  fastmap             1.1.1    2023-02-24 [1] CRAN (R 4.2.0)
##  generics            0.1.3    2022-07-05 [1] CRAN (R 4.2.0)
##  glue                1.6.2    2022-02-24 [1] CRAN (R 4.2.0)
##  htmltools           0.5.5    2023-03-23 [1] CRAN (R 4.2.0)
##  knitr               1.43     2023-05-25 [1] CRAN (R 4.2.0)
##  LakeMetabolizer     1.5.5    2022-11-15 [1] CRAN (R 4.2.0)
##  lazyeval            0.2.2    2019-03-15 [1] CRAN (R 4.2.0)
##  lifecycle           1.0.3    2022-10-07 [1] CRAN (R 4.2.0)
##  lubridate           1.9.2    2023-02-10 [1] CRAN (R 4.2.0)
##  magrittr            2.0.3    2022-03-30 [1] CRAN (R 4.2.0)
##  pillar              1.9.0    2023-03-22 [1] CRAN (R 4.2.0)
##  pkgconfig           2.0.3    2019-09-22 [1] CRAN (R 4.2.0)
##  plyr                1.8.8    2022-11-11 [1] CRAN (R 4.2.0)
##  purrr               1.0.1    2023-01-10 [1] CRAN (R 4.2.0)
##  R6                  2.5.1    2021-08-19 [1] CRAN (R 4.2.0)
##  Rcpp                1.0.10   2023-01-22 [1] CRAN (R 4.2.0)
##  rLakeAnalyzer       1.11.4.1 2019-06-09 [1] CRAN (R 4.2.0)
##  rlang               1.1.1    2023-04-28 [1] CRAN (R 4.2.0)
##  rmarkdown           2.22     2023-06-01 [1] CRAN (R 4.2.0)
##  rstudioapi          0.14     2022-08-22 [1] CRAN (R 4.2.0)
##  sessioninfo         1.2.2    2021-12-06 [1] CRAN (R 4.2.0)
##  streamMetabolizer * 0.12.1   2023-07-02 [1] local
##  tibble              3.2.1    2023-03-20 [1] CRAN (R 4.2.0)
##  tidyr               1.3.0    2023-01-24 [1] CRAN (R 4.2.0)
##  tidyselect          1.2.0    2022-10-10 [1] CRAN (R 4.2.0)
##  timechange          0.2.0    2023-01-11 [1] CRAN (R 4.2.0)
##  unitted             0.2.9    2023-06-05 [1] Github (appling/unitted@d1f1172)
##  utf8                1.2.3    2023-01-31 [1] CRAN (R 4.2.0)
##  vctrs               0.6.3    2023-06-14 [1] CRAN (R 4.2.0)
##  xfun                0.39     2023-04-20 [1] CRAN (R 4.2.0)
##  yaml                2.3.7    2023-01-23 [1] CRAN (R 4.2.0)
## 
##  [1] /Library/Frameworks/R.framework/Versions/4.2/Resources/library
```

## Disclaimer

This software is preliminary or provisional and is subject to revision.
It is being provided to meet the need for timely best science. The
software has not received final approval by the U.S. Geological Survey
(USGS). No warranty, expressed or implied, is made by the USGS or the
U.S. Government as to the functionality of the software and related
material nor shall the fact of release constitute any such warranty. The
software is provided on the condition that neither the USGS nor the U.S.
Government shall be held liable for any damages resulting from the
authorized or unauthorized use of the software.
