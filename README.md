# MLZ: Mean Length-based Z Estimators

[![Build Status](https://travis-ci.org/quang-huynh/MLZ.svg?branch=master)](https://travis-ci.org/quang-huynh/MLZ)
[![CRAN_Status_Badge](http://www.r-pkg.org/badges/version/MLZ)](https://CRAN.R-project.org/package=MLZ)

## Introduction
MLZ is a R package that facilitates data preparation and estimation of mortality with statistical diagnostics using the mean length-based mortality estimator and several extensions.

After installation, an overview of the main features of the package can be found in the vignette:
```r
library(MLZ)
vignette("MLZ")
```

## Version
Currently, the package is in active development, version 0.1.1. Testing, reporting issues, and pull requests are welcome.

## Installation
The package is available on CRAN:
```r
install.packages("MLZ")
```

The development source code can also be downloaded and installed with the `devtools` package:
```r
devtools::install_github("quang-huynh/MLZ")
```
Since the estimation models are written in TMB, Rtools is required on a Windows machine for installation from the source package. Rtools can be downloaded here: https://cran.r-project.org/bin/windows/Rtools/index.html.

## Citations
Gedamke, T. and Hoenig, J.M. 2006. Estimating mortality from mean length data in nonequilibrium situations, with application to the assessment of goosefish. Transactions of the American Fisheries Society 135:476-487.

Huynh, Q.C, Gedamke, T., Hoenig, J.M, and Porch C. 2017. Multispecies Extensions to a Nonequilibrium Length-Based Mortality Estimator. Marine and Coastal Fisheries 9:68-78.

Huynh, Q.C., Gedamke, T., Porch, C.E., Hoenig, J.M., Walter, J.F, Bryan, M., and Brodziak, J. 2017. Estimating Total Mortality Rates from Mean Lengths and Catch Rates in Non-equilibrium Situations. Transactions of the American Fisheries Society 146:803-815.

Then, A.Y, Hoenig, J.M, and Huynh, Q.C. 2018. Estimating fishing and natural mortality rates, and catchability coefficient, from a series of observations on mean length and fishing effort. ICES Journal of Marine Science 75: 610-620.
