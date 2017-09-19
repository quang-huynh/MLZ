# MLZ: Mean Length-based Z Estimators

## Introduction
MLZ is a package that facilitates data preparation and estimation of mortality with statistical diagnostics using the mean length-based mortality estimator and several extensions.

## Version
Currently, the package is in development, version 0.0.0.9001. Testing, reporting issues, and pull requests are welcome.

## Installation
In version 0.0.0.90001, models are written in TMB. As a result, Rtools is required if running the package on a Windows machine. Rtools can be downloaded here: https://cran.r-project.org/bin/windows/Rtools/index.html

The MLZ package can be installed using the `devtools` package:

```r
### install missing dependencies first
depends <- list("methods", "stats", "graphics", "grDevices", "dplyr", "gplots", "ggplot2", "reshape2", 
                "parallel", "TMB")
lapply(depends, function(x) if(!x %in% installed.packages()) install.packages(x))

devtools::install_github("quang-huynh/MLZ")
```

An overview of the main features of the package is described in the vignette:
```r
vignette("MLZ")
```
## Citations
Gedamke, T. and Hoenig, J.M. 2006. Estimating mortality from mean length data in nonequilibrium situations, with application to the assessment of goosefish. Transactions of the American Fisheries Society 135:476-487.

Huynh, Q.C, Gedamke, T., Hoenig, J.M, and Porch C. 2017. Multispecies Extensions to a Nonequilibrium Length-Based Mortality Estimator. Marine and Coastal Fisheries 9:68-78.

Huynh, Q.C., Gedamke, T., Porch, C.E., Hoenig, J.M., Walter, J.F, Bryan, M., and Brodziak, J. 2017. Estimating Total Mortality Rates from Mean Lengths and Catch Rates in Non-equilibrium Situations. Transactions of the American Fisheries Society 146:803-815.

Then, A.Y, Hoenig, J.M, and Huynh, Q.C. In revision. Estimating fishing and natural mortality rates, and catchability coefficient, from a series of observations on mean length and fishing effort. ICES Journal of Marine Science.
