# MLZ: Mean Length-based Z Estimators

## Introduction
MLZ is a package that facilitates data preparation and estimation of mortality with statistical diagnostics using the mean length-based mortality estimator and several extensions.

## Installation
The package can be installed using the `devtools` package:

```r
### install missing dependencies first
depends <- list("methods", "stats", "graphics", "grDevices", "dplyr", "gplots", "ggplot2", "reshape2", "numDeriv", "parallel", "msm", "Rcpp", "RcppArmadillo")
lapply(depends, function(x) if(!x %in% installed.packages()) install.packages(x))

devtools::install_github("quang-huynh/MLZ")
```

An overview of the main features of the package is described in the vignette:
```r
vignette("MLZ")
```

## Version
Currently, the package is in development in version 0.0.0.9000. Testing, reporting issues, and pull requests are welcome.
