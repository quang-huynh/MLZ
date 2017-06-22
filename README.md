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
