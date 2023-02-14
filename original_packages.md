# Overview of packages

Felix sent the code in January, 2023 along with a list of 28 packages. Robert successfully ran the analysis in RStudio. This file describes the way in which each package needs to be loaded into a `conda` environment, and (eventually) also a Binder image.

## Original List

This list is the original set of 28 packages that were provided by Felix and loaded in RStudio by importing all as `library*)`.

- base
- copula
- datasets
- DISTRIB
- distributions3
- fitdistrplus
- graphics
- grDevices
- invgamma
- lattice
- lestat
- lme4
- lmomco
- Lmoments
- lsei
- MASS
- Matrix
- methods
- mvtnorm
- npsurv
- pracma
- rvinecopulib
- stats
- stats4
- survival
- utils
- VineCopula
- xlsx

## Available in `conda` 

`conda` has a set of "essential" packages for R, which can be loaded in a `requirements.txt` file (with `r-*` prefix):

- r-base
- r-copula
- r-DISTRIB
- r-distributions3
- r-invgamma
- r-lattice
- r-lestat
- r-lme4
- r-MASS
- r-Matrix
- r-mvtnorm
- r-pracma
- r-survival
- r-xlsx

## Not available in `conda` 

These were not available in `conda`:
- stats4
- lsei
- fitdistrplus
- utils
- methods
- Lmoments
- datasets
- VineCopula
- grDevices
- stats
- lmomco
- rvinecopulib
- npsurv
- graphics

### Included in basic R installation

Of the packages not available in `conda`, some are part of the 15 packages distributed with base R:
- stats4
- utils
- methods
- datasets
- grDevices
- stats
- graphics

## Direct download

These are the packages that need to be obtained manually from CRAN as `install.packages("npsurv", repos="<URL>")`
- lsei
- fitdistrplus
- Lmoments
- VineCopula
- lmomco
- rvinecopulib
- npsurv
