# Unreleased

# OutSeekR 1.1.0 - 2026-01-15

## Major Changes
* **Breaking change**: Replaced GAMLSS-based distribution fitting with Gaussian Mixture Models (GMM) using the `mclust` package for more flexible and robust distribution modeling.
* **Breaking change**: Removed `distributions` output from `detect.outliers()` function as individual distribution fitting is no longer performed.
* **Breaking change**: Removed `distribution` parameters from `calculate.residuals()`, `simulate.null()`, and `outlier.detection.cosine()` functions.

## Improvements
* Enhanced null data simulation using GMM-based approach that better captures complex data distributions.

## Dependencies
* Added `mclust` package dependency for Gaussian mixture modeling.
* Removed dependencies on `gamlss` and `gamlss.dist` packages.

## Documentation
* Updated methodology descriptions in README and vignettes to reflect GMM-based approach.
* Revised examples to remove distribution-specific parameters and outputs.

## Function Changes
* `detect.outliers()`: Removed `distributions` from output list.
* `calculate.residuals()`: Removed `distribution` parameter, now uses GMM fitting.
* `simulate.null()`: Removed `x.distribution` and `r.distribution` parameters, now uses GMM sampling.
* `outlier.detection.cosine()`: Removed `distribution` parameter, now uses GMM-based quantiles.
* `calculate.p.values()`: Removed `x.distribution` parameter.
