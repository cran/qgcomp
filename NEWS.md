# qgcomp v2.3.0
## Major changes
- Documentation improvements
- Added ability to perform weighted analysis for glm, survival, and zero inflated objects
- Implement sample splitting to estimate partial effects via `qgcomp.partials`.

## Bug fixes
- Fixed bug in print function that missed linear fits
- Fixed visual bug in plots introduced after updating ggplot 3.3.0  (instituted forced update of ggplot2 for users, sorry!)

# qgcomp v2.2.2
## Major changes
- Documentation improvements
- Added `pointwisebound.boot()` and  `pointwisebound.noboot()` function to get confidence intervals for pointwise comparisons (e.g. mean outcome at all exposures below first quartile compared to mean outcome at all exposures above the third quartile. Used after qgcomp fits.
- Added support for Poisson link in `qgcomp.boot()` and `qgcomp.noboot()` functions to fit quantile g-computation for count data 
- Added `generics` and `tibble` package dependencies and optional installs of `broom`
- Added functions to facilitate multiple imputation through the `mice` package.
- Added `weight` parameters to base, zi, and cox (not yet used) in preparation for v2.3.0, which will implement weighting
- Added `vc_comb()` to better calculate covariance matrices.

## Bug fixes
- Fixed an error in confidence bounds in plot functions introduced in v2.0.0 (binary outcomes in `qgcomp.boot`).
- `id` parameter caused intermittent error when used in  `qgcomp.boot`
- `MCsize` parameter caused intermittent error when used in  `qgcomp.boot`
-  `vcov()` didn't yield a covariance matrix when applied to `qgcomp.*.noboot` objects

# qgcomp v2.1.2
## Major changes
- Documentation improvements
- Added `qgcomp.zi.boot` function to fit non-linear/non-additive versions of quantile g-computation for zero-inflated Poisson or negative binomial models 
- Added two functions to estimate confidence intervals for regression line fits and pointwise comparisons: modelbound.boot, pointwisebound.boot 

## Bug fixes
- MAJOR: Fixed error in confidence bounds in plot functions introduced in v2.0.0

# qgcomp v2.0.0
## Major changes
- Documentation improvements
- Added `qgcomp.zi.noboot` function to fit (linear/additive only) versions of quantile g-computation for zero-inflated Poisson or negative binomial models (adds additional dependency on `pscl` package)
- Major version change to account for possible breaking of old code in v.1.3.0

## Bug fixes
- MAJOR: Fixed error in print function that printed incorrect confidence intervals introduced in v1.3.0

# qgcomp v1.3.0
## Major changes
- Documentation improvements
- Improved error messages
- Added intercept to main output (this may break some code that depends on qgcomp)
- Added pointwise confidence bounds and better control over plot output for `qgcomp.boot` fits.

## Bug fixes
- Fixed fatal error for plotting output of `qgcomp.cox.boot` when model was adjusted for covariates
- Fixed variable typing issue that arose when fitting `qgcomp.boot` with a single exposure variable with non-linear fit


# qgcomp v1.2.0
## Major changes
- Added bootstrapped version of Cox proportional hazards model, which estimates parameters of a **marginal structural Cox model** that uses a simulation-based approach to characterize the population average effect of increasing all exposures of interest by one quantile. Survival analysis now available from `qgcomp.cox.boot` and `qgcomp.cox.noboot` functions.
- Added (semi-)**Bayesian** approaches to linear, logistic models via the `arm` package. Available for `qgcomp.boot`, `qgcomp.noboot`, and `qgcomp` functions
- Allow parallel processing of bootstrap variance estimation via the `future` and `future.apply` packages
- Added default plotting of survival curves for the `qgcomp.cox.boot` function
- Miscellaneous error checking to improve readability of error messages

## Bug fixes
- Better guessing for `qgcomp` function

# qgcomp v1.1.0
## Major changes
- Added bootstrap covariance matrix to output of qgcomp.boot, which allows plotting of MSM intervals
- Plots now include confidence bands, when appropriate, and offer more control
- Fixed standardization in data
- Added Cox proportional hazards variants of qgcomp models
- Allow tibbles to be used as data objects

## Bug fixes
- Fixed issue where MSM fit would default to 4 levels if `q` is unspecified and breaks are user defined (rare, major)
- Fixed error that occurred when a tibble was used as the data object, rather than a data frame (common, minor)

# qgcomp v1.0.0
## Major changes
- First release on CRAN, benchmark v0.8.14 as release version 1.0

## Bug fixes
- v0.8.12: Breaks now properly reported as NULL using `qgcomp.boot` when q = NULL, rather than causing a fatal error (rare, minor)
- v0.8.13: Fixed error if specifying nonlinear fit in `qgcomp.boot` when q = NULL by specifying fixed values for "interventions" and warning users (rare, minor)
- v0.8.14: Documentation fixes for CRAN release

# qgcomp v0.8.13
## Major changes

## Bug fixes
- v0.8.12: Breaks now properly reported as NULL using `qgcomp.boot` when q = NULL, rather than causing a fatal error
- v0.8.13: Fixed error if specifying nonlinear fit in `qgcomp.boot` when q = NULL by specifying fixed values for "interventions" and warning users

# qgcomp v0.8.11
## Major changes
- Code finalized for release 
