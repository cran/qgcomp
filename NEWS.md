
# qgcomp v2.18.7
## Major changes
- None

## Minor changes
- export .qgcomp_object function and .qgcomp_object_add function. These function facilitate transportable methods that can be used in external extensions to qgcomp (.e.g)

## Bug fixes
- bug fix for weights, cluster in Cox models


# qgcomp v2.18.4
## Major changes
- Added modelbound.ee function to generate confidence intervals for regression lines with estimating equation fits
- Revamped underlying code to plots to facilitate integration/overlays

## Minor changes
- Implementing a few functions from `qgcompint` package to support consistency across packages
- Documentation fixes (cch functions and missing plots/formulae in vignettes)

## Bug fixes
- fix bug in df.residual function that impacted integration with `mice` package


# qgcomp v2.17.3
## Major changes
- Added glm.tobit.noboot function to fit a Tobit model with qgcomp
 
## Minor changes
- Documentation improvements
- Code hygeine to prevent issues with edge cases

## Bug fixes
- None


# qgcomp v2.16.4
## Major changes
- Added glm.ee function to use estimating equations approach (allows non-linear and longitudinal models without bootstrapping)

## Minor changes
- Split vignette into two vignettes and added examples

## Bug fixes
- Mice function (mice.impute.leftcenslognorm) to allow empty prediction sets: fixed a bug that led to hidden failures resulting in very low values for imputations
- Fixed glm.ee for Tibbles, non-linear models (which never made it to CRAN release)


# qgcomp v2.15.2
## Major changes
- Added multinomial models (qgcomp.multinomial.noboot and qgcomp.multinomial.boot)
- Plot for multinomial models (qgcomp.multinomial.noboot only)
- Changed qgcomp.boot to qgcomp.glm.boot and qgcomp.noboot to qgcomp.glm.noboot for clarity. The old functions still work (for now).

## Minor changes
- Mice function (mice.impute.leftcenslognorm) to allow empty prediction sets 
- mice.impute.leftcenslognorm has been given an alias mice.impute.tobit
- Added simulation of multinomial outcomes
- Code-base hygiene

## Bug fixes
- Fixed some miscellaneous bugs in the `print` function for non-standard link functions.
- Fixing bug in model pointwisebound.noboot that gave non-sense confidence intervals

# qgcomp v2.10.0
## Major changes
- None

## Minor changes
- qgcomp.partials now allows quantile definitions based on the training and validation data, which treats the quantiles as fixed values across both datasets and leads to more stable results in small samples (set via .globalbreaks = TRUE in qgcomp.partials)

## Bug fixes
- Fixed major bug in qgcomp.partials: https://github.com/alexpkeil1/qgcomp/issues/28

# qgcomp v2.9.0
## Major changes
- None

## Minor changes
- Replace qgcomp.partials with re-factored version with new quantization defaults

## Bug fixes
- Fixed error https://github.com/alexpkeil1/qgcomp/issues/26
- Clarifed error message in https://github.com/alexpkeil1/qgcomp/issues/14
- Created new hidden functions for testing that will eventually replace qgcomp.partials


# qgcomp v2.8.6
## Major changes
- None

## Minor changes
- None

## Bug fixes
- Fixed misuse of "all.equal" in tests

# qgcomp v2.8.5
## Major changes
- None

## Minor changes
- Handling of parallel processing through futures: bringing in line with package recommendations and planned deprecations.
- Documentation changes

# qgcomp v2.8.0
## Major changes
- Exposed simulation functions to the user
- Enabled no-intercept models in qgcomp.boot and qgcomp.noboot

## Minor changes
- Parallel code now enables external setting of if(parplan) future::plan to facilitate broader parallelization schemes
- Documentation improvements
- Included FAQs in vignette

# qgcomp v2.7.0
## Major changes
- None

## Minor changes
- Changed underlying parallel code to adapt to deprecations from future package
- Documentation improvements
- Enabled by-observation limit of detection/right censoring point for mice.impute.leftcensorlognorm

## Bug fixes
- None

# qgcomp v2.6.0
## Major changes
- Added `cox.survcurve.boot` function to estimate survival curves from a qgcomp.cox.boot fit
- Refactored plot code (changes should be invisible to users except bug fixes)

## Minor changes
- documentation improvements
- Changed "pointwisebound.\*" and "modelbounds.\*" to give more sensible output

## Bug fixes
- Fix bug in tidy when using Cox model
- Bug fixes in bounds calculation with binomial outcomes, which lead to incorrect bounds on plots
- Fixed long running harmless and annoying print bug that confused "t" with "z"


# qgcomp v2.5.0
## Major changes
- Documentation improvements
- Added `split_data` function to assist with inference from data adaptive procedures (such as qgcomp.partials) and model selection

## Minor changes
- Removed "experimental" tags on Poisson family, Zero inflated models
- Added smarter "expnms" default to `predict.qgcompfit` function with new data

# qgcomp v2.4.0
## Major changes
- Documentation improvements
- Added hurdle models (qgcomp.hurdle.*boot functions), which address the excess in
zeroes in a distribution slightly differently from zero inflated models (based on
`pscl` package functions.)

## Bug fixes
- Package `broom` deprecates a dependency in an upcoming release. This dependency was removed
- Bug with pointwisebound.noboot when non-exposure variables are encoded as factors
- fixed bug in .split.cluster.data function that failed with most ID variable types
- fixed R 4.0.0 bug where MCsize was ignored
- fixed rare bug in hurdle/zi methods with non-quantized exposures, resulted in ridiculous runtimes if not errors

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
