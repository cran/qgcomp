# qgcomp v1.1.0
## Major changes
- Added bootstrap covariance matrix to output of qgcomp.boot, which allows plotting of MSM intervals
- Plots now include confidence bands, when appropriate, and offer more control
- Fixed standardization in data
- Added Cox proportional hazards variants of qgcomp models

## Bug fixes
- Fixed issue where MSM fit would default to 4 levels if `q` is unspecified and breaks are user defined (rare, major)

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
