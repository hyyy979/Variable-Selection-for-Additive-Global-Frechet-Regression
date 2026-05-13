# Variable Selection for Additive Global Fréchet Regression

This repository contains reproducible example code for the paper:

> "Variable Selection for Additive Global Fréchet Regression"

The repository includes four simulation examples corresponding to:
- distribution-valued responses
- SPD-valued responses
- linear settings
- nonlinear settings

## Files

`distribution_linear_example.R`  
Distribution-valued response under the linear setting.

`distribution_nonlinear_example.R`  
Distribution-valued response under the nonlinear setting.

`spd_linear_example.R`  
SPD-valued response under the linear setting.

`spd_nonlinear_example.R`  
SPD-valued response under the nonlinear setting.

## Methods

All examples implement:
- Elastic Net regularization
- Adaptive SCAD refinement
- ADMM optimization
- Validation-based tuning parameter selection

## Required Packages

```r
library(MASS)
library(Matrix)
