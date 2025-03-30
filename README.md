Documentation For hdjmtTest
================

<!-- README.md is generated from README.Rmd. Please edit that file -->

# Introduction

The `hdjmtTest` package performs high-dimensional mediator screening and
testing in the presence of high correlation among mediators using
Ridge-HOLP and Approximate Orthogonalization methods. In addition to
accounting for correlation, the method is designed to detect mediators
that exhibit joint effects. This vignette demonstrates how to use the
primary functions with a simulated dataset, `ExampleData`.

# Structure of ExampleData

The example data provided with the package contains:

- information on `n = 200` samples and `p = 2000` potential active
  mediators.
- **Variables**:
  - `M`: A 200x2000 matrix of mediators generated using the compound
    symmetry covariance structure to introduce high correlation among
    the mediators, with $\rho= 0.8$.
  - `x`: A vector of length 200 representing exposures.
  - `y`: A vector of length 200 representing outcomes.
  - `alp_vec`: A parameter vector of length 200 that relates the
    exposure variable to the mediators. The first 8 values are
    `-0.5127127, 0.6597036, 0.6756640, 0.5235137, 0.9305369, -0.9827865, -0.8941141, -0.9230220`
    with the rest being zeros.
  - `beta_vec`: A parameter vector of length 200 that relates the
    mediators to the outcome variable. The first 12 values are
    `-0.8033093, 0.9360975, 0.8185305, -0.7951502, -0.8783739, 0.8940459, -0.6911509, -0.8524771, -0.6812097, -0.8285034, -0.5986530, -0.9639383`
    with the rest being zeros.
- **Active Mediators**: given the non-zero values in `alp_vec` and
  `beta_vec`, only the first 8 mediators are truly active.

# How to install package from github

``` r
#library(devtools)

#devtools::install_github('samuelOsarfo/hdjmtTest')
```

If package “qvalue” is not found, please first install “qvalue” package
through Bioconductor:
<https://www.bioconductor.org/packages/release/bioc/html/qvalue.html>

``` r
#if (!require("BiocManager", quietly = TRUE))
#    install.packages("BiocManager")
#BiocManager::install("qvalue")
```

# Summary of Functions

## `medsc_holp`

- **Description**: Performs mediator screening using the Ridge-HOLP
  method.
- **Usage**: `medsc_holp(y, x, M, d = NULL, r = 1)`
- **Arguments**:
  - `y`: Outcome vector.
  - `x`: Exposure vector.
  - `M`: Matrix of mediators.
  - `COV.S`: a `data.frame` or `matrix` of covariates (optional).
  - `d`: Desired number of mediators to select (optional). Default value
    is $d= n/\log(n)$.
  - `r`: Ridge-HOLP penalty parameter (default is 1).

## `app_orth`

- **Description**: Fits the Approximate Orthogonalization model to test
  the significance of mediators.
- **Usage**:  
  `app_orth(y, x, chosen_M, COV.S=NULL, k = 1)`  
- **Arguments**:
  - `y`: Outcome vector.
  - `x`: Exposure vector.
  - `chosen_M`: Matrix of mediators selected during screening.
  - `COV.S`: a `data.frame` or `matrix` of covariates (optional).
  - `k`: Scalar used to compute projection directions (default is 1).

## `get_active_med.hima`, `get_active_med.mod`, `get_active_med.hdmt`

- **Description**: These are wrapper functions that combine screening
  and testing to identify active mediators using different
  joint-significance test variants.
- **Usage**:  
  `get_active_med.<variant>(y, x, M, COV.S=NULL, pval.adjust='HDMT', d=NULL, r=1, k=1)`  
- **Variants**:
  - `get_active_med.hima`: Uses the `null_estimation` function from the
    original HIMA repository.
  - `get_active_med.mod`: Uses a modified version of the
    `null_estimation` function.
  - `get_active_med.hdmt`: Uses the implementation in the HDMT package.
  - `get_active_med.fast`: Supports all three methods above and allows
    specifying the `variant` for null estimation.
- **Arguments**:
  - `y`: Outcome vector.
  - `x`: Exposure vector.
  - `M`: Matrix of mediators.
  - `COV.S`: A `data.frame` or `matrix` of covariates (optional).
  - `pval.adjust`: Specifies which method to use for controlling FWER or
    FDR in the joint significance testing. Either `HDMT` (default) or
    `Bonferroni`.
  - `d`: Desired number of mediators to select (optional). Default is
    $d = n / \log(n)$.
  - `r`: Ridge penalty for HOLP (default = 1).
  - `k`: Scalar for computing projection directions for AO (default =
    1).

``` r
# Load the ExampleData
library(hdjmtTest)
data(ExampleData)

# Extract the components
y <- ExampleData$y
x <- ExampleData$x
M <- ExampleData$M
```

# Mediator Screening using Ridge-HOLP

``` r
# Perform Ridge-HOLP screening
chosen_ind <- medsc_holp(y, x, M)

print("Indexes of selected mediators:")
#> [1] "Indexes of selected mediators:"
print(chosen_ind)
#>  [1]    1    2    3    4    5    6    7    8    9   13   14   15   16  122  143
#> [16]  158  225  231  335  508  564  617  621  735  770  841  874  913 1121 1132
#> [31] 1141 1198 1199 1261 1283 1312 1341 1713
```

# Fitting Approximate Orthogonalization (AO)

``` r
# Apply AO on the chosen mediators
chosen_med <- M[, medsc_holp(y, x, M)]
ao_result <- app_orth(y, x, chosen_med)

print("Test statistics for selected mediators:")
#> [1] "Test statistics for selected mediators:"
print(ao_result$ts)
#>  [1] -3.74992060  2.69746840  3.02959821 -3.75784237 -4.00899067  3.12368544
#>  [7] -6.04295950 -5.15851686 -4.34526882 -1.38353473 -0.09319805  1.58176699
#> [13] -0.94060837 -0.63105086 -0.30072448 -2.73986510 -2.96439066 -2.54683352
#> [19]  0.65730656 -0.49086227  1.05970206  0.96275885 -2.41355155  0.08815957
#> [25] -0.29510390 -0.52217751  0.54425397  0.85680171 -1.87447185  1.09829248
#> [31]  0.62027384  0.66017099  1.48312019  0.61257594 -1.09832611  1.57864991
#> [37]  0.16413432  0.65825643

print("P-values for selected mediators:")
#> [1] "P-values for selected mediators:"
print(ao_result$pval)
#>  [1] 1.768906e-04 6.986892e-03 2.448793e-03 1.713848e-04 6.097883e-05
#>  [6] 1.786013e-03 1.513126e-09 2.489137e-07 1.391052e-05 1.665010e-01
#> [11] 9.257462e-01 1.137028e-01 3.469056e-01 5.280073e-01 7.636246e-01
#> [16] 6.146441e-03 3.032830e-03 1.087053e-02 5.109838e-01 6.235239e-01
#> [21] 2.892802e-01 3.356686e-01 1.579789e-02 9.297499e-01 7.679145e-01
#> [26] 6.015467e-01 5.862667e-01 3.915545e-01 6.086542e-02 2.720768e-01
#> [31] 5.350775e-01 5.091441e-01 1.380425e-01 5.401568e-01 2.720621e-01
#> [36] 1.144164e-01 8.696254e-01 5.103734e-01
```

# Identifying Active Mediators

``` r
# Using HDMT (HIMA null_estimation)
active_mediators.hima <- get_active_med.hima(y, x, M)
#> Step 1: Ridge-HOLP Screening   -----  03:25:58 PM
#> Step 2: Approximate Orthogonalization Estimates   -----  03:25:58 PM
#> Step 3: Joint Significance Testing   -----  03:25:59 PM
#> Complete!!   03:25:59 PM
print(active_mediators.hima)
#> [1] 1 2 3 4 5 6 7 8

# Using HDMT (modified null_estimation)
active_mediators.mod <- get_active_med.mod(y, x, M)
#> Step 1: Ridge-HOLP Screening   -----  03:25:59 PM
#> Step 2: Approximate Orthogonalization Estimates   -----  03:25:59 PM
#> Step 3: Joint Significance Testing   -----  03:26:00 PM
#> Complete!!   03:26:00 PM
print(active_mediators.mod)
#> [1] 1 2 3 4 5 6 7 8

# Using HDMT package implementation
active_mediators.hdmt <- get_active_med.hdmt(y, x, M)
#> Step 1: Ridge-HOLP Screening   -----  03:26:00 PM
#> Step 2: Approximate Orthogonalization Estimates   -----  03:26:01 PM
#> Step 3: Joint Significance Testing   -----  03:26:02 PM
#> Complete!!   03:26:02 PM
print(active_mediators.hdmt)
#> [1] 1 2 3 4 5 6 7 8

# Using Bonferroni correction
active_mediators_Bonferroni <- get_active_med.hima(y, x, M, pval.adjust='bonferroni')
#> Step 1: Ridge-HOLP Screening   -----  03:26:02 PM
#> Step 2: Approximate Orthogonalization Estimates   -----  03:26:02 PM
#> Step 3: Joint Significance Testing   -----  03:26:03 PM
#> Complete!!   03:26:03 PM
print(active_mediators_Bonferroni)
#> [1] 1 4 5 7 8
```

# Competing packages

``` r
## Install the HIMA package
# install.packages('HIMA')

suppressMessages(library(HIMA))
#> Warning: package 'HIMA' was built under R version 4.3.3
#> Warning: package 'ncvreg' was built under R version 4.3.3
#> Warning: package 'glmnet' was built under R version 4.3.3
#> Warning in .recacheSubclasses(def@className, def, env): undefined subclass
#> "ndiMatrix" of class "replValueSp"; definition not updated

HIMA::hima_dblasso(x, M, y)
#> Step 1: Sure Independent Screening ...  (3:26:03 PM)
#> Step 2: De-biased Lasso Estimates ...   (3:26:03 PM)
#> Step 3: Joint significance test ...     (3:26:10 PM)
#> Done!     (3:26:10 PM)
#>   Index  alpha_hat   alpha_se   beta_hat   beta_se        IDE      rimp
#> 1     1 -0.4989049 0.06159061 -0.7713661 0.2273450  0.3848383  7.465352
#> 2     3  0.6351005 0.05489418  1.3000239 0.2661458  0.8256458 16.016431
#> 3     4  0.5432349 0.05966641 -0.7974223 0.2396099 -0.4331876  8.403263
#> 4     5  0.7310222 0.04849277 -1.0652445 0.2750990 -0.7787174 15.106081
#> 5     6 -0.7606745 0.04613191  0.9062336 0.2796277 -0.6893488 13.372450
#> 6     7 -0.7180654 0.04946084 -1.4829216 0.2617912  1.0648347 20.656377
#> 7     8 -0.7294041 0.04861566 -1.3413962 0.2703725  0.9784199 18.980045
#>           pmax
#> 1 6.914881e-04
#> 2 1.036232e-06
#> 3 8.747125e-04
#> 4 1.078474e-04
#> 5 1.191708e-03
#> 6 1.474369e-08
#> 7 7.002678e-07
```

Out of the 8 active mediators,
`M[, 1], M[,3], M[,4], M[,5], M[,6], M[,7], M[,8]` were identified as
active by HIMA (specifically HIMA2).

# Reference

- Wang, X., & Leng, C. (2016). High dimensional ordinary least squares
  projection for screening variables. *Journal of the Royal Statistical
  Society Series B: Statistical Methodology*, 78(3), 589–611.
- Battey, H. S., & Reid, N. (2023). On inference in high-dimensional
  regression. *Journal of the Royal Statistical Society Series B:
  Statistical Methodology*, 85(1), 149–175.
