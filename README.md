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

#devtools::install_github('samuelOsarfo/hdjmtTest',force=TRUE)
```

If package “qvalue” is not found, please first install “qvalue” package
through Bioconductor:
<https://www.bioconductor.org/packages/release/bioc/html/qvalue.html>

``` r
#if (!require("BiocManager", quietly = TRUE))
#    install.packages("BiocManager")
#BiocManager::install("qvalue")

#devtools::install_github('samuelOsarfo/hdjmtTest',force=TRUE)
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
#>  [1] -3.76107863  2.70549481  3.03861289 -3.76902396 -4.02091957  3.13298008
#>  [7] -6.06094055 -5.17386622 -4.35819833 -1.38765149 -0.09347536  1.58647359
#> [13] -0.94340719 -0.63292858 -0.30161929 -2.74801767 -2.97321132 -2.55441172
#> [19]  0.65926240 -0.49232285  1.06285525  0.96562358 -2.42073316  0.08842189
#> [25] -0.29598199 -0.52373127  0.54587342  0.85935116 -1.88004941  1.10156049
#> [31]  0.62211949  0.66213536  1.48753327  0.61439869 -1.10159422  1.58334725
#> [37]  0.16462271  0.66021510

print("P-values for selected mediators:")
#> [1] "P-values for selected mediators:"
print(ao_result$pval)
#>  [1] 1.691822e-04 6.820271e-03 2.376700e-03 1.638872e-04 5.797140e-05
#>  [6] 1.730412e-03 1.353279e-09 2.292990e-07 1.311375e-05 1.652432e-01
#> [11] 9.255259e-01 1.126319e-01 3.454727e-01 5.267803e-01 7.629423e-01
#> [16] 5.995678e-03 2.947015e-03 1.063674e-02 5.097273e-01 6.224911e-01
#> [21] 2.878476e-01 3.342326e-01 1.548924e-02 9.295414e-01 7.672438e-01
#> [26] 6.004655e-01 5.851530e-01 3.901468e-01 6.010134e-02 2.706528e-01
#> [31] 5.338633e-01 5.078845e-01 1.368740e-01 5.389519e-01 2.706381e-01
#> [36] 1.133423e-01 8.692410e-01 5.091158e-01
```

# Identifying Active Mediators

``` r
# Using HDMT (HIMA null_estimation)
active_mediators.hima <- get_active_med.hima(y, x, M)
#> Step 1: Ridge-HOLP Screening   -----  03:28:44 PM
#> Step 2: Approximate Orthogonalization Estimates   -----  03:28:44 PM
#> Step 3: Joint Significance Testing   -----  03:28:45 PM
#> Complete!!   03:28:48 PM
print(active_mediators.hima)
#> [1] 1 2 3 4 5 6 7 8

# Using HDMT (modified null_estimation)
active_mediators.mod <- get_active_med.mod(y, x, M)
#> Step 1: Ridge-HOLP Screening   -----  03:28:48 PM
#> Step 2: Approximate Orthogonalization Estimates   -----  03:28:48 PM
#> Step 3: Joint Significance Testing   -----  03:28:49 PM
#> Complete!!   03:28:49 PM
print(active_mediators.mod)
#> [1] 1 2 3 4 5 6 7 8

# Using HDMT package implementation
active_mediators.hdmt <- get_active_med.hdmt(y, x, M)
#> Step 1: Ridge-HOLP Screening   -----  03:28:49 PM
#> Step 2: Approximate Orthogonalization Estimates   -----  03:28:50 PM
#> Step 3: Joint Significance Testing   -----  03:28:51 PM
#> Complete!!   03:28:51 PM
print(active_mediators.hdmt)
#> [1] 1 2 3 4 5 6 7 8

# Using Bonferroni correction
active_mediators_Bonferroni <- get_active_med.hima(y, x, M, pval.adjust='bonferroni')
#> Step 1: Ridge-HOLP Screening   -----  03:28:51 PM
#> Step 2: Approximate Orthogonalization Estimates   -----  03:28:51 PM
#> Step 3: Joint Significance Testing   -----  03:28:52 PM
#> Complete!!   03:28:52 PM
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
#> Step 1: Sure Independent Screening ...  (3:28:53 PM)
#> Step 2: De-biased Lasso Estimates ...   (3:28:53 PM)
#> Step 3: Joint significance test ...     (3:29:00 PM)
#> Done!     (3:29:00 PM)
#>   Index  alpha_hat   alpha_se   beta_hat   beta_se        IDE      rimp
#> 1     1 -0.4989049 0.06159061 -0.7674050 0.2296871  0.3828621  7.467072
#> 2     3  0.6351005 0.05489418  1.2737283 0.2688877  0.8089454 15.777100
#> 3     4  0.5432349 0.05966641 -0.7968484 0.2420784 -0.4328758  8.442504
#> 4     5  0.7310222 0.04849277 -1.0458150 0.2779331 -0.7645140 14.910540
#> 5     6 -0.7606745 0.04613191  0.9114397 0.2825085 -0.6933089 13.521807
#> 6     7 -0.7180654 0.04946084 -1.4696024 0.2644882  1.0552706 20.581252
#> 7     8 -0.7294041 0.04861566 -1.3566723 0.2731579  0.9895624 19.299725
#>           pmax
#> 1 8.345056e-04
#> 2 2.168760e-06
#> 3 9.958531e-04
#> 4 1.680011e-04
#> 5 1.254289e-03
#> 6 2.753950e-08
#> 7 6.812877e-07
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
