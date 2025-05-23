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
  and testing to identify active mediators using different variants of the
  joint-significance test.
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
#>  [1] -3.75826498  2.70347084  3.03633972 -3.76620437 -4.01791153  3.13063631
#>  [7] -6.05640638 -5.16999566 -4.35493798 -1.38661339 -0.09340543  1.58528676
#> [13] -0.94270143 -0.63245508 -0.30139365 -2.74596189 -2.97098707 -2.55250077
#> [19]  0.65876921 -0.49195454  1.06206013  0.96490120 -2.41892222  0.08835574
#> [25] -0.29576057 -0.52333947  0.54546505  0.85870828 -1.87864295  1.10073641
#> [31]  0.62165408  0.66164001  1.48642045  0.61393906 -1.10077012  1.58216275
#> [37]  0.16449956  0.65972119

print("P-values for selected mediators:")
#> [1] "P-values for selected mediators:"
print(ao_result$pval)
#>  [1] 1.710956e-04 6.861947e-03 2.394694e-03 1.657481e-04 5.871623e-05
#>  [6] 1.744280e-03 1.391961e-09 2.340994e-07 1.331046e-05 1.655597e-01
#> [11] 9.255815e-01 1.129012e-01 3.458336e-01 5.270895e-01 7.631143e-01
#> [16] 6.033377e-03 2.968443e-03 1.069527e-02 5.100440e-01 6.227515e-01
#> [21] 2.882084e-01 3.345943e-01 1.556657e-02 9.295939e-01 7.674129e-01
#> [26] 6.007380e-01 5.854337e-01 3.905015e-01 6.029326e-02 2.710114e-01
#> [31] 5.341694e-01 5.082020e-01 1.371679e-01 5.392556e-01 2.709967e-01
#> [36] 1.136124e-01 8.693379e-01 5.094328e-01
```

# Identifying Active Mediators

``` r
# Using HDMT (HIMA null_estimation)
active_mediators.hima <- get_active_med.hima(y, x, M)
#> Step 1: Ridge-HOLP Screening   -----  04:06:23 PM
#> Step 2: Approximate Orthogonalization Estimates   -----  04:06:23 PM
#> Step 3: Joint Significance Testing   -----  04:06:24 PM
#> Complete!!   04:06:25 PM
print(active_mediators.hima)
#> [1] 1 2 3 4 5 6 7 8

# Using HDMT (modified null_estimation)
active_mediators.mod <- get_active_med.mod(y, x, M)
#> Step 1: Ridge-HOLP Screening   -----  04:06:25 PM
#> Step 2: Approximate Orthogonalization Estimates   -----  04:06:25 PM
#> Step 3: Joint Significance Testing   -----  04:06:26 PM
#> Complete!!   04:06:26 PM
print(active_mediators.mod)
#> [1] 1 2 3 4 5 6 7 8

# Using HDMT package implementation
active_mediators.hdmt <- get_active_med.hdmt(y, x, M)
#> Step 1: Ridge-HOLP Screening   -----  04:06:26 PM
#> Step 2: Approximate Orthogonalization Estimates   -----  04:06:26 PM
#> Step 3: Joint Significance Testing   -----  04:06:27 PM
#> Complete!!   04:06:27 PM
print(active_mediators.hdmt)
#> [1] 1 2 3 4 5 6 7 8

# Using Bonferroni correction
active_mediators_Bonferroni <- get_active_med.hima(y, x, M, pval.adjust='bonferroni')
#> Step 1: Ridge-HOLP Screening   -----  04:06:27 PM
#> Step 2: Approximate Orthogonalization Estimates   -----  04:06:27 PM
#> Step 3: Joint Significance Testing   -----  04:06:28 PM
#> Complete!!   04:06:28 PM
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
#> Step 1: Sure Independent Screening ...  (4:06:28 PM)
#> Step 2: De-biased Lasso Estimates ...   (4:06:28 PM)
#> Step 3: Joint significance test ...     (4:06:35 PM)
#> Done!     (4:06:35 PM)
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
