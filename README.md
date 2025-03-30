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
#>  [1] -3.76599233  2.70902943  3.04258272 -3.77394804 -4.02617274  3.13707320
#>  [7] -6.06885892 -5.18062567 -4.36389214 -1.38946440 -0.09359749  1.58854626
#> [13] -0.94463971 -0.63375547 -0.30201335 -2.75160785 -2.97709570 -2.55774896
#> [19]  0.66012370 -0.49296605  1.06424382  0.96688513 -2.42389575  0.08853741
#> [25] -0.29636868 -0.52441550  0.54658658  0.86047386 -1.88250562  1.10299963
#> [31]  0.62293226  0.66300041  1.48947667  0.61520137 -1.10303341  1.58541583
#> [37]  0.16483778  0.66107764

print("P-values for selected mediators:")
#> [1] "P-values for selected mediators:"
print(ao_result$pval)
#>  [1] 1.658889e-04 6.748035e-03 2.345573e-03 1.606842e-04 5.669205e-05
#>  [6] 1.706435e-03 1.288222e-09 2.211429e-07 1.277686e-05 1.646916e-01
#> [11] 9.254289e-01 1.121629e-01 3.448428e-01 5.262404e-01 7.626419e-01
#> [16] 5.930349e-03 2.909932e-03 1.053521e-02 5.091745e-01 6.220366e-01
#> [21] 2.872183e-01 3.336015e-01 1.535501e-02 9.294496e-01 7.669485e-01
#> [26] 5.999896e-01 5.846628e-01 3.895279e-01 5.976740e-02 2.700273e-01
#> [31] 5.333290e-01 5.073303e-01 1.363619e-01 5.384217e-01 2.700127e-01
#> [36] 1.128719e-01 8.690717e-01 5.085625e-01
```

# Identifying Active Mediators

``` r
# Using HDMT (HIMA null_estimation)
active_mediators.hima <- get_active_med.hima(y, x, M)
#> Step 1: Ridge-HOLP Screening   -----  04:00:07 PM
#> Step 2: Approximate Orthogonalization Estimates   -----  04:00:07 PM
#> Step 3: Joint Significance Testing   -----  04:00:08 PM
#> Complete!!   04:00:09 PM
print(active_mediators.hima)
#> [1] 1 2 3 4 5 6 7 8

# Using HDMT (modified null_estimation)
active_mediators.mod <- get_active_med.mod(y, x, M)
#> Step 1: Ridge-HOLP Screening   -----  04:00:09 PM
#> Step 2: Approximate Orthogonalization Estimates   -----  04:00:09 PM
#> Step 3: Joint Significance Testing   -----  04:00:10 PM
#> Complete!!   04:00:10 PM
print(active_mediators.mod)
#> [1] 1 2 3 4 5 6 7 8

# Using HDMT package implementation
active_mediators.hdmt <- get_active_med.hdmt(y, x, M)
#> Step 1: Ridge-HOLP Screening   -----  04:00:10 PM
#> Step 2: Approximate Orthogonalization Estimates   -----  04:00:10 PM
#> Step 3: Joint Significance Testing   -----  04:00:11 PM
#> Complete!!   04:00:11 PM
print(active_mediators.hdmt)
#> [1] 1 2 3 4 5 6 7 8

# Using Bonferroni correction
active_mediators_Bonferroni <- get_active_med.hima(y, x, M, pval.adjust='bonferroni')
#> Step 1: Ridge-HOLP Screening   -----  04:00:11 PM
#> Step 2: Approximate Orthogonalization Estimates   -----  04:00:11 PM
#> Step 3: Joint Significance Testing   -----  04:00:12 PM
#> Complete!!   04:00:12 PM
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
#> Step 1: Sure Independent Screening ...  (4:00:12 PM)
#> Step 2: De-biased Lasso Estimates ...   (4:00:12 PM)
#> Step 3: Joint significance test ...     (4:00:19 PM)
#> Done!     (4:00:19 PM)
#>   Index  alpha_hat   alpha_se   beta_hat   beta_se        IDE      rimp
#> 1     1 -0.4989049 0.06159061 -0.7688478 0.2455030  0.3835819  7.539413
#> 2     3  0.6351005 0.05489418  1.1656284 0.2893627  0.7402912 14.550634
#> 3     4  0.5432349 0.05966641 -0.7943372 0.2586624 -0.4315117  8.481486
#> 4     5  0.7310222 0.04849277 -1.0409096 0.2995266 -0.7609280 14.956257
#> 5     6 -0.7606745 0.04613191  0.9308149 0.3042549 -0.7080472 13.916870
#> 6     7 -0.7180654 0.04946084 -1.4157387 0.2846038  1.0165929 19.981424
#> 7     8 -0.7294041 0.04861566 -1.4350577 0.2941913  1.0467370 20.573915
#>           pmax
#> 1 1.737825e-03
#> 2 5.619099e-05
#> 3 2.133843e-03
#> 4 5.105058e-04
#> 5 2.218360e-03
#> 6 6.544336e-07
#> 7 1.071809e-06
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
