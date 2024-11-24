Documentation For hdjmtTest
================

<!-- README.md is generated from README.Rmd. Please edit that file -->

# Introduction

The `hdjmtTest` package performs high-dimensional mediator screening and
testing in the presence of high correlation among mediators using
Ridge-HOLP and Approximate Orthogonalization methods. This vignette
demonstrates how to use the primary functions with a simulated dataset,
`ExampleData`.

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
    `-0.8033093, 0.9360975, 0.8185305, -0.7951502, -0.8783739, 0.8940459, -0.6911509, -0.8524771, -0.6812097, -0.8285034,           -0.5986530, -0.9639383`
    with the rest being zeros.
- **Active Mediators**: given the non-zero values in `alp_vec` and
  `beta_vec`, only the first 8 mediators are truly active

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
    is $d= 0.5\cdot n/\log(n)$.
  - `r`: Ridge-HOLP penalty parameter (default is 1).

## `app_orth`

- **Description**: Fits the Approximate Orthogonalization model to test
  the significance of mediators.
- **Usage**: `app_orth(y, x, chosen_M, COV.S=NULL, k = 1)`
- **Arguments**:
  - `y`: Outcome vector.
  - `x`: Exposure vector.
  - `chosen_M`: Matrix of mediators selected during screening.
  - `COV.S`: a `data.frame` or `matrix` of covariates (optional).
  - `k`: Scalar used to compute projection directions (default is 1).

## `get_active_med`

- **Description**: A wrapper function that combines screening and
  testing to identify active mediators.
- **Usage**:
  `get_active_med(y, x, M, COV.S=NULL,pval.adjust='HDMT',d=NULL, r=1, k=1)`
- **Arguments**:
  - `y`: Outcome vector.
  - `x`: Exposure vector.
  - `M`: Matrix of mediators.
  - `COV.S`: A `data.frame` or `matrix` of covariates (optional).
  - `pval.adjust`: Specifies which method to use for controlling FWER in
    the joint significance testing. Either `HDMT` (default) or
    `Bonferroni`.
  - `d`: Desired number of mediators to select (optional). Default value
    is $d= 0.5\cdot n/\log(n)$.
  - `r`: Penalty parameter for the Ridge-HOLP. Default value is `1`
  - `k` : Scalar used to compute projection directions for AO. Default
    value is `1`.

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

The `medsc_holp` function performs Ridge-HOLP screening to select
potential active mediators based on their association with the outcome
and exposure.

``` r
# Perform Ridge-HOLP screening
chosen_ind <- medsc_holp(y, x, M)


print("Indexes of selected mediators:")
#> [1] "Indexes of selected mediators:"
print(chosen_ind)
#>  [1]    1    2    3    4    5    6    7    8   13   14   15   16  158  564  841
#> [16]  874 1199 1261 1713
```

# Fitting Approximate Orthogonalization (AO)

With the selected mediators, we can use `app_orth` to fit the
Approximate Orthogonalization model. This model computes test statistics
and p-values for the selected mediators.

``` r
# Apply AO on the chosen mediators
chosen_med <- M[, medsc_holp(y, x, M)]
ao_result <- app_orth(y, x, chosen_med)


# Print the test statistics and p-values
print("Test statistics for selected mediators:")
#> [1] "Test statistics for selected mediators:"
print(ao_result$ts)
#>  [1] -4.1440145  2.4370064  2.9433927 -4.7247547 -5.0984107  2.7990726
#>  [7] -6.1042222 -5.1643208 -2.3633393 -0.2800171  0.2995497 -1.3588344
#> [13] -3.0376229  0.4060232 -1.1739012  0.4506148  1.6929801 -0.4317281
#> [19]  0.0411413

print("P-values for selected mediators:")
#> [1] "P-values for selected mediators:"
print(ao_result$pval)
#>  [1] 3.412782e-05 1.480942e-02 3.246363e-03 2.303933e-06 3.425170e-07
#>  [6] 5.124961e-03 1.033025e-09 2.413136e-07 1.811108e-02 7.794644e-01
#> [11] 7.645206e-01 1.741991e-01 2.384521e-03 6.847255e-01 2.404346e-01
#> [16] 6.522672e-01 9.045926e-02 6.659390e-01 9.671833e-01
```

# Identifying Active Mediators

The `get_active_med` function integrates screening and testing,
identifying active mediators and their corresponding indexes. The steps
involved are as follows:

### 1. Mediator Screening with Ridge-HOLP

The function starts by using the Ridge-HOLP method to reduce the
dimensionality of the mediator matrix (`M`). This step ensures that only
potential mediators with strong associations with both the exposure
(`x`) and the outcome (`y`) are considered for further testing.

### 2. Fitting the AO Model

After initial screening, the function fits the Approximate
Orthogonalization (AO) model to the selected mediators. This fitting
produces test statistics and p-values that help determine the
significance of each mediator.

### 3. Joint-Significance Test

- A joint-significance test is performed to identify active mediators
  using the p-values obtained from both the AO model and the mediator
  model.

- To control for the family-wise error rate due to multiple testing,
  Bonferroni/HDMT correction is applied to the p-values.

- **Output**: Indexes of the mediators that are deemed significant after
  the joint-significance test.

``` r
# Identify active mediators

#Using HDMT in joint significance testing (Default)
active_mediators_HDMT <- get_active_med(y, x, M) 
#> Step 1: Ridge-HOLP Screening   -----  08:45:23 PM
#> Step 2: Approximate Orthogonalization Estimates   -----  08:45:23 PM
#> Step 3: Joint Significance Testing   -----  08:45:24 PM
#> Complete!!   08:45:25 PM

#Indexes of active mediators identified using HDMT:"
print(active_mediators_HDMT)
#> [1]  1  2  3  4  5  6  7  8 13



#Using Bonferroni in joint significance testing
active_mediators_Bonferroni <- get_active_med(y, x, M, pval.adjust='bonferroni') 
#> Step 1: Ridge-HOLP Screening   -----  08:45:25 PM
#> Step 2: Approximate Orthogonalization Estimates   -----  08:45:26 PM
#> Step 3: Joint Significance Testing   -----  08:45:27 PM
#> Complete!!   08:45:27 PM

#Indexes of active mediators identified using Bonferroni:"
print(active_mediators_Bonferroni)
#> [1] 1 4 5 7 8
```

# Competing packages

We can perform this test using `HIMA` package as well.

``` r
## Install the HIMA package
#install.packages('HIMA')

suppressMessages(library(HIMA))
#> Warning: package 'HIMA' was built under R version 4.3.3
#> Warning: package 'ncvreg' was built under R version 4.3.3
#> Warning: package 'glmnet' was built under R version 4.3.3
#> Warning in .recacheSubclasses(def@className, def, env): undefined subclass
#> "ndiMatrix" of class "replValueSp"; definition not updated
```

The `dblassoHIMA` function in the HIMA package combines the sure
independence screening method and de-biased lasso to identify the active
mediators. Consider the rcode below for identifying the active mediators
in the simulated data.

``` r

HIMA::dblassoHIMA(x,M,y)
#> Step 1: Sure Independent Screening ...  (8:45:27 PM)
#> Step 2: De-biased Lasso Estimates ...   (8:45:28 PM)
#> Step 3: Joint significance test ...     (8:45:34 PM)
#> Done!     (8:45:34 PM)
#>   Index  alpha_hat   alpha_se   beta_hat   beta_se        IDE      rimp
#> 1     1 -0.4989049 0.06159061 -0.7624210 0.2317676  0.3803756  7.463652
#> 2     3  0.6351005 0.05489418  1.2455445 0.2713233  0.7910459 15.521741
#> 3     4  0.5432349 0.05966641 -0.7964028 0.2442711 -0.4326338  8.489052
#> 4     5  0.7310222 0.04849277 -1.0260188 0.2804506 -0.7500425 14.717181
#> 5     6 -0.7606745 0.04613191  0.9155677 0.2850675 -0.6964490 13.665580
#> 6     7 -0.7180654 0.04946084 -1.4551693 0.2668840  1.0449067 20.502945
#> 7     8 -0.7294041 0.04861566 -1.3722434 0.2756322  1.0009200 19.639848
#>           pmax
#> 1 1.003326e-03
#> 2 4.419123e-06
#> 3 1.112854e-03
#> 4 2.537300e-04
#> 5 1.319254e-03
#> 6 4.968293e-08
#> 7 6.406853e-07
```

Out of the 8 active mediators,
`M[, 1], M[,3], M[,4], M[,5], M[,6], M[,7], M[,8]` were identified as
active my HIMA (specifically HIMA2).

# Conclusion

The `hdjmtTest` package simplifies high-dimensional mediation analysis,
providing tools for selecting active mediators under high correlation
scenarios. Refer to the function documentations for more details.

# Reference

- Wang, X., & Leng, C. (2016). High dimensional ordinary least squares
  projection for screening variables. Journal of the Royal Statistical
  Society Series B: Statistical Methodology, 78(3), 589-611.

- Battey, H. S., & Reid, N. (2023). On inference in high-dimensional
  regression. Journal of the Royal Statistical Society Series B:
  Statistical Methodology, 85(1), 149-175.
