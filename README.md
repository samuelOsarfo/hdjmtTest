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

The example data provided with the package. The data contains:

- `M`: A 200x2000 matrix of mediators generated using the compound
  symmetry covariance structure to introduce high correlation among the
  mediators. Only the first 8 mediators are active.
- `y`: A vector of length 200 representing outcomes.
- `x`: A vector of length 200 representing exposures.
- `alp_vec`: A parameter vector of length 200 that relates the exposure
  variable to the mediators. The first 8 values are
  `-0.5127127, 0.6597036, 0.6756640, 0.5235137, 0.9305369, -0.9827865, -0.8941141, -0.9230220`
  with the rest being zeros.
- `beta_vec`:A parameter vector of length 200 that relates the mediators
  to the outcome variable. The first 12 values are
  `-0.8033093, 0.9360975, 0.8185305, -0.7951502, -0.8783739, 0.8940459, -0.6911509, -0.8524771, -0.6812097, -0.8285034, -0.5986530, -0.9639383`
  with the rest being zeros.

# How to install package from github

``` r
#library(devtools)

#devtools::install_github('samuelOsarfo/hdjmtTest')
```

If package “qvalue” is not found, please first install “qvalue” package
through Bioconductor:
<https://www.bioconductor.org/packages/release/bioc/html/qvalue.html>

``` r
if (!require("BiocManager", quietly = TRUE))
    install.packages("BiocManager")
#> Warning: package 'BiocManager' was built under R version 4.3.3
#> Bioconductor version '3.18' is out-of-date; the current release version '3.20'
#>   is available with R version '4.4'; see https://bioconductor.org/install
BiocManager::install("qvalue")
#> 'getOption("repos")' replaces Bioconductor standard repositories, see
#> 'help("repositories", package = "BiocManager")' for details.
#> Replacement repositories:
#>     CRAN: https://cran.rstudio.com/
#> Bioconductor version 3.18 (BiocManager 1.30.25), R 4.3.1 (2023-06-16 ucrt)
#> Warning: package(s) not installed when version(s) same as or greater than current; use
#>   `force = TRUE` to re-install: 'qvalue'
#> Installation paths not writeable, unable to update packages
#>   path: C:/Program Files/R/R-4.3.1/library
#>   packages:
#>     boot, cluster, codetools, foreign, KernSmooth, lattice, mgcv, nlme, rpart,
#>     spatial, survival
#> Old packages: 'cli', 'fontawesome', 'rlang'
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
  - `d`: Desired number of mediators to select (optional).
  - `r`: Ridge penalty parameter (default is 1).

## `app_orth`

- **Description**: Fits the Approximate Orthogonalization model to test
  the significance of mediators.
- **Usage**: `app_orth(y, x, chosen_M, k = 1)`
- **Arguments**:
  - `y`: Outcome vector.
  - `x`: Exposure vector.
  - `chosen_M`: Matrix of mediators selected during screening.
  - `k`: Scalar factor used for projection direction (default is 1).

## `get_active_med`

- **Description**: A wrapper function that combines screening and
  testing to identify active mediators.
- **Usage**: `get_active_med(y, x, M)`
- **Arguments**:
  - `y`: Outcome vector.
  - `x`: Exposure vector.
  - `M`: Matrix of mediators.

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
medsc_holp(y, x, M)$chosen_ind_approved  # new implementation following HIMA2
#>  [1]    1    2    3    4    5    6    7    8   13   14   15   16  158  564  841
#> [16]  874 1199 1261 1713
medsc_holp(y, x, M)$chosen_ind  #old implementation 
#>  [1]    1    2    3    4    5    6    7    8   12   13   14   15   16  303  608
#> [16]  991 1199 1435 1713

# Print the indices of selected mediators
#print("Indexes of selected mediators:")
#print(chosen_ind)
```

# Fitting Approximate Orthogonalization (AO)

With the selected mediators, we can use `app_orth` to fit the
Approximate Orthogonalization model. This model computes test statistics
and p-values for the selected mediators.

``` r
# Apply AO on the chosen mediators
chosen_med <- M[, medsc_holp(y, x, M)$chosen_ind_approved]
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
  Bonferroni correction is applied to the p-values.

- **Output**: Indexes of the mediators that are deemed significant after
  the joint-significance test.

``` r
# Identify active mediators
active_mediators <- get_active_med(y, x, M)

# Display indices of active mediators and the corresponding matrix
print("Indexes of active mediators:")
#> [1] "Indexes of active mediators:"
print(active_mediators)
#> $which_sig_benferroni
#> [1] 1 4 5 7 8
#> 
#> $which_with_fdr_in_HDMT_package_approved
#> [1] 1 2 3 4 5 6 7 8 9
```

# Competing packages

We can perform this test using `HIMA` package as well.

``` r
## Install the HIMA package
install.packages('HIMA')
#> Installing package into 'C:/Users/fsosa/AppData/Local/Temp/RtmpKQDP3I/temp_libpath6e6c43336c93'
#> (as 'lib' is unspecified)
#> package 'HIMA' successfully unpacked and MD5 sums checked
#> 
#> The downloaded binary packages are in
#>  C:\Users\fsosa\AppData\Local\Temp\RtmpyuZzgS\downloaded_packages

library(HIMA)
#> Warning: package 'HIMA' was built under R version 4.3.3
#> Loading required package: ncvreg
#> Warning: package 'ncvreg' was built under R version 4.3.3
#> Loading required package: glmnet
#> Warning: package 'glmnet' was built under R version 4.3.3
#> Loading required package: Matrix
#> Loaded glmnet 4.1-8
#> Warning in .recacheSubclasses(def@className, def, env): undefined subclass
#> "ndiMatrix" of class "replValueSp"; definition not updated
#> ************************************************************************************************************************
#> HIMA version 2.3.0
#> To access full functionality of HIMA, please make sure this version is current.
#> 
#> Citation:
#>   1. Zhang H, Zheng Y, Zhang Z, Gao T, Joyce B, Yoon G, Zhang W, Schwartz J,
#>      Just A, Colicino E, Vokonas P, Zhao L, Lv J, Baccarelli A, Hou L, Liu L.
#>      Estimating and Testing High-dimensional Mediation Effects in Epigenetic Studies.
#>      Bioinformatics. 2016.
#>      PMID: 27357171; PMCID: PMC5048064.
#> 
#>   2. Zhang H, Zheng Y, Hou L, Zheng C, Liu L.
#>      Mediation Analysis for Survival Data with High-Dimensional Mediators.
#>      Bioinformatics. 2021.
#>      PMID: 34343267; PMCID: PMC8570823.
#> 
#>   3. Zhang H, Chen J, Feng Y, Wang C, Li H, Liu L.
#>      Mediation effect selection in high-dimensional and compositional microbiome data.
#>      Stat Med. 2021.
#>      PMID: 33205470; PMCID: PMC7855955.
#> 
#>   4. Perera C, Zhang H, Zheng Y, Hou L, Qu A, Zheng C, Xie K, Liu L.
#>      HIMA2: high-dimensional mediation analysis and its application in epigenome-wide DNA methylation data.
#>      BMC Bioinformatics. 2022.
#>      PMID: 35879655; PMCID: PMC9310002.
#> 
#>   5. Zhang H, Hong X, Zheng Y, Hou L, Zheng C, Wang X, Liu L.
#>      High-Dimensional Quantile Mediation Analysis with Application to a Birth Cohort Study of Mother-Newborn Pairs.
#>      Bioinformatics. 2024.
#>      PMID: 38290773; PMCID: PMC10873903.
#> 
#>   6. Bai X, Zheng Y, Hou L, Zheng C, Liu L, Zhang H.
#>      An Efficient Testing Procedure for High-dimensional Mediators with FDR Control.
#>      Statistics in Biosciences. 2024.
#> ************************************************************************************************************************
```

The `dblassoHIMA` function in the HIMA package combines the sure
independence screening method and de-biased lasso to identify the active
mediators. Consider the rcode below for identifying the active mediators
in the simulated data.

``` r

HIMA::dblassoHIMA(x,M,y)
#> Step 1: Sure Independent Screening ...  (10:58:59 PM)
#> Step 2: De-biased Lasso Estimates ...   (10:58:59 PM)
#> Step 3: Joint significance test ...     (10:59:06 PM)
#> Done!     (10:59:06 PM)
#>   Index  alpha_hat   alpha_se   beta_hat   beta_se        IDE      rimp
#> 1     1 -0.4989049 0.06159061 -0.7553381 0.2351410  0.3768419  7.447952
#> 2     3  0.6351005 0.05489418  1.2139274 0.2752724  0.7709659 15.237469
#> 3     4  0.5432349 0.05966641 -0.7969289 0.2478265 -0.4329196  8.556278
#> 4     5  0.7310222 0.04849277 -1.0000087 0.2845326 -0.7310286 14.448143
#> 5     6 -0.7606745 0.04613191  0.9201336 0.2892166 -0.6999222 13.833353
#> 6     7 -0.7180654 0.04946084 -1.4384565 0.2707685  1.0329058 20.414485
#> 7     8 -0.7294041 0.04861566 -1.3916668 0.2796440  1.0150875 20.062321
#>           pmax
#> 1 1.316873e-03
#> 2 1.034121e-05
#> 3 1.301393e-03
#> 4 4.404721e-04
#> 5 1.465304e-03
#> 6 1.081337e-07
#> 7 6.472224e-07
```

Out of the 8 active mediators,
`M[, 1], M[,3], M[,4], M[,5], M[,6], M[,7], M[,8]` were identified as
active my HIMA (specifically HIMA2).

# Conclusion

The `hdjmtTest` package simplifies high-dimensional mediation analysis,
providing tools for selecting active mediators under high correlation
scenarios.

Refer to the function documentation for more details.
