#' Get Active Mediators :: Null Estimation Function From the HDMT package
#'
#' @param y  vector of outcomes
#' @param x vector of exposures
#' @param M a \code{data.frame} or \code{matrix} of mediators
#' @param COV.S a \code{data.frame} or \code{matrix} of covariates
#' @param pval.adjust specifies which method to use for controlling FWER/FDR in the joint significance testing. Either \code{'HDMT'} (default) or \code{'Bonferroni'}
#' @param d the number of screened mediators. Default value is \eqn{d =  n/\log(n)}.
#' @param r  a penalty parameter for the Ridge-HOLP. Default value is `1`
#' @param k  a scalar for computing projection directions for AO. Default value is `1`.
#'
#' @return A vector of the indexes for active mediators.
#' @export
#'
#' @examples
#' data(ExampleData) # Load ExampleData from the package
#' y <- ExampleData$y # 'y' is a vector of outcomes in ExampleData
#' x <- ExampleData$x #  'x' is a vector of exposures in ExampleData
#' M <- ExampleData$M #  'M' is a matrix of mediators in ExampleData
#'
#' # Get active mediators
#' #using HDMT to control FDR (default)
#' active_mediators_index <- get_active_med.hdmt(y, x, M)
#'
#' # Print the indexes of active mediators
#' print(active_mediators_index)
#'
#'
#'#using Bonferroni correction to control FWER
#' active_mediators_index <- get_active_med.hdmt(y, x, M, pval.adjust='bonferroni')
#'
#' #Print the indexes of active mediators
#'print(active_mediators_index)
get_active_med.hdmt <- function(y, x, M, COV.S=NULL, pval.adjust='HDMT', d=NULL, r=1,  k=1){

  #screen mediators
  message("Step 1: Ridge-HOLP Screening   -----  ", format(Sys.time(), "%I:%M:%S %p"))
     chosen_ind <- medsc_holp(y, x, M,  COV.S, d ,r)

  #get p-values and test statistics for \beta_j's in outcome model
     message("Step 2: Approximate Orthogonalization Estimates   -----  ", format(Sys.time(), "%I:%M:%S %p"))
     ao_obj <- app_orth(y, x, M[, chosen_ind], COV.S, k)


  #get p-values and test statistics for \alpha_j's in mediator model
     alp_all <- comp_alpha(x, M[,chosen_ind], COV.S)


  #get index for active mediators in chosen mediators
     message("Step 3: Joint Significance Testing   -----  ", format(Sys.time(), "%I:%M:%S %p"))
     active_index <- js_test_hdmt(chosen_ind, alp_all$pval, ao_obj$pval, method=pval.adjust)

  #results
     message("Complete!!   ", format(Sys.time(), "%I:%M:%S %p"))

    return(active_index)
}
