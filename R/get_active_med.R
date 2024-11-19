#' Get Active Mediators
#'
#' @param y  vector of outcome
#' @param x vector of exposure
#' @param M a \code{data.frame}\\code{matrix} of mediators
#' @param COV.S a \code{data.frame}\\code{matrix} of covariates
#' @param pval.adjust specifies which method to use for controlling FWER in the joint significance testing. Either \code{'HDMT'} (default) or \code{'bonferroni'}
#' @param d the number of screened mediators
#' @param r  scalar multiplied by the identity matrix for the Ridge-HOLP. Default value is 1
#' @param k a scalar multiplied by identity matrix to compute projection direction for AO. Default value 1.
#'
#'
#' @return A vector of the indexes for active mediators.
#' @export
#'
#' @examples
#' data(ExampleData) # Load ExampleData from the package
#' y <- ExampleData$y # 'y' is a vector of outcomes in ExampleData
#' x <- ExampleData$x # assuming 'x' is a vector of exposures in ExampleData
#' M <- ExampleData$M # assuming 'M' is a matrix of mediators in ExampleData
#'
#' # Get active mediators
#' active_mediators_index <- get_active_med(y, x, M)
#' print(active_mediators_index) # Print the indexes of active mediators
#'
get_active_med <- function(y, x, M, COV.S=NULL, pval.adjust='HDMT', d=NULL, r=1,  k=1){

  #screen mediators
  message("Step 1: Ridge-HOLP Screening   -----  ", format(Sys.time(), "%I:%M %p"))
     chosen_ind <- medsc_holp(y, x, M,  COV.S, d ,r)

  #get p-values and test statistics for \beta_j's in outcome model
     message("Step 2: Approximate Orthogonalization Estimates   -----  ", format(Sys.time(), "%I:%M %p"))
     ao_obj <- app_orth(y, x, M[, chosen_ind], COV.S, k)


  #get p-values and test statistics for \alpha_j's in mediator model
     alp_all <- comp_alpha(x, M)

  #get index for active mediators in chosen mediators
     message("Step 3: Joint Significance Testing   -----  ", format(Sys.time(), "%I:%M %p"))
     active_index <- js_test(chosen_ind, alp_all$pval[chosen_ind], ao_obj$pval, method=pval.adjust)

  #results
     message("Complete!!", format(Sys.time(), "%I:%M %p"))

    return(active_index)
}
