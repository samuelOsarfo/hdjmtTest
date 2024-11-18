#' Get Active Mediators
#'
#' @param y  vector of outcome
#' @param x vector of exposure
#' @param M a \code{data.frame}\\code{matrix} of mediators
#' @param COV.S a \code{data.frame}\\code{matrix} of covariates
#' @param d the number of screened mediators
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
#' active_index <- get_active_med(y, x, M)
#' print(active_index) # Print the indexes of active mediators
#'
get_active_med <- function(y, x, M, COV.S=NULL, d=NULL){

  #screen mediators
     chosen_ind <- medsc_holp(y, x, M,  COV.S, d ,r = 1)$chosen_ind_approved

  #get p-values and test statistics for \beta_j's in outcome model
     ao_obj <- app_orth(y, x, M[, chosen_ind], COV.S)

  #get p-values and test statistics for \alpha_j's in mediator model
     alp_all <- comp_alpha(x, M)

  #get index for active mediators in chosen mediators
     active_index <- js_test(chosen_ind, alp_all$pval[chosen_ind], ao_obj$pval)

  #results
    return(active_index)
}
