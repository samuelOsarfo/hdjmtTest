#' Fit Approximate Orthogonalization by Battey and Reid (2023)
#'
#' @param y  An n-dimensional vector of outcomes
#' @param x  An n-dimensional vector of exposure
#' @param chosen_M An n by q matrix of mediators chosen by some screening method
#' @param COV.S a \code{data.frame}\\code{matrix} of covariates
#' @param k A scalar multiplied by identity matrix to compute projection direction.
#'    k = 1 or 1/n were used in the original paper by Battey and Reid (2023). Default value is 1.
#'
#' @return A a vector of test statistics (ts) each of which asymptotically follow std. normal under the null and
#'    a vector of the corresponding p-values
#' @export

#' @examples
#' data(ExampleData) # Load ExampleData1 from the package
#' y <- ExampleData$y # assuming 'y' is a vector/column in ExampleData
#' x <- ExampleData$x # assuming 'x' is a vector/column in ExampleData
#' M <- ExampleData$M # assuming 'M' is a matrix of mediators in ExampleData
#' # Perform screening using Ridge-HOLP
#' chosen_ind <- medsc_holp(y, x, M)
#' chosen_M <- M[, chosen_ind]
#'
#' # Apply Approximate Orthogonalization
#' orth_result <- app_orth(y, x, chosen_M)
#' print(orth_result$ts) # Display test statistics
#' print(orth_result$pval) # Display p-values
#'
app_orth <- function(y, x, chosen_M, COV.S=NULL, k = 1){

  chosen_M <- scale(chosen_M)
  n <- length(y)
  XC <- scale(x)
  MCX <- cbind(chosen_M, XC)
  q1 <- ncol(MCX)
  q2 <- 0
  ts <- pval <- rep(NA, q1-1)




  if(!is.null(COV.S)){
    COV.S <- scale(COV.S)
    MCX  <- cbind(chosen_M, COV.S, XC)
    q1 <- ncol(MCX)
    q2 <- ncol(COV.S)
    ts <- pval <- rep(NA, (q1-(q2 + 1)))
  }

  sig_hat <- selectiveInference::estimateSigma(MCX, y)$sigmahat # compute the estimates of standard deviation of random errors


  # Compute the projection vector as described in Battey and Reid (2023)
  # and the resulting test statistics and p-values for each selected mediator
  for(j in 1:(q1-(q2 + 1))){
    # the proposed projection direction in Battey and Reid (2023)
    proj_vec <- solve(k*diag(n) + tcrossprod(MCX[,-j]))%*%MCX[,j] # dimension : n by 1
    # compute the test statistic which asymptotically follows a standard normal dist'n
    ts[j] <- sum(proj_vec*y)/(sig_hat*sqrt(sum(proj_vec**2)))
    # compute the p-value
    pval[j] <- 2*stats::pnorm(abs(ts[j]), lower.tail = F)
  }

  return(list(ts = ts, pval = pval))
}
