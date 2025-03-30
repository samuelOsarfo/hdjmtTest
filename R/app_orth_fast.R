#' Fast Approximate Orthogonalization
#'
#' Performs approximate orthogonalization using Rcpp for improved performance.
#' This implements the method described by Battey and Reid (2023).
#' @name app_orth_fast
#' @param y Response vector (numeric)
#' @param x Main predictors (numeric matrix)
#' @param chosen_M Selected mediators (numeric matrix)
#' @param COV.S Optional covariates (numeric matrix, optional)
#' @param k Regularization parameter (numeric, default=1)
#' @return List containing:
#' \itemize{
#'   \item \code{ts}: Test statistics
#'   \item \code{pval}: P-values
#' }
#'
#' @examples
#' data(ExampleData) # Load ExampleData1 from the package
#' y <- ExampleData$y # 'y' is a vector of outcomes in ExampleData
#' x <- ExampleData$x # 'x' is a vector of exposures in ExampleData
#' M <- ExampleData$M # 'M' is a matrix of mediators in ExampleData
#'
#' # Perform screening using Ridge-HOLP
#' chosen_ind <- medsc_holp(y, x, M)
#' chosen_M <- M[, chosen_ind]
#'
#' # Apply Approximate Orthogonalization
#' orth_result <- app_orth_fast(y, x, chosen_M)
#' print(orth_result$ts) # Display test statistics
#' print(orth_result$pval) # Display p-values
#'
#' @export
#' @importFrom Rcpp evalCpp
#' @importFrom selectiveInference estimateSigma
#' @useDynLib hdjmtTest, .registration=TRUE
app_orth_fast <- function(y, x, chosen_M, COV.S = NULL, k = 1) {


  # Scale matrices
  chosen_M <- scale(chosen_M)
  XC <- scale(x)

  # Construct MCX matrix
  if (!is.null(COV.S)) {
    COV.S <- scale(COV.S)
    MCX <- cbind(chosen_M, COV.S, XC)
    q2 <- ncol(COV.S)
  } else {
    MCX <- cbind(chosen_M, XC)
    q2 <- 0
  }

  # Compute sigma hat
  sig_hat <- selectiveInference::estimateSigma(MCX, y)$sigmahat

  # Call Rcpp function
  results <- app_orth_rcpp(
    y = as.numeric(y),
    MCX = MCX,
    sig_hat = sig_hat,
    k = k,
    q2 = q2
  )

  return(results)
}
