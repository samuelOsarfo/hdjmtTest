#' Screen mediators by Ridge-HOLP method
#'
#' @param y  n-dimensional vector of continuous outcomes.
#' @param x  n-dimensional vector of exposures
#' @param M  n by p matrix of mediators
#' @param COV.S a \code{data.frame} or \code{matrix} of covariates
#' @param d  desired number of chosen mediators. Default value is \eqn{d = 0.5 \cdot n/\log(n)}.
#' @param r  scalar for the penalty parameter. Default value is 1.
#'
#' @return A vector for indexes of  the selected mediators.
#' @export
#'
#' @examples
#' data(ExampleData) # Load ExampleData from the data from package
#' y <- ExampleData$y #   'y' is a vector of continuous outcomes in ExampleData
#' x <- ExampleData$x #   'x' is a vector of exposures in ExampleData
#' M <- ExampleData$M #   'M' is a matrix  of mediators in ExampleData
#'
#' chosen_ind <- medsc_holp(y, x, M)
#' print(chosen_ind) # print indexes for the selected mediators
#'
medsc_holp <- function(y, x, M, COV.S=NULL, d = NULL, r = 1){

  M <- scale(M)
  XC <- scale(x)
  n <- length(y);
  p <- ncol(M)
  q <- 0
  MCX <- cbind(M, XC)



  if(!is.null(COV.S)){
    COV.S <- scale(COV.S)
    MCX  <- cbind(M, COV.S, XC)
    XC <- cbind(XC, COV.S)
    q <- ncol(COV.S)
  }



  if(is.null(d)){
    d <- ceiling(0.5* n/log(n))
  }



  # if d > p select all mediators
  num_chosen <- min(p, d)


  # Compute the Ridge-HOLP estimates for beta
  holp_est_beta <- t(MCX) %*% solve(tcrossprod(MCX) + diag(n) * r) %*% y # must be (p + 1) dimensional


  #compute the estimates for alpha
    alpha_est <- matrix(NA, p, 1)
    for(i in 1:p){
      est <- stats::lsfit(XC, M[, i], intercept = FALSE) #removes the intercept  #rem to add stats:: if approved by Dr. Yi
      alpha_est[i, 1] <- est$coefficients[1]

          }


    sort_obj <- sort(abs(holp_est_beta[-c((p + 1):(p + q + 1)),1]* alpha_est[,1]), index.return = T, decreasing = T)
    chosen_ind <- sort(sort_obj$ix[c(1:num_chosen)])




  return(chosen_ind)
}
