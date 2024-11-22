###########################################################################################
#### 16. function to compute the p-values to test H_{0}: \alpha_{j}=0 based on the OLS ####
###########################################################################################

# Inputs
# x : n-dimensional vectors of exposure
# chosen_M : n by q matrix of mediators chosen by some screening method

# Output
# ts: a vector of test statistics each of which asymptotically follow std. normal under the null
# pval: a vector of p-values to test \alpha_{j}=0


comp_alpha <- function(x, chosen_M){

  x <- scale(x)
  chosen_M <- scale(chosen_M)
  q <- ncol(chosen_M)
  ts <- pval <- rep(NA, q)

  for(j in 1:q){
    res <- stats::coef(summary(stats::lm(chosen_M[,j]~x)))
    ts[j] <- res[2,3] ; pval[j] <- res[2,4]
  }

  return(list(ts = ts, pval = pval))
}

#####################################################
####### Code from HDMT for js test ############
#####################################################
null_estimation <- function(input_pvalues, eps = 10^-8){
  ## updated function that automatically choose best lambda that result in better behave q-q plot
  if (is.null(ncol(input_pvalues)))
    stop("input_pvalues should be a matrix or data frame")
  if (ncol(input_pvalues) != 2)
    stop("inpute_pvalues should have 2 column")
  input_pvalues <- matrix(as.numeric(input_pvalues), nrow = nrow(input_pvalues))
  if (sum(stats::complete.cases(input_pvalues)) < nrow(input_pvalues))
    warning("input_pvalues contains NAs to be removed from analysis")
  input_pvalues <- input_pvalues[stats::complete.cases(input_pvalues),
  ]
  if (!is.null(nrow(input_pvalues)) & nrow(input_pvalues) <
      1)
    stop("input_pvalues doesn't have valid p-values")

  ### first identify features that may come from H11 (alternatives for both hypotheses) ###
  #library(qvalue)
  #ish11 <- qvalue(input_pvalues[,1])$qvalue<0.25 & qvalue(input_pvalues[,2])$qvalue<0.25

  pcut <- seq(0.1, 0.8, 0.1)
  frac1 <- rep(0, 8)
  frac2 <- rep(0, 8)
  frac12 <- rep(0, 8)
  for (i in 1:8) {
    frac1[i] <- mean(input_pvalues[, 1] >= pcut[i])/(1 - pcut[i])
    frac2[i] <- mean(input_pvalues[, 2] >= pcut[i])/(1 - pcut[i])
    frac12[i] <- mean(input_pvalues[, 2] >= pcut[i] & input_pvalues[,1] >= pcut[i])/(1 - pcut[i])^2
  }

  alphaout <- matrix(0,4,5)
  ll <- 1
  qqslope <- rep(0,4)
  for (lambda in c(0.5,0.6,0.7,0.8)) {
    alpha00 <- min(frac12[pcut >= lambda][1], 1)
    if (stats::ks.test(input_pvalues[, 1], "punif", 0, 1, alternative = "greater")$p > 0.05)
      alpha1 <- 1 else alpha1 <- min(frac1[pcut >= lambda][1], 1)
      if (stats::ks.test(input_pvalues[, 2], "punif", 0, 1, alternative = "greater")$p > 0.05)
        alpha2 <- 1 else alpha2 <- min(frac2[pcut >= lambda][1], 1)
        if (alpha00 == 1) {
          alpha01 <- 0
          alpha10 <- 0
          alpha11 <- 0
        } else {
          if (alpha1 == 1 & alpha2 == 1) {
            alpha01 <- 0
            alpha10 <- 0
            alpha11 <- 0
            alpha00 <- 1
          }
          if (alpha1 == 1 & alpha2 != 1) {
            alpha10 <- 0
            alpha11 <- 0
            alpha01 <- alpha1 - alpha00
            alpha01 <- max(0, alpha01)
            alpha00 <- 1 - alpha01
          }
          if (alpha1 != 1 & alpha2 == 1) {
            alpha01 <- 0
            alpha11 <- 0
            alpha10 <- alpha2 - alpha00
            alpha10 <- max(0, alpha10)
            alpha00 <- 1 - alpha10
          }
          if (alpha1 != 1 & alpha2 != 1) {
            alpha10 <- alpha2 - alpha00
            alpha10 <- max(0, alpha10)
            alpha01 <- alpha1 - alpha00
            alpha01 <- max(0, alpha01)
            if ((1 - alpha00 - alpha01 - alpha10) < 0) {
              alpha11 <- 0
              alpha10 <- 1 - alpha1
              alpha01 <- 1 - alpha2
              alpha00 <- 1 - alpha10 - alpha01
            }
            else {
              alpha11 <- 1 - alpha00 - alpha01 - alpha10
            }
          }
        }

        pmax <- apply(input_pvalues,1,max)
        pmax <- pmax[order(pmax)]
        nnulls <- sum(pmax>0.8)
        nmed <- nrow(input_pvalues)
        pexp <- rep(0,nnulls)
        for (i in 1:nmed) {
          c <- (-i/nmed)
          b <- alpha01+alpha10-eps
          a <- 1-b
          pexp[i] <- (-b+sqrt(b^2-4*a*c))/(2*a)
        }
        xx <- -log(pexp[(nmed-nnulls+1):nmed],base=10)
        yy <- -log(pmax[(nmed-nnulls+1):nmed],base=10)
        fit1 <- stats::lm(yy~xx-1)

        qqslope[ll]<- fit1$coef[1]
        alphaout[ll,1] <- alpha10
        alphaout[ll,2] <- alpha01
        alphaout[ll,3] <- alpha00
        alphaout[ll,4] <- alpha1
        alphaout[ll,5] <- alpha2

        ll <- ll+1

  }

  bestslope <- which.min(qqslope)
  alpha.null <- list(alpha10 = alphaout[bestslope,1], alpha01 = alphaout[bestslope,2],alpha00 = alphaout[bestslope,3], alpha1 = alphaout[bestslope,4], alpha2 = alphaout[bestslope,5])
  return(alpha.null)
}




###########################################################################################
#### 17. function to conduct the joint significance test as the existing approaches do ####
###########################################################################################

# Inputs
# chosen_ind: a vector of indices for chosen mediators
# pval_alp: a vector of RAW p-values to test H_{0}: \alpha_{j} = 0 where \alpha_{j} is the effect of the exposure on the j-th mediator
# pval_beta: a vector of RAW p-values to test H_{0}: \beta_{j} = 0 where \beta_{j} is the effect of the j-th mediator on the response
# method: a method used to adjust each p-value that will be plugged as an argument of the p.adjust function, see ?p.adjust for more details. HDMT is the default
# alpha: a significance level
# output
# which_sig: an index set for active mediators that turn out to be significant

js_test <- function(chosen_ind, pval_alp, pval_beta, method=NULL, alpha = 0.05){


      PA <- cbind(pval_alp, pval_beta)
      P_value <- apply(PA, 1, max) # The joint p-values for SIS variable

      N0 <- dim(PA)[1] * dim(PA)[2]
      input_pvalues <- PA + matrix(stats::runif(N0, 0, 10^{-10}), dim(PA)[1], 2)

      # Estimate the proportions of the three component nulls
      nullprop <- null_estimation(input_pvalues)

      fdrcut <- HDMT::fdr_est(nullprop$alpha00,
                              nullprop$alpha01,
                              nullprop$alpha10,
                              nullprop$alpha1,
                              nullprop$alpha2,
                              input_pvalues,
                              exact = 0
      )

      which_sig <- which(fdrcut <= alpha)


  if(method == 'bonferroni'){
       m <- length(pval_beta)
       adjp_alp <- stats::p.adjust(pval_alp, method = method)
       adjp_beta <- stats::p.adjust(pval_beta, method = method)
       max_pval <- sapply(c(1:m), FUN = function(j){max(adjp_alp[j], adjp_beta[j])})
       which_sig <- chosen_ind[max_pval <= alpha]
  }



  return(which_sig)
}

