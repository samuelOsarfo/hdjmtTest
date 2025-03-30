#include <RcppArmadillo.h>
// [[Rcpp::depends(RcppArmadillo)]]

//' @title Fast Approximate Orthogonalization :: Using C++
 //' @name app_orth_rcpp
 //' @description Rcpp implementation of approximate orthogonalization
 //' @param y Response vector
 //' @param MCX Combined matrix of predictors
 //' @param sig_hat Estimated sigma
 //' @param k Regularization parameter
 //' @param q2 Number of covariates
 //' @return List with test statistics and p-values
 //' @keywords internal
 //' @export
 // [[Rcpp::export]]
 Rcpp::List app_orth_rcpp(Rcpp::NumericVector y,
                          Rcpp::NumericMatrix MCX,
                          double sig_hat,
                          double k,
                          int q2) {
   arma::vec y_vec(y.begin(), y.size(), false);
   arma::mat MCX_mat(MCX.begin(), MCX.nrow(), MCX.ncol(), false);

   int n = y_vec.n_elem;
   int q1 = MCX_mat.n_cols;
   int n_tests = q1 - (q2 + 1);

   arma::vec ts(n_tests);
   arma::vec pval(n_tests);
   arma::mat kI = k * arma::eye(n, n);

   for(int j = 0; j < n_tests; j++) {
     arma::vec X_j = MCX_mat.col(j);
     arma::mat MCX_j = MCX_mat;
     MCX_j.shed_col(j);

     //arma::mat proj_vec = arma::solve(kI + MCX_j * MCX_j.t(), X_j,
      //                                arma::solve_opts::fast + arma::solve_opts::likely_sympd);
      arma::mat proj_vec = arma::solve(kI + MCX_j * MCX_j.t(), X_j, arma::solve_opts::none);

     double proj_norm = arma::norm(proj_vec, 2);
     ts(j) = arma::dot(proj_vec, y_vec) / (sig_hat * proj_norm);
     pval(j) = 2 * R::pnorm(std::abs(ts(j)), 0.0, 1.0, false, false);
   }

   return Rcpp::List::create(
     Rcpp::Named("ts") = ts,
     Rcpp::Named("pval") = pval
   );
 }
