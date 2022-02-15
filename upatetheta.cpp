#include <RcppArmadillo.h>

// [[Rcpp::depends(RcppArmadillo)]]

// [[Rcpp::export]]

arma::mat updateTheta (arma::mat X,
               arma::mat Y,
               arma::mat P, 
               double rho,
               arma::vec Q, 
               arma::vec U
) {
  arma::mat betahat ;
  betahat = (X.t() * X + rho*P).i() * (X.t() * Y+rho*(Q-U)) ;
  return(betahat) ;
}

// You can include R code blocks in C++ files processed with sourceCpp
// (useful for testing and development). The R code will be automatically 
// run after the compilation.
//

/*** R
set.seed(100)
X <- matrix(rnorm(4*4), 4, 4)
Y <- matrix(1,4,1)
P = diag(4)
rho = 1
Q <-  seq(1,4)
U<- rep(1,4)
updateTheta(X= X, Y = Y,P=P,rho = rho,Q = Q,U = U)
  */
