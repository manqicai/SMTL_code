#include <RcppArmadillo.h>

// [[Rcpp::depends(RcppArmadillo)]]

// [[Rcpp::export]]

Rcpp::List admm(Rcpp::List X, Rcpp::List Y,
          arma::mat P, 
          double rho,
          arma::mat Q, 
          arma::mat U, int n
            
) {
  Rcpp::List betahat(n) ;
  for (int k=0; k < n; ++k) {
    
    //arma::mat betahat[k] ;
    //arma::mat Xk = Rcpp::as<arma::mat>(X[k]);
    //arma::mat Yk = Rcpp::as<arma::mat>(Y[k]);
    //arma::vec Qk = Rcpp::as<arma::vec>(Q[k]);
   // arma::vec Uk = Rcpp::as<arma::vec>(U[k]);
   Rcpp::NumericMatrix x = X[k];
   arma::mat Xk = arma::mat(x.begin(), x.nrow(), x.ncol(), false);
   Rcpp::NumericMatrix y = Y[k];
   arma::mat Yk = arma::mat(y.begin(), y.nrow(), y.ncol(), false);
   
   betahat[k] = (Xk.t() * Xk + rho*P).i() * (Xk.t() * Yk+rho*(Q.col(k)-U.col(k))) ;
    //betahat[k] = (Xk.t() * Xk + rho*P).i() * (Xk.t() * Yk+rho*(Qk-Uk)) ;
  }
  
  return(betahat) ;
}



/*** R
set.seed(100)
X <- list(a = matrix(rnorm(4*4), 4, 4), b =matrix(rnorm(4*4), 4, 4))
Y <- list(matrix(1,4,1),matrix(1,4,1))
P = diag(4)
rho = 1
Q <-  cbind(seq(1,4),seq(1,4))
U<- cbind(rep(1,4),seq(1,4))
K = 2
admm(X= X, Y = Y,P=P,rho = rho,Q = Q,U = U,n = K)
*/