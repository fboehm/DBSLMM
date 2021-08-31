#define ARMA_DONT_USE_WRAPPER

#include <iostream>
#include <vector>
#include <armadillo>

//' Calculate the asymptotic variance for the beta hat estimator
//' 
//' @param Xl genotypes matrix for large effects markers
//' @param Xs genotypes matrix for small effects markers
//' @param sigma2_s estimated value of sigma^2_s
//' @param y trait values vector
//' @return variance of predicted y values

arma::mat calc_asymptotic_variance(arma::mat Xl_training, 
                                   arma::mat Xs_training, 
                                   arma::mat Xl_test, 
                                   arma::mat Xs_test,
                                   double sigma2_s,
                                   arma::vec y){
  arma::mat Hinv = calc_Hinv(Xs_training, sigma2_s);
  
  return(result);
}

//' Calculate H inverse matrix (with Woodbury identity)
//' 
//' @param Xs_training
//' @param sigma2_s estimated value of sigma^2_s
//' @return H inverse matrix

arma::mat calc_Hinv(arma::mat Xs_training, 
                    double sigma2_s){
  int n = Xs_training.n_rows;
  int ms = Xs_training.n_cols;
  arma::mat result = arma::eye(n, n) - Xs_training * arma::inv(arma::eye(ms, ms) / sigma2_s + Xs_training.t() * Xs_training) * Xs_training.t();
  return(result);
}

//' Calculate variance of coefficient estimator for large effects
//' 
//' @param Xl matrix of genotypes for large effect markers
//' @param Hinv inverse of H matrix
//' @param y trait values vector
//' @return covariance matrix

arma::mat calc_var_betal(arma::mat Xl, 
                      arma::mat Hinv, 
                      arma::vec y){
  //var y
  arma::mat vy = y * y.t(); //check this!! only true if mean(y) = 0 vector.
  // (Xl^T Hinv Xl)^{-1}
  arma::mat arg1 = arma::inv(Xl.t() * Hinv * Xl);
  arma::mat result = arg1 * Xl.t() * Hinv * vy * Hinv * Xl * arg1;
  return result;
}
  
