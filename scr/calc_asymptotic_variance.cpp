#define ARMA_DONT_USE_WRAPPER

#include <iostream>
#include <vector>
#include <armadillo>
#include "calc_asymptotic_variance.h"

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
                                   arma::vec y_training){
  arma::mat Hinv = calc_Hinv(Xs_training, sigma2_s);
  arma::mat var_bl = calc_var_betal(Xl = Xl_training, 
                                    Hinv = Hinv, 
                                    y = y_training);
  arma::mat var_bs = calc_var_betas(Xl = Xl_training, 
                                    Xs = Xs_training, 
                                    Hinv = Hinv, 
                                    sigma2_s = sigma2_s,
                                    y = y_training, 
                                    var_bl = var_bl);
  arma::mat result = Xl_test.t() * var_bl * Xl_test + Xs_test.t() * var_bs * Xs_test;
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
  
//' Calculate variance of coefficient estimator for small effects
//' 
//' @param Xs matrix of genotypes for small effect markers
//' @param Xl matrix of genotypes for large effect markers
//' @param Hinv inverse of H matrix
//' @param sigma2_s estimated value of sigma^2_s
//' @param y trait values vector
//' @param var_bl variance of beta hat l
//' @return covariance matrix
  
arma::mat calc_var_betas(arma::mat Xl, 
                         arma::mat Xs,
                         arma::mat Hinv,
                         double sigma2_s,
                         arma::vec y,
                         arma::mat var_bl){
  //var y
  arma::mat vy = y * y.t(); 
  // 
  arma::mat arg1 = sigma2_s * sigma2_s * Xs.t() * Hinv * vy * Hinv * Xs;
  arma::mat result = arg1 + sigma2_s * sigma2_s * Xs.t() * Hinv * Xl * var_bl * Xl.t() * Hinv * Xs;
  return result;
}

