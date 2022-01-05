
#include <armadillo>
#include "calc_asymptotic_variance.hpp"

using namespace std;
using namespace arma;

//' Calculate for one LD block the n tilde by n tilde matrix that contributes to the asymptotic variance for the predicted y values
//' 
//' @param Sigma_ll Sigma_ll matrix for a single LD block
//' @param Sigma_ls Sigma_ls matrix for a single LD block
//' @param Sigma_ss Sigma_ss matrix for a single LD block
//' @param sigma2_s estimated value of sigma^2_s
//' @param n_training sample size for training data
//' @param Xl_test genotypes matrix - for the LD block - for large effect SNPs for test subjects
//' @param Xs_test genotypes matrix - for the LD block - for small effect SNPs for test subjects
//' @return n tilde by n tilde matrix that contributes to the asymptotic variance

arma::mat calc_nt_by_nt_matrix(arma::mat Sigma_ll, 
                                   arma::mat Sigma_ls, 
                                   arma::mat Sigma_ss,
                                   double sigma2_s, 
                                   unsigned int n_training,
                                   arma::mat Xl_test, 
                                   arma::mat Xs_test){
  arma::mat Ainv = calc_A_inverse(Sigma_ss, sigma2_s, n_training);
  arma::mat var_bl = calc_var_betal(Sigma_ll, 
                                    Sigma_ls, 
                                    Ainv, 
                                    n_training);
  arma::mat var_bs = calc_var_betas(Sigma_ss, 
                                    Sigma_ls,
                                    Ainv,
                                    sigma2_s,
                                    n_training,
                                    var_bl);
  arma::mat result = Xl_test * var_bl * arma::trans(Xl_test) + Xs_test * var_bs * arma::trans(Xs_test);
  return(result);
}


// for blocks without large effects, via C++ "overloading" 
arma::mat calc_nt_by_nt_matrix(arma::mat Sigma_ss,
                               double sigma2_s, 
                               unsigned int n_training,
                               arma::mat Xs_test){
  arma::mat Ainv = calc_A_inverse(Sigma_ss, sigma2_s, n_training);
  arma::mat var_bs = calc_var_betas(Sigma_ss, 
                                    Ainv,
                                    sigma2_s,
                                    n_training);
  arma::mat result = Xs_test * var_bs * arma::trans(Xs_test);
  return(result);
}


//' Calculate A inverse matrix
//' 
//' @details (sigma^{-2}n^{-1} I_ms + Sigma_ss) = A. 
//' @param Sigma_ss
//' @param sigma2_s estimate of sigma^2_s
//' @param n sample size
//' @return an Armadillo matrix, the inverse of A

arma::mat calc_A_inverse(arma::mat Sigma_ss, 
                         double sigma2_s, 
                         unsigned int n)
{
  unsigned int m_s = Sigma_ss.n_rows;
  arma::mat result = arma::inv_sympd(arma::eye(m_s, m_s) / (n * sigma2_s) + Sigma_ss);
  return(result);
}

//' Calculate variance of coefficient estimator for large effects
//' 
//' @param Sigma_ll Sigma_ll constructed for one LD block
//' @param Sigma_ls Sigma_ls constructed for one LD block
//' @param Sigma_ss Sigma_ss constructed for one LD block
//' @param A_inverse inverse of (sigma^{-2}n^{-1} I_ms + Sigma_ss)
//' @param n sample size (of training set)
//' @return covariance matrix

arma::mat calc_var_betal(arma::mat Sigma_ll, 
                         arma::mat Sigma_ls, 
                         arma::mat A_inverse,
                         unsigned int n){

  //calculate big matrix
  arma::mat big = Sigma_ll - Sigma_ls * A_inverse * arma::trans(Sigma_ls);
  //invert and divide by n
  arma::mat result = arma::inv_sympd(big) / n; // we invert a ml by ml matrix - no problem!
  return (result);
}

//' Calculate variance of coefficient estimator for small effects
//' 
//' @param Sigma_ss Sigma_ss matrix 
//' @param Sigma_sl Sigma_sl matrix 
//' @param A_inverse A inverse matrix 
//' @param sigma2_s estimated value of sigma^2_s
//' @param n sample size
//' @param var_bl variance of beta hat l
//' @return covariance matrix

arma::mat calc_var_betas(arma::mat Sigma_ss, 
                         arma::mat Sigma_sl,
                         arma::mat A_inverse,
                         double sigma2_s,
                         unsigned int n,
                         arma::mat var_bl){
  arma::mat mat1 = Sigma_ss - Sigma_ss * A_inverse * Sigma_ss;
  arma::mat mat2 = Sigma_sl - Sigma_ss * A_inverse * Sigma_sl;
  arma::mat result = n * sigma2_s * sigma2_s * (mat1 + mat2 * n * var_bl * arma::trans(mat2));
  return (result);
}

// for LD blocks without large effects:
arma::mat calc_var_betas(arma::mat Sigma_ss, 
                         
                         arma::mat A_inverse,
                         double sigma2_s,
                         unsigned int n
                         ){
  arma::mat mat1 = Sigma_ss - Sigma_ss * A_inverse * Sigma_ss;
  arma::mat result = n * sigma2_s * sigma2_s * mat1;
  return (result);
}

