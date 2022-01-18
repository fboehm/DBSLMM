
#include <armadillo>

arma::mat calc_nt_by_nt_matrix(arma::mat Sigma_ll, 
                               arma::mat Sigma_ls, 
                               arma::mat Sigma_ss,
                               double sigma2_s, 
                               unsigned int n_training,
                               arma::mat Xl_test, 
                               arma::mat Xs_test);
arma::mat calc_nt_by_nt_matrix(arma::mat Sigma_ss,
                               double sigma2_s, 
                               unsigned int n_training,
                               arma::mat Xs_test);
arma::mat calc_A_inverse(arma::mat Sigma_ss, 
                         double sigma2_s, 
                         unsigned int n);
arma::mat calc_var_betal(arma::mat Sigma_ll, 
                         arma::mat Sigma_ls, 
                         arma::mat A_inverse,
                         unsigned int n);
arma::mat calc_var_betas(arma::mat Sigma_ss, 
                         arma::mat Sigma_sl, //check this!
                         arma::mat A_inverse,
                         double sigma2_s,
                         unsigned int n,
                         arma::mat var_bl);
arma::mat calc_var_betas(arma::mat Sigma_ss, 
                         arma::mat A_inverse,
                         double sigma2_s,
                         unsigned int n);
  