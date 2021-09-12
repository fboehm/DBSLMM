#define ARMA_DONT_USE_WRAPPER

#include <iostream>
#include <vector>
#include <armadillo>
#include <string>
#include <boost/algorithm/string.hpp> // split
#include <tuple> // std::tuple, std::get, std::tie, std::ignore


arma::mat calc_asymptotic_variance(arma::mat Xl_training, 
                                   arma::mat Xs_training, 
                                   arma::mat Xl_test, 
                                   arma::mat Xs_test,
                                   double sigma2_s,
                                   arma::vec y_training);
arma::mat calc_Hinv(arma::mat Xs_training, 
                    double sigma2_s);
arma::mat calc_var_betal(arma::mat Xl, 
                      arma::mat Hinv, 
                      arma::vec y);
  

arma::mat calc_var_betas(arma::mat Xl, 
                         arma::mat Xs,
                         arma::mat Hinv,
                         double sigma2_s,
                         arma::vec y,
                         arma::mat var_bl);


std::tuple<vector<string>, vector<string> > read_pheno(std::string file_path, 
                                                       int column_number = 6);

std::vector<double> convert_string_vector_to_double_vector(vector<string> string_vector);

