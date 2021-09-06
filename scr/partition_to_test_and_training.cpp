#include <armadillo>
#include <containers>
#include <algorithm> //for std::sample
#include <math.h>       /* floor */
#include <RcppArmadilloExtensions/sample.h>


//' (Pseudo-)Randomly sample indices to determine a test set
//' 
//' @param n_obs total number of subjects (test plus training)
//' @param test_proportion proportion of subjects to place in test set
//' @return arma::vec containing integers to indicate test set membership

arma::vec get_test_indices(int n_obs, double test_proportion){
  // calculate number of subjects to put into test set
  int n_test = floor(test_proportion * n_obs);
  arma::vec n_vec = arma::regspace(0, n_obs - 1); //create vector 0,1,2,..., n_obs - 1
  //sample the n_test subjects
  arma::vec test_indices = Rcpp::RcppArmadillo::sample(n_vec, n_test, false); // vector, size, replace
  return(test_indices);
}

//' Subset entire set of subjects' genotype data 
//' 
//' @param mat a matrix, eg., of genotypes, for the entire cohort, with one subject per row
//' @param test_indices vector with subject indices to go into test set
//' @return 

arma::mat subset_matrix(arma::mat matrix, arma::vec indices){
  arma::mat result = matrix.rows(indices);
  return(result);
}

//' Subset a vector by indices
//' 
//' @param vector a vector
//' @param indices arma::vec of indices to indicate which entries to extract 
//' @return an arma::vec with a subset of the entries in vector

arma::vec subset_vector(arma::vec vector, arma::vec indices){
  arma::vec result = vector.elem(indices);
  return(result);
}
