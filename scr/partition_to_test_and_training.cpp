#include <armadillo>
#include <containers>
#include <algorithm> //for std::sample
#include <math.h>       /* floor */
#include <RcppArmadilloExtensions/sample.h>


//' (Pseudo-)Randomly sample indices to determine a test set
//' 
//' @param n_obs total number of subjects (test plus training)
//' @param test_proportion proportion of subjects to place in test set
//' @return  integers to indicate test set membership

arma::Col<uword> get_test_indices(int n_obs, double test_proportion){
  // calculate number of subjects to put into test set
  int n_test = floor(test_proportion * n_obs);
  arma::Col<uword> result = arma::randperm(n_obs, n_test); //randomly sample without replacement from the integers 0,1,...,n_obs - 1 and return n_test of them.
  return(result);
}

//' Subset a matrix's rows by indices 
//' 
//' @param mat a matrix, eg., of genotypes, for the entire cohort, with one subject per row
//' @param test_indices vector with subject indices to go into test set
//' @return matrix of genotypes for the subsetted collection of subjects

arma::mat subset(arma::mat matrix, arma::Col<uword> indices){
  arma::mat result = matrix.rows(indices);
  return(result);
}

//' Subset a vector by indices
//' 
//' @param vector a vector
//' @param indices arma::vec of indices to indicate which entries to extract 
//' @return vector of values for the subsetted collection of subjects

arma::vec subset(arma::vec vector, arma::Col<uword> indices){
  arma::vec result = vector.elem(indices);
  return(result);
}
