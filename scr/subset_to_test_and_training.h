#include <armadillo>
#include <math.h>       /* floor */

std::vector<int> fisher_yates_shuffle(std::size_t size, 
                                      std::size_t max_size, 
                                      std::mt19937& gen);

arma::Col<arma::uword> get_test_indices(int n_obs, double test_proportion);

arma::mat subset(arma::mat matrix, arma::Col<arma::uword> indices);

arma::vec subset(arma::vec vector, arma::Col<arma::uword> indices);

arma::Col<arma::uword> get_training_indices(arma::Col<arma::uword> test_indices, int sample_size);  




