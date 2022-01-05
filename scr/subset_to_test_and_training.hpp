#include <armadillo>
#include <math.h>       /* floor */

std::vector<int> fisher_yates_shuffle(std::size_t size, 
                                      std::size_t max_size, 
                                      std::mt19937& gen);

arma::Col<arma::uword> get_test_indices(int n_obs, 
                                        double test_proportion, 
                                        unsigned int seed);

arma::mat subset(arma::mat matrix, arma::Col<arma::uword> indices);

arma::vec subset(arma::vec vector, arma::Col<arma::uword> indices);

arma::Col<arma::uword> get_complementary_indices(arma::Col<arma::uword> indices, int sample_size);  

arma::Col<arma::uword> convert_string_to_Col(std::vector<std::string> string);

arma::Col <arma::uword> read_indices_file(const char filepath);


