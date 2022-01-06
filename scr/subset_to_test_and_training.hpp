#include <armadillo>
#include <math.h>       /* floor */
#include <algorithm> /* std::sort, std::set_difference */
#include <fstream> //std::ifstream
#include <string>




arma::mat subset(arma::mat matrix, arma::Col<uint> indices);

arma::vec subset(arma::vec vector, arma::Col<uint> indices);

arma::Col<uint> get_complementary_indices(arma::Col<uint> indices, int sample_size);  

arma::Col<uint> convert_string_to_Col(std::vector<std::string> string);

arma::Col <uint> read_indices_file(const std::string filepath);


