#include <armadillo>
#include <math.h>       /* floor */
#include <algorithm> /* std::sort, std::set_difference */
#include <fstream> //std::ifstream
#include <string>




arma::mat subset(arma::mat matrix, arma::uvec indices);

arma::vec subset(arma::vec vector, arma::uvec indices);

arma::uvec get_complementary_indices(arma::uvec indices, int sample_size);  

arma::uvec convert_string_to_indices(std::string in_string);

arma::uvec read_indices_file(const std::string filepath);

template<typename C1, typename C2>
void castContainer(const C1& source, C2& destination);

template<typename T, typename T2>
std::vector<T>& operator<<(std::vector<T>& v, T2 t);
