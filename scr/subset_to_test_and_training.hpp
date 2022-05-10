#include <armadillo>
#include <math.h>       /* floor */
#include <algorithm> /* std::sort, std::set_difference */
#include <fstream> //std::ifstream
#include <string>
#include <boost/lexical_cast.hpp>
#include <iterator>
#include <regex>
#include <vector>

std::vector<int> make_ones_and_zeroes_vec(arma::uvec ones_positions, unsigned int length);

arma::mat subset(arma::mat matrix, arma::uvec indices);

arma::vec subset(arma::vec vector, arma::uvec indices);

arma::uvec get_complementary_indices(arma::uvec indices, int sample_size);  

arma::uvec convert_string_to_indices(std::string in_string);

std::vector<std::string> split(const std::string str, const std::string regex_str);

arma::uvec read_indices_file(const std::string filepath);

template<typename C1, typename C2>
void castContainer(const C1& source, C2& destination);

template<typename T, typename T2>
std::vector<T>& operator<<(std::vector<T>& v, T2 t);



std::vector<int> findItems(std::vector<int> const &v, int target);

