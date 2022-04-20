#include <armadillo>


#include "dtpr.hpp"

using namespace std;
using namespace arma;

arma::vec count_snps_per_block(int num_block, vector <INFO> info);

unsigned int sum_vec(vector<int> vv);

vector <INFO> subset_info_by_block(vector <INFO> info, int block_num);

vector <INFO> populate_vector_vector_info(vector <INFO> info, int max_length);

INFO make_pseudo_info();

