#include <armadillo>


#include "dtpr.hpp"

using namespace std;
using namespace arma;

arma::vec count_snps_per_block(int num_block, vector <INFO> info);

unsigned int sum_vec(vector<int> vv);


vector <INFO> populate_vector_vector_info(vector <INFO> info, int max_length);

INFO make_pseudo_info();

EFF make_pseudo_eff();

arma::mat populate_geno(string bed_str, 
                        vector<int> idv, 
                        vector <INFO> info_block,
                        double maf = 0.0);

vector<INFO> populate_info_block(vector <INFO> info_block_full, int num_block);


