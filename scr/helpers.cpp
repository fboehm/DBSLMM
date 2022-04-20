#include <armadillo>

#include "helpers.hpp"
#include "dtpr.hpp"

using namespace std;
using namespace arma;


//' Count the number of SNPs in each block on one chromosome
//' 
//' @param num_block number of blocks on the chromosome
//' @param info a vector <INFO> for small or large effect SNPs on that chromosome
//' @return an armadillo vector with one entry per block; that entry is the number of SNPs in the block

arma::vec count_snps_per_block(int num_block, vector <INFO> info){
  int count = 0;
  arma::vec num = zeros<vec>(num_block); 
  for (int i = 0; i < num_block; i++) {
    for (size_t j = count; j < info.size(); j++) {
      if(info[j].block == i){ 
        num(i) += 1; 
        count++;
      }else{
        break;
      }
    }
  }
  return num;
}

//' Sum entries in a vector with only nonnegative integer entries 
//' 
//' @param vv vector<int> to be summed
//' @return unsigned int that is the sum of the entries in the vector

unsigned int sum_vec(vector<int> vv){
  unsigned int result = std::accumulate(vv.begin(), vv.end(),
                  decltype(vv)::value_type(0));
  return result;
}




//' Populate vector of vector of INFO
//' 
//' @param info a vector <INFO> object
//' @return 
//' 
vector <INFO> populate_vector_vector_info(vector <INFO> info, int max_length){
  vector <INFO>  result(max_length);
  for (size_t k = 0; k < info.size(); k++)
    // info_s_Block[B][k] = &info_s_block[k]; 
    result[k] = info[k]; 
  return result;
}
  
//' Prepare pseudo INFO object
//' 
//' @return a pseudo INFO object

INFO make_pseudo_info(){
  INFO info_pseudo; 
  info_pseudo.snp = "rs"; 
  info_pseudo.ps = 0; 
  info_pseudo.pos = 0; 
  info_pseudo.block = 0; 
  info_pseudo.a1 = "Y"; 
  info_pseudo.maf = 0; 
  info_pseudo.z = 0; 
  info_pseudo.P = 0; 
  return info_pseudo;
}
