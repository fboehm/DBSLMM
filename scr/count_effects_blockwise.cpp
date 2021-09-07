#define ARMA_DONT_USE_WRAPPER

#include <armadillo>
#include <vector> // for counting block memberships
#include <unordered_map> // for counting block memberships
#include "dbslmmfit.hpp"
#include "dtpr.hpp"



//' Count number of SNPs with effects - large or small (choose one per call) - across the genome. 
//' 
//' @param info an info_s or info_l object for the small or large effect snps
//' @param num_block a nonnegative integer specifying the number of genomewide LD blocks
//' @return a vector with one entry per block. Each entry is the nonnegative integer number of effects (small or large) for each block

arma::Col<arma::uword> count_snps_by_effect_size(vector <INFO> info, int num_block){
  arma::Col<arma::uword> block_membership = extract_block_membership(info);
  // convert block_membership to std::vector
  std::vector bm_vec ;
  // count element occurrences in bm_vec
  std::unordered_map<int, int> freq;
  for (int const &i: bm_vec) {
    freq[i]++;
  }
  //get a vector of counts from block_membership
  return(result);
} 

//' Extract block memberships from info
//' 
//' @param info an info object for small or large effects (choose one per call)
//' @return a vector, with one integer entry per SNP. It specifies the LD block membership of each SNP 

arma::Col<arma::uword> extract_block_membership(vector <INFO> info){
  arma::Col<arma::uword> result;
  for (size_t j = 0; j < info.size(); j++) {
    result[j] = info[j].block;
  }
  return(result);
} 


