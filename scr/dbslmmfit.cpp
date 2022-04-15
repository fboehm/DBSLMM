/*
Deterministic Bayesian Sparse Linear Mixed Model (DBSLMM)
Copyright (C) 2019  Sheng Yang and Xiang Zhou

This program is free software: you can redistribute it and/or modify
it under the terms of the GNU General Public License as published by
the Free Software Foundation, either version 3 of the License, or
(at your option) any later version.

This program is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
GNU General Public License for more details.

You should have received a copy of the GNU General Public License
along with this program.  If not, see <http://www.gnu.org/licenses/>.
*/

#define ARMA_DONT_USE_WRAPPER

#include <iostream>
#include <vector>
#include <armadillo>
#include <string>
#include <numeric>
#include "omp.h"

#include "dtpr.hpp"
#include "dbslmmfit.hpp"
#include "calc_asymptotic_variance.hpp"
#include "subset_to_test_and_training.hpp"

using namespace std;
using namespace arma;

//' Estimate large and small effects
//' 
//' @param n_ref sample size of the reference panel
//' @param n_obs sample size of the nonreference data
//' @param sigma_s the estimate for sigma_s^2
//' @param num_block number of blocks in the genome
//' @param idv indicator of trait missingness for the reference panel data. Since we don't care about trait values for the reference data, this is typically a vector of all 1's to indicate no missingness 
//' @param bed_str file path for the plink bed file
//' @param info_s small effect SNP info object
//' @param info_l large effect SNP info object
//' @param thread number of threads to use
//' @param eff_s small effects SNP effects object
//' @param eff_l large effects SNP effects object
//' @param test_indicator vector containing ones and zeroes for indicating which subjects are in the test set.
//' @param genotypes_str a string containing the file path to the file containing genotypes for test and training subjects
//' @return zero is returned
// estimate large and small effect
int  DBSLMMFIT::est(int n_ref, 
                    int n_obs, 
                    double sigma_s, 
                    int num_block, 
                    vector<int> idv, 
                    string bed_str,
				            vector <INFO> info_s, 
				            vector <INFO> info_l, 
				            int thread, 
                    vector <EFF> &eff_s, 
                    vector <EFF> &eff_l,
                    vector<int> test_indicator,
                    string genotypes_str, 
                    vector <INFO> test_info_s, 
                    vector <INFO> test_info_l, 
                    bool badsnp_s[], 
                    bool badsnp_l[]){
	
	// get the maximum number of each block
	int count_s = 0, count_l = 0;
	vec num_s = zeros<vec>(num_block), num_l = zeros<vec>(num_block); 
	for (int i = 0; i < num_block; i++) {
		for (size_t j = count_s; j < info_s.size(); j++) {
			if(info_s[j].block == i){ 
				num_s(i) += 1; 
				count_s++;
			}else{
				break;
			}
		}
		for (size_t j = count_l; j < info_l.size(); j++) {
			if(info_l[j].block == i){ 
				num_l(i) += 1; 
				count_l++;
			}else{
				break;
			}
		}
	}// end of counting and populating num_l and num_s, ie, iterates over i
	count_l = count_s = 0; // reset
	// get the maximum number of each block
	vec test_num_s = zeros<vec>(num_block), test_num_l = zeros<vec>(num_block); 
	for (int i = 0; i < num_block; i++) {
	  for (size_t j = count_s; j < test_info_s.size(); j++) {
	    if(test_info_s[j].block == i){ 
	      test_num_s(i) += 1; 
	      count_s++;
	    }else{
	      break;
	    }
	  }
	  for (size_t j = count_l; j < test_info_l.size(); j++) {
	    if(test_info_l[j].block == i){ 
	      test_num_l(i) += 1; 
	      count_l++;
	    }else{
	      break;
	    }
	  }
	}

	count_l = count_s = 0; // reset
	
	double len_l = num_l.max(); //len_l is the maximum, across blocks, of the per-block number of large effect SNPs
	double len_s = num_s.max(); //len_s is the max of the per-block number of small effect SNPs
	double test_len_l = test_num_l.max();
	double test_len_s = test_num_s.max();

	
	int B = 0;
	int B_MAX = 60; // number of blocks to batch at one time
	if (num_block < 60){ //num_block is the number of blocks across the genome.
		B_MAX = num_block; 
	}

	// pseudo INFO
	INFO info_pseudo; 
	info_pseudo.snp = "rs"; 
	info_pseudo.ps = 0; 
	info_pseudo.pos = 0; 
	info_pseudo.block = 0; 
	info_pseudo.a1 = "Y"; 
	info_pseudo.maf = 0; 
	info_pseudo.z = 0; 
	info_pseudo.P = 0; 
	
	// loop 
	// vector < vector <INFO*> > info_s_Block(B_MAX, vector <INFO*> ((int)len_s)), info_l_Block(B_MAX, vector <INFO*> ((int)len_l));

	vector < vector <INFO> > info_s_Block(B_MAX, vector <INFO> ((int)len_s)), info_l_Block(B_MAX, vector <INFO> ((int)len_l));
	vector < vector <INFO> > test_info_s_Block(B_MAX, vector <INFO> ((int)test_len_s)), test_info_l_Block(B_MAX, vector <INFO> ((int)test_len_l));
	vector < vector <EFF> > eff_s_Block(B_MAX, vector <EFF> ((int)len_s)), eff_l_Block(B_MAX, vector <EFF> ((int)len_l));
	vector <int> num_s_vec, num_l_vec;
	vector <int> test_num_s_vec, test_num_l_vec;
	
  unsigned int n_test = std::accumulate(test_indicator.begin(), test_indicator.end(),
                                        decltype(test_indicator)::value_type(0));
	arma::vec diags = zeros(n_test); //specify length of diags
	

	for (int i = 0; i < num_block; ++i) {
		// small effect SNP information
		vector <INFO> info_s_block; 
	  vector <INFO> test_info_s_block; 
		for (size_t j = count_s; j < info_s.size(); j++) {
			if(info_s[j].block == i){ 
				info_s_block.push_back(info_s[j]);
				count_s++;
			}else{
				break;
			}
		} //end loop for populating info_s_block
		count_s = 0; //reset
		for (size_t j = count_s; j < test_info_s.size(); j++) {
		  if(test_info_s[j].block == i){ 
		    test_info_s_block.push_back(test_info_s[j]);
		    count_s++;
		  }else{
		    break;
		  }
  } //end loop for populating test_info_s_block
		for (size_t k = 0; k < info_s_block.size(); k++)
			// info_s_Block[B][k] = &info_s_block[k]; 
			info_s_Block[B][k] = info_s_block[k]; 
		num_s_vec.push_back((int)num_s(i));
		for (size_t k = 0; k < test_info_s_block.size(); k++)
		  // info_s_Block[B][k] = &info_s_block[k]; 
		  test_info_s_Block[B][k] = test_info_s_block[k]; 
		test_num_s_vec.push_back((int)test_num_s(i)); 
		// large effect SNP information
		if (num_l(i) == 0){
			// info_l_Block[B][0] = &info_pseudo;
			info_l_Block[B][0] = info_pseudo;
		}else{
			vector <INFO> info_l_block;
			for (size_t j = count_l; j < info_l.size(); j++) {
				if(info_l[j].block == i){ 
					info_l_block.push_back(info_l[j]);
					count_l++;
				}else{
					break;
				}
			}
			for (size_t k = 0; k < info_l_block.size(); k++)
				// info_l_Block[B][k] = &info_l_block[k]; 
				info_l_Block[B][k] = info_l_block[k]; 
		}
		num_l_vec.push_back((int)num_l(i));
		count_l = 0; //reset

    if (test_num_l(i) == 0){
      // info_l_Block[B][0] = &info_pseudo;
      test_info_l_Block[B][0] = info_pseudo;
    }else{
      vector <INFO> test_info_l_block;
      for (size_t j = count_l; j < test_info_l.size(); j++) {
        if(test_info_l[j].block == i){ 
          test_info_l_block.push_back(test_info_l[j]);
          count_l++;
        }else{
          break;
        }
      }
      for (size_t k = 0; k < test_info_l_block.size(); k++)
        // info_l_Block[B][k] = &info_l_block[k]; 
        test_info_l_Block[B][k] = test_info_l_block[k]; 
    }
    test_num_l_vec.push_back((int)test_num_l(i));

//		
		B++;
		if (B == B_MAX || i + 1 == num_block) { // process the block of SNPs using multi-threading
			omp_set_num_threads(thread);
#pragma omp parallel for schedule(dynamic)
			for (int b = 0; b < B; b++){
			  cout << "call calcBlock... b has value:"<< b << endl;
			  diags += calcBlock(n_ref, 
                                n_obs, 
                                sigma_s, 
                                idv, 
                                bed_str, 
                                info_s_Block[b], 
                                info_l_Block[b],
                                num_s_vec[b], 
                                num_l_vec[b], 
                                eff_s_Block[b], 
                                eff_l_Block[b], 
                                 test_indicator,
                                 genotypes_str, 
                                 test_info_s_Block[b], 
                                test_info_l_Block[b],
                                test_num_s_vec[b], 
                                test_num_l_vec[b], 
                                 badsnp_s, 
                                 badsnp_l);

			}
			// eff of small effect SNPs
			for (int r = 0; r < B; r++) {
				for (int l = 0; l < num_s_vec[r]; l++){
					eff_s.push_back(eff_s_Block[r][l]);
				}
			}
			// eff of large effect SNPs
			for (int r = 0; r < B; r++){
				for (int l = 0; l < num_l_vec[r]; l++){
					eff_l.push_back(eff_l_Block[r][l]);
				}
			}
			B = 0;
			num_l_vec.clear(); 
			num_s_vec.clear();
		}
	}
	//write the diags object to a file
	diags.save("variance.txt", arma_ascii);
	return 0;
}

// estimate only small effect
int DBSLMMFIT::est(int n_ref, 
                   int n_obs, 
                   double sigma_s, 
                   int num_block, 
                   vector<int> idv, 
                   string bed_str,
				            vector <INFO> info_s, 
				            int thread, 
				            vector <EFF> &eff_s,
				            vector<int> test_indicator,
				            string genotypes_str,  
				            vector <INFO> test_info_s, 
				            bool badsnp_s[]){
	
	// get the maximum number of each block
	int count_s = 0;
	vec num_s = zeros<vec>(num_block); 
	for (int i = 0; i < num_block; i++) {
		for (size_t j = count_s; j < info_s.size(); j++) {
			if(info_s[j].block == i){ 
				num_s(i) += 1; 
				count_s++;
			}else{
				break;
			}
		}
	}
	count_s = 0; // reset
	
	double len_s = num_s.max(); 
	//get counts for test_info_s
	vec test_num_s = zeros<vec>(num_block); 
	for (int i = 0; i < num_block; i++) {
	  for (size_t j = count_s; j < test_info_s.size(); j++) {
	    if(test_info_s[j].block == i){ 
	      test_num_s(i) += 1; 
	      count_s++;
	    }else{
	      break;
	    }
	  }
	} 
  double test_len_s = test_num_s.max(); 
	
	count_s = 0; // reset
	int B = 0;
	int B_MAX = 60;
	if (num_block < 60){
		B_MAX = num_block; 
	}

	// pseudo INFO
	INFO info_pseudo; 
	info_pseudo.snp = "rs"; 
	info_pseudo.ps = 0; 
	info_pseudo.pos = 0; 
	info_pseudo.block = 0; 
	info_pseudo.a1 = "Y"; 
	info_pseudo.maf = 0; 
	info_pseudo.z = 0; 
	info_pseudo.P = 0; 
	
	// loop 
	// vector < vector <INFO*> > info_s_Block(B_MAX, vector <INFO*> ((int)len_s));
	vector < vector <INFO> > info_s_Block(B_MAX, vector <INFO> ((int)len_s));
	vector < vector <INFO> > test_info_s_Block(B_MAX, vector <INFO> ((int)test_len_s));
	vector < vector <EFF> > eff_s_Block(B_MAX, vector <EFF> ((int)len_s));
	vector <int> num_s_vec;
	vector <int> test_num_s_vec;
	
	unsigned int n_test = std::accumulate(test_indicator.begin(), test_indicator.end(),
                                       decltype(test_indicator)::value_type(0));
	cout << "n_test: " << n_test << endl;
	arma::vec diags = zeros(n_test); //specify length of diags & set all to zeros

	for (int i = 0; i < num_block; ++i) {
		// small effect SNP information
		vector <INFO> info_s_block; 
	  vector <INFO> test_info_s_block; 
		for (size_t j = count_s; j < info_s.size(); j++) {
			if(info_s[j].block == i){ 
				info_s_block.push_back(info_s[j]);
				count_s++;
			}else{
				break;
			}
		}
		count_s = 0; //reset
		for (size_t j = count_s; j < test_info_s.size(); j++) {
		  if(test_info_s[j].block == i){ 
		    test_info_s_block.push_back(test_info_s[j]);
		    count_s++;
		  }else{
		    break;
		  }
} 
		// loop to populate info_s_Block
		for (size_t k = 0; k < info_s_block.size(); k++)
			// info_s_Block[B][k] = &info_s_block[k]; 
			info_s_Block[B][k] = info_s_block[k]; 
		for (size_t k = 0; k < test_info_s_block.size(); k++)
		  // info_s_Block[B][k] = &info_s_block[k]; 
		  test_info_s_Block[B][k] = test_info_s_block[k]; 
		num_s_vec.push_back((int)num_s(i));
		test_num_s_vec.push_back((int)test_num_s(i));
		
		B++;
		if (B == B_MAX || i + 1 == num_block) { // process the block of SNPs using multi-threading
		  omp_set_num_threads(thread);
#pragma omp parallel for schedule(dynamic)
			for (int b = 0; b < B; b++){
			  diags += calcBlock(n_ref, 
                            n_obs, 
                            sigma_s, 
                            idv, 
                            bed_str, 
                            info_s_Block[b],
				                    num_s_vec[b], 
                            eff_s_Block[b], 
                           test_indicator, 
                           genotypes_str, 
                           test_info_s_Block[b],
                            test_num_s_vec[b], 
                            badsnp_s);
			  //cout <<"index: " << index << endl; 
			}
			// eff of small effect SNPs
			for (int r = 0; r < B; r++) {
				for (int l = 0; l < num_s_vec[r]; l++){
					eff_s.push_back(eff_s_Block[r][l]);
				}
			}
			B = 0;
			num_s_vec.clear();
		}
	}
	//write the diags object to a file
	diags.save("variance.txt", arma_ascii);
	return 0;
}

// estimate large and small effect for each block
arma::vec DBSLMMFIT::calcBlock(int n_ref, 
                               int n_obs, 
                               double sigma_s, 
                               vector<int> idv, 
                               string bed_str, 
						                    vector <INFO> info_s_block_full, 
						                    vector <INFO> info_l_block_full, 
						                    int num_s_block, 
						                    int num_l_block, 
						                    vector <EFF> &eff_s_block, 
						                    vector <EFF> &eff_l_block,
						                    vector <int> test_indicator,
						                    string genotypes_str,  
						                    vector <INFO> test_info_s_block_full, 
						                    vector <INFO> test_info_l_block_full, 
						                    int test_num_s_block, 
						                    int test_num_l_block, 
						                    bool badsnp_s[],  
						                    bool badsnp_l[]){
	SNPPROC cSP;
	IO cIO; 
	ifstream bed_in(bed_str.c_str(), ios::binary);
	//open stream for observed data (not the reference panel)
	ifstream dat_in(genotypes_str.c_str(), ios::binary);
	
	// INFO small effect SNPs 
	// vector <INFO*> info_s_block(num_s_block);
	// for (int i = 0; i < num_s_block; i++)
		// info_s_block[i] = info_s_block_full[i];
	vector <INFO> info_s_block(num_s_block);
	vector <INFO> test_info_s_block(test_num_s_block);
	for (int i = 0; i < num_s_block; i++)
	  info_s_block[i] = info_s_block_full[i];
	for (int i = 0; i < test_num_s_block; i++)
	  test_info_s_block[i] = test_info_s_block_full[i]; 
	// z_s
	vec z_s = zeros<vec>(num_s_block); 
	for (int i = 0; i < num_s_block; i++) 
		// z_s(i) = info_s_block[i]->z;
		z_s(i) = info_s_block[i].z;
	// small effect genotype matrix
	mat geno_s = zeros<mat>(n_ref, num_s_block); // from reference data
	cout << "number of rows in geno_s: " << geno_s.n_rows << endl;
	cout << "number of columns in geno_s: " << geno_s.n_cols << endl;
	

	//make a test_indicator indicator vector
	cout << "length of test_indicator: " << test_indicator.size() << endl;
	unsigned int n_test = std::accumulate(test_indicator.begin(), test_indicator.end(),
                                       decltype(test_indicator)::value_type(0)); //https://stackoverflow.com/questions/3221812/how-to-sum-up-elements-of-a-c-vector
	//https://stackoverflow.com/questions/3221812/how-to-sum-up-elements-of-a-c-vector
	//initialize a matrix for reading test genotype data
	arma::mat X_s = zeros<mat>(n_test, num_s_block);
	cout << "Number of columns in X_s: " << X_s.n_cols << endl;
	cout << "Number of rows in X_s: " << X_s.n_rows << endl;
  int second_index = 0;
	for (int i = 0; i < num_s_block; ++i) {
	  vec geno = zeros<vec>(n_ref);
	  arma::vec gg = zeros <vec>(n_test);
	  double maf = 0.0; 
	  // cIO.readSNPIm(info_s_block[i]->pos, n_ref, idv, bed_in, geno, maf);
	  cIO.readSNPIm(info_s_block[i].pos, n_ref, idv, bed_in, geno, maf);
	  cSP.nomalizeVec(geno);
	  geno_s.col(i) = geno;
	  //		int pp_ref = info_s_block[i].pos;
	  string ss_ref = info_s_block[i].snp;
	  //check that ss_ref and pp_ref equal test_info_s_block[j].snp and test_info_s_block[j].pos for some j >= i
	  
	  while (badsnp_s[second_index] == true){
	    ++second_index;
	  }
	  cIO.readSNPIm(test_info_s_block[second_index].pos, n_test, test_indicator, dat_in, gg, maf);
	  cSP.nomalizeVec(gg);
	  X_s.col(i) = gg;
	  ++second_index;
	}
	// pseudo EFF
	EFF eff_pseudo; 
	eff_pseudo.snp = "rs"; 
	eff_pseudo.a1 = "Y"; 
	eff_pseudo.maf = 0.0; 
	eff_pseudo.beta = 0.0; 
	arma::mat result;
	vec beta_s = zeros<vec>(num_s_block); 
	// INFO large effect SNPs 
	if (num_l_block != 0){
		// vector <INFO*> info_l_block(num_l_block);
		vector <INFO> info_l_block(num_l_block);
		vector <INFO> test_info_l_block(test_num_l_block);
		for (int i = 0; i < num_l_block; i++) 
			info_l_block[i] = info_l_block_full[i];
		for (int i = 0; i < test_num_l_block; i++) 
		  test_info_l_block[i] = test_info_l_block_full[i];
		// z_l
		vec z_l = zeros<vec>(num_l_block); 
		for(int i = 0; i < num_l_block; i++) 
			// z_l(i) = info_l_block[i]->z;
			z_l(i) = info_l_block[i].z;
		// large effect matrix
		mat geno_l = zeros<mat>(n_ref, num_l_block);
		arma::mat X_l = zeros<mat>(n_test, num_l_block);
		int second_index = 0;
		for (int i = 0; i < num_l_block; ++i) {//num_l_block is the number of large effect SNPs in the block
			vec geno = zeros<vec>(n_ref);
		  arma::vec gg = zeros<vec>(n_test);
			double maf = 0.0; 
			// cIO.readSNPIm(info_l_block[i]->pos, n_ref, idv, bed_in, geno, maf);
			cIO.readSNPIm(info_l_block[i].pos, n_ref, idv, bed_in, geno, maf);
			cSP.nomalizeVec(geno);
			geno_l.col(i) = geno;
			string ss_ref = info_l_block[i].snp;
			//check that ss_ref and pp_ref equal test_info_s_block[j].snp and test_info_s_block[j].pos for some j >= i
			while (badsnp_l[second_index] == true){
			  ++second_index;
			}
			cIO.readSNPIm(test_info_l_block[second_index].pos, n_test, test_indicator, dat_in, gg, maf);
			cSP.nomalizeVec(gg);
			X_l.col(i) = gg;
			++second_index;
		}//end populating of geno_l & of X_l. 
		//X_l is the genotypes data for the non-reference (ie "summary")subjects
		//partition subjects - for X_s and X_l - into training and test
		cout << "dimensions of geno_s: " << geno_s.n_rows << " rows and " << geno_s.n_cols << " columns" << endl;

		// estimation
		vec beta_l = zeros<vec>(num_l_block); 
		arma::field <arma::mat> out = estBlock(n_ref, 
                                           n_obs, 
                                           sigma_s, 
                                           geno_s, //ref data
                                           geno_l, //ref data
                                           z_s, 
                                           z_l, 
                                           beta_s, 
                                           beta_l);
		//variance calcs
		result = calc_nt_by_nt_matrix(out(2), //Sigma_ss 
                                          out(1), //Sigma_sl - no need transpose because estBlock outputs Sigma_sl
                                          out(0), //Sigma_ll
                                          sigma_s, 
                                          n_obs, // SHOULD THIS BE n_test??? 
                                          X_l, 
                                          X_s);
		
		
		// summary 
		for(int i = 0; i < num_l_block; i++) {
			EFF eff_l; 
			// eff_l.snp = info_l_block[i]->snp;
			// eff_l.a1 = info_l_block[i]->a1;
			// eff_l.maf = info_l_block[i]->maf;
			eff_l.snp = info_l_block[i].snp;
			eff_l.a1 = info_l_block[i].a1;
			eff_l.maf = info_l_block[i].maf;
			eff_l.beta = beta_l(i);
			eff_l_block[i] = eff_l;
		}
	}
	else{
	  // estimation
		
		arma::field <arma::mat> out = estBlock(n_ref, 
                                           n_obs, //size of training set
                                           sigma_s, 
                                           geno_s, 
                                           z_s, 
                                           beta_s);
		//variance calcs
		result = calc_nt_by_nt_matrix(out(0), //Sigma_ss 
                                          sigma_s, 
                                          n_obs, // SHOULD THIS BE n_test???
                                          X_s);
		
		eff_l_block[0].snp = eff_pseudo.snp;
		eff_l_block[0].a1 = eff_pseudo.a1;
		eff_l_block[0].maf = eff_pseudo.maf;
		eff_l_block[0].beta = eff_pseudo.beta;
	}
	// output small effect
	for(int i = 0; i < num_s_block; i++) {
		EFF eff_s; 
		// eff_s.snp = info_s_block[i]->snp;
		// eff_s.a1 = info_s_block[i]->a1;
		// eff_s.maf = info_s_block[i]->maf;
		eff_s.snp = info_s_block[i].snp;
		eff_s.a1 = info_s_block[i].a1;
		eff_s.maf = info_s_block[i].maf;
		eff_s.beta = beta_s(i); 
		eff_s_block[i] = eff_s;
	}
	return result.diag(); 
}

// estimate only small effect for each block
arma::vec DBSLMMFIT::calcBlock(int n_ref, 
                               int n_obs, 
                               double sigma_s, 
                               vector<int> idv, //indicator for missingness in reference data - typically this is a vector of all zeroes to indicate no missingness, since we don't care about reference data's trait values
                               string bed_str, //for reference data
                               vector <INFO> info_s_block_full, 
                               int num_s_block, 
                    					 vector <EFF> &eff_s_block,
                    					 vector <int> test_indicator,
                    					 string genotypes_str, 
                    					 vector <INFO> test_info_s_block_full,
                    					 int test_num_s_block,
                    					 bool badsnp_s[]){
	SNPPROC cSP;
	IO cIO; 
	ifstream bed_in(bed_str.c_str(), ios::binary);
	//open stream for observed data (not the reference panel)
	ifstream dat_in(genotypes_str.c_str(), ios::binary);
	// INFO small effect SNPs 
	// vector <INFO*> info_s_block(num_s_block);
	// for (int i = 0; i < num_s_block; i++)
		// info_s_block[i] = info_s_block_full[i];
	vector <INFO> info_s_block(num_s_block);
	vector <INFO> test_info_s_block(test_num_s_block);
	for (int i = 0; i < num_s_block; i++)
		info_s_block[i] = info_s_block_full[i];
	for (int i = 0; i < test_num_s_block; i++)
	  test_info_s_block[i] = test_info_s_block_full[i];
	// z_s
	vec z_s = zeros<vec>(num_s_block); 
	for (int i = 0; i < num_s_block; i++) 
		// z_s(i) = info_s_block[i]->z;
		z_s(i) = info_s_block[i].z;
	cout << "length of test_indicator: " << test_indicator.size() << endl;
	
	// 337129 because that is the number of UKB subjects in the fam file
	unsigned int n_test = std::accumulate(test_indicator.begin(), test_indicator.end(),
                                       decltype(test_indicator)::value_type(0)); //https://stackoverflow.com/questions/3221812/how-to-sum-up-elements-of-a-c-vector
	// small effect genotype matrix
	mat geno_s = zeros<mat>(n_ref, num_s_block);
	arma::mat X_s = zeros <mat>(n_test, num_s_block);
	int second_index = 0;
	for (int i = 0; i < num_s_block; ++i) {
		vec geno = zeros<vec>(n_ref);
	  arma::vec gg = zeros <vec>(n_test);
		double maf = 0.0; 
		// cIO.readSNPIm(info_s_block[i]->pos, n_ref, idv, bed_in, geno, maf);
		cIO.readSNPIm(info_s_block[i].pos, n_ref, idv, bed_in, geno, maf);
		cSP.nomalizeVec(geno);
		geno_s.col(i) = geno;
//		int pp_ref = info_s_block[i].pos;
		string ss_ref = info_s_block[i].snp;
		//check that ss_ref and pp_ref equal test_info_s_block[j].snp and test_info_s_block[j].pos for some j >= i

		while (badsnp_s[second_index] == true){
		  ++second_index;
		}
		cIO.readSNPIm(test_info_s_block[second_index].pos, n_test, test_indicator, dat_in, gg, maf);
		cSP.nomalizeVec(gg);
		X_s.col(i) = gg;
		++second_index;
	}
	
	// estimation
	vec beta_s = zeros<vec>(num_s_block); 
	// partition subjects into training and test sets

	//call estBlock on training data
	arma::field <arma::mat> out = estBlock(n_ref, 
                                        n_obs, 
                                        sigma_s, 
                                        geno_s, 
                                        z_s, 
                                        beta_s);
  //variance calcs
  arma::mat result = calc_nt_by_nt_matrix(out(0), 
                                          sigma_s, 
                                          n_obs, //// SHOULD THIS BE n_test??? 
                                          X_s);
	
	// output small effect
	for(int i = 0; i < num_s_block; i++) {
		EFF eff_s; 
		// eff_s.snp = info_s_block[i]->snp;
		// eff_s.a1 = info_s_block[i]->a1;
		// eff_s.maf = info_s_block[i]->maf;
		eff_s.snp = info_s_block[i].snp;
		eff_s.a1 = info_s_block[i].a1;
		eff_s.maf = info_s_block[i].maf;
		eff_s.beta = beta_s(i); 
		eff_s_block[i] = eff_s;
	}
	return result.diag(); 
}

// solve the equation Ax=b, x is a variables
vec DBSLMMFIT::PCGv(mat A, vec b, size_t maxiter, const double tol){
	vec dA = A.diag();
	// stable checking, using available func to speed up
	for(size_t i = 0; i< dA.n_elem; i++){
		if(dA[i] == 0)
			dA[i] = 1e-4;
	}// end for i
	vec Minv = 1.0/dA;
	// initialize
	vec x = zeros<vec>(b.n_elem);
	vec r = zeros<vec>(b.n_elem);
	vec r1 = zeros<vec>(b.n_elem);
	vec z1 = zeros<vec>(b.n_elem);
	r = b;
	vec z = Minv % r;
	vec p = z;
	size_t iter = 0;
	double sumr2 = norm(r, 2);
	// PCG main loop 
	while( sumr2 > tol && iter < maxiter){
		iter += 1;
		// move direction
		vec Ap = A*p;
		// step size
		double a = dot(r,z)/dot(p,Ap);
		// move
		x = x + a * p;
		r1 = r - a * Ap;
		z1 = Minv % r1;
		double bet = dot(z1, r1)/dot(z, r);
		p = z1 + bet * p;
		z = z1;
		r = r1;
		sumr2 = norm(r, 2);
	}// end while loop
	if (iter >= maxiter){
		cerr << "ERROR: Matrix is Singular!" << endl;
	}
	return(x);
}// end function

mat DBSLMMFIT::PCGm(mat A, mat B, size_t maxiter, const double tol){

	size_t n_iter = B.n_cols;
	mat x = zeros<mat>(A.n_rows, n_iter);
	for (size_t i = 0; i < n_iter; i++){
		x.col(i) = PCGv(A, B.col(i), maxiter, tol);
	}// end for loop
	return(x);
}// end function

arma::field< arma::mat > DBSLMMFIT::estBlock(int n_ref, 
                                             int n_obs, 
                                             double sigma_s, 
                                             mat geno_s, 
                                             mat geno_l, 
                                             vec z_s, 
                                             vec z_l, 
                                             vec &beta_s, 
                                             vec &beta_l) {
	// LD matrix 
	// mat SIGMA_ls = geno_l.t() * geno_s; 
	// SIGMA_ls /= (double)n_ref; 
	// mat SIGMA_ll = geno_l.t() * geno_l; 
	// SIGMA_ll /= (double)n_ref;
	// mat SIGMA_ss = geno_s.t() * geno_s; 
	// SIGMA_ss /= (double)n_ref;
	
	double tau = 0.8;
	mat SIGMA_ls = geno_l.t() * geno_s;
	SIGMA_ls *= tau/(double)n_ref;
	mat SIGMA_ll = geno_l.t() * geno_l; 
	SIGMA_ll *= tau/(double)n_ref;
	vec DIAG_l(geno_l.n_cols, fill::ones);
	DIAG_l *= (1.0-tau);
	SIGMA_ll += diagmat(DIAG_l);
	mat SIGMA_ss = geno_s.t() * geno_s; 
	SIGMA_ss *= tau/(double)n_ref;
	vec DIAG_s(geno_s.n_cols, fill::ones);
	DIAG_s *= (1.0-tau);
	SIGMA_ss += diagmat(DIAG_s);

	// beta_l
	SIGMA_ss.diag() += 1.0 / (sigma_s * (double)n_obs);
	mat SIGMA_ss_inv_SIGMA_sl = PCGm(SIGMA_ss, SIGMA_ls.t(), 1000, 1e-7);
	mat SIGMA_ls_SIGMA_ss_inv_SIGMA_sl = - SIGMA_ls * SIGMA_ss_inv_SIGMA_sl;
	SIGMA_ls_SIGMA_ss_inv_SIGMA_sl += SIGMA_ll;
	vec SIGMA_ss_inv_z_s = PCGv(SIGMA_ss, z_s, 1000, 1e-7);
	vec SIGMA_ls_SIGMA_ss_inv_z_s = -SIGMA_ls * SIGMA_ss_inv_z_s;
	SIGMA_ls_SIGMA_ss_inv_z_s += z_l;
	beta_l = PCGv(SIGMA_ls_SIGMA_ss_inv_SIGMA_sl, SIGMA_ls_SIGMA_ss_inv_z_s, 1000, 1e-7);
	beta_l /= sqrt(n_obs);

	// beta_s
	SIGMA_ss_inv_z_s *= sqrt(n_obs);
	vec SIGMA_ss_inv_SIGMA_sl_beta_l = (double)n_obs * SIGMA_ss_inv_SIGMA_sl * beta_l; 
	vec SIGMA_ss_inv_z_s_SIGMA_sl_beta_l = SIGMA_ss_inv_z_s - SIGMA_ss_inv_SIGMA_sl_beta_l; 
	SIGMA_ss.diag() -= 1.0 / (sigma_s * (double)n_obs);
	vec SIGMA_ss_z_s_SIGMA_sl_beta_l = SIGMA_ss * SIGMA_ss_inv_z_s_SIGMA_sl_beta_l; 
	beta_s = sqrt(n_obs) * z_s - (double)n_obs * SIGMA_ls.t() * beta_l - SIGMA_ss_z_s_SIGMA_sl_beta_l; 
	beta_s *= sigma_s;
	arma::field<arma::mat> result(3);
	result(0) = SIGMA_ss;
	result(1) = arma::trans(SIGMA_ls);//ie, Sigma_sl
//	cout << "number of rows in SIGMA_ls: " << SIGMA_ls.n_rows << endl;
//	cout << "number of columns in SIGMA_ls: " << SIGMA_ls.n_cols << endl; 
	result(2) = SIGMA_ll;
	
	return result; 
}

arma::field <arma::mat> DBSLMMFIT::estBlock(int n_ref, 
                                            int n_obs, 
                                            double sigma_s, 
                                            mat geno_s, 
                                            vec z_s, 
                                            vec &beta_s) {

	// LD matrix 
	// mat SIGMA_ss = geno_s.t() * geno_s; 
	// SIGMA_ss /= (double)n_ref;

	double tau = 0.8;
	mat SIGMA_ss = geno_s.t() * geno_s; 
	SIGMA_ss *= tau/(double)n_ref;
	vec DIAG_s(geno_s.n_cols, fill::ones);
	DIAG_s *= (1.0-tau);
	SIGMA_ss += diagmat(DIAG_s);
	
	// beta_s
	SIGMA_ss.diag() += 1.0 / (sigma_s * (double)n_obs);
	vec SIGMA_ss_inv_z_s = PCGv(SIGMA_ss, z_s, 1000, 1e-7);
	SIGMA_ss.diag() -= 1.0 / (sigma_s * (double)n_obs);
	vec SIGMA_ss_SIGMA_ss_inv_z_s = SIGMA_ss * SIGMA_ss_inv_z_s;
	vec z_s_SIGMA_ss_SIGMA_ss_inv_SIGMA_sl = z_s - SIGMA_ss_SIGMA_ss_inv_z_s; 
	beta_s = sqrt(n_obs) * sigma_s * z_s_SIGMA_ss_SIGMA_ss_inv_SIGMA_sl; 
	
	arma::field <arma::mat> result(3);//length is 3 for compatibility with overloading
	result(0) = SIGMA_ss;

	return result; 
}
