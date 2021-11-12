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
#include <tuple> // std::tuple, std::get, std::tie, std::ignore

#include "omp.h"

#include "dtpr.hpp"
#include "dbslmmfit.hpp"
#include "calc_asymptotic_variance.h"
#include "subset_to_test_and_training.h"

using namespace std;
using namespace arma;

//' Estimate large and small effects
//' 
//' @param n_ref sample size of the reference panel
//' @param n_obs sample size of the data
//' @param sigma_s the estimate for sigma_s^2
//' @param num_block number of blocks in the genome
//' @param idv 
//' @param bed_str file path for the plink bed file
//' @param info_s small effect SNP info object
//' @param info_l large effect SNP info object
//' @param thread number of threads to use
//' @param eff_s small effects SNP effects object
//' @param eff_l large effects SNP effects object
//' @param fam_file path to the plink fam file
//' @param seed a seed (positive integer) for the pseudo RNG
//' @return zero is returned
// estimate large and small effect
int DBSLMMFIT::est(int n_ref, 
                   int n_obs, 
                   double sigma_s, 
                   int num_block, 
                   vector<int> idv, 
                   string bed_str,
                   vector <INFO> info_s, // one entry per small effect SNP
                   vector <INFO> info_l, //one entry per large effect SNP
                   int thread, 
                   vector <EFF> &eff_s, 
                   vector <EFF> &eff_l){ 
  //return to Sheng's code
	// get the maximum number of each block
	int count_s = 0, count_l = 0; //set counters at zero
	vec num_s = zeros<vec>(num_block), num_l = zeros<vec>(num_block); //num_s & num_l have one entry per block 
	for (int i = 0; i < num_block; i++) {
		for (size_t j = count_s; j < info_s.size(); j++) {
			if(info_s[j].block == i){ //within info_s, block is an integer. I suspect that it indicates block membership for a SNP, ie, if a SNP is in block 10, then this value is 10.
				num_s(i) += 1; //num_s is a vector with one entry per block. It will contain the number of small effect SNPs in each block.
				count_s++; //count_s becomes the number of small effect SNPs in the whole genome
			}else{
				break;
			}
		}//loop above populates num_s. Each entry of num_s is the number of small effect SNPs in a block
		for (size_t j = count_l; j < info_l.size(); j++) {
			if(info_l[j].block == i){ //again, I think this indicates block membership
				num_l(i) += 1;  //num_l becomes the vector containing per-block counts of large effect SNPs
				count_l++; //count_l becomes the number of genome-wide large effect SNPs
			}else{
				break;
			}
		}
	}// end of counting and populating num_l and num_s, ie, iterates over i
	count_l = count_s = 0; // reset
	
	double len_l = num_l.max(); //len_l is the maximum, across blocks, of the per-block number of large effect SNPs
	double len_s = num_s.max(); //len_s is the max of the per-block number of small effect SNPs
	
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
	vector < vector <INFO> > info_s_Block(B_MAX, 
                                       vector <INFO> ((int)len_s)), 
                           info_l_Block(B_MAX, 
                                        vector <INFO> ((int)len_l)); //declare info_l_Block & info_s_Block
	vector < vector <EFF> > eff_s_Block(B_MAX, 
                                     vector <EFF> ((int)len_s)), 
                          eff_l_Block(B_MAX, 
                                      vector <EFF> ((int)len_l)); //declare eff_l_Block & eff_s_Block
	vector <int> num_s_vec, num_l_vec; //declare num_s_vec & num_l_vec
	for (int i = 0; i < num_block; ++i) {//iterate over blocks, ie, i indexes block number
		// small effect SNP information
		vector <INFO> info_s_block; //declare info_s_block
		for (size_t j = count_s; j < info_s.size(); j++) {
			if(info_s[j].block == i){ 
				info_s_block.push_back(info_s[j]);//appends info_s[j] to info_s_block, ie, does the work to create INFO for the snps in the block of interest, ie because the snp has snpinfo block value matching
			  //https://www.cplusplus.com/reference/vector/vector/push_back/
				count_s++; //count_s records the number of small SNPs in the genome. it starts at zero per line above
				
			}else{
				break;
			}
		}
		for (size_t k = 0; k < info_s_block.size(); k++)
			// info_s_Block[B][k] = &info_s_block[k]; 
			info_s_Block[B][k] = info_s_block[k]; //Organizes the contents of info_s_block into the larger info_s_Block. info_s_Block will hold info for all blocks, it seems
		num_s_vec.push_back((int)num_s(i)); //num_s was determined in lines 64-67. it's the vector with length equal to the number of blocks and entries equal to the number of small effects snps per block.
		
		// large effect SNP information
		if (num_l(i) == 0){
			// info_l_Block[B][0] = &info_pseudo;
			info_l_Block[B][0] = info_pseudo; //put the pseudo info object in for those blocks without large effect snps.
		}else{ //when there is one or more large effect snp in the block
			vector <INFO> info_l_block;//declare info_l_block
			for (size_t j = count_l; j < info_l.size(); j++) {//count_l keeps a tally of genomewide number of large effect SNPs
				if(info_l[j].block == i){ 
					info_l_block.push_back(info_l[j]);
					count_l++;
				}else{
					break;
				}
			}
			for (size_t k = 0; k < info_l_block.size(); k++)
				// info_l_Block[B][k] = &info_l_block[k]; 
				info_l_Block[B][k] = info_l_block[k]; //info_l_Block is like a vector of vectors. B, the first index, indexes the blocks, while k indexes the markers within a block
		}
		num_l_vec.push_back((int)num_l(i));

		B++; //increment B 
		if (B == B_MAX || i + 1 == num_block) { // process the block of SNPs using multi-threading
			
			omp_set_num_threads(thread);
#pragma omp parallel for schedule(dynamic)
			for (int b = 0; b < B; b++){
				calcBlock(n_ref, 
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
              y_training, 
              training_indices, 
              test_indices, 
              i);
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
			B = 0;// reset B to zero
			num_l_vec.clear(); 
			num_s_vec.clear();
		}//end if statement starting on line 154, if (B == B_MAX...
	}//end loop for i
	return 0;
}//end function

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
          				 string fam_file,
          				 unsigned int seed,
          				 double test_proportion){
  // split subjects into training and test sets
  //specify proportion of n_obs that goes into test set
  arma::Col<arma::uword> test_indices = get_test_indices(n_obs, 
                                                         test_proportion,
                                                         seed);
  arma::Col<arma::uword> training_indices = get_training_indices(test_indices, 
                                                                 n_obs);
  // read phenotype data
  std::tuple<vector<string>, vector<string> > pheno_struct = read_pheno(fam_file, 6);
  // extract id and pheno from pheno_struct
  std::vector<string> id = std::get<0>(pheno_struct);
  std::vector<string> pheno_string = std::get<1>(pheno_struct);
  //convert pheno to numeric vector
  std::vector<double> pheno_numeric = convert_string_vector_to_double_vector(pheno_string);
  // convert to arma::vec
  arma::vec pheno_arma = arma::conv_to<arma::vec>::from(pheno_numeric);
  //mean center pheno
  arma::vec y = center_vector(pheno_arma);
  arma::vec y_training = subset(y, training_indices);
  arma::vec y_test = subset(y, test_indices);
  //save y_test as csv
//  y_test.save("y_test.csv", arma_ascii);
  //return to Sheng's code
  
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
	vector < vector <EFF> > eff_s_Block(B_MAX, vector <EFF> ((int)len_s));
	vector <int> num_s_vec;
	for (int i = 0; i < num_block; ++i) {
		// small effect SNP information
		vector <INFO> info_s_block; // declare object for small effect SNP info
		for (size_t j = count_s; j < info_s.size(); j++) { //info_s.size is the number of SNPs in info_s
			if(info_s[j].block == i){ //if jth SNP in info_s is in block i
				info_s_block.push_back(info_s[j]);
				count_s++; //increase count_s by 1, ie, count of number of small effect SNPs increases by 1
			}else{ //for SNPs not in block i, just leave the loop
				break;
			}
		}
		for (size_t k = 0; k < info_s_block.size(); k++)
			// info_s_Block[B][k] = &info_s_block[k]; 
			info_s_Block[B][k] = info_s_block[k]; 
		num_s_vec.push_back((int)num_s(i));
		
		B++;
		if (B == B_MAX || i + 1 == num_block) { // process the block of SNPs using multi-threading
			
			omp_set_num_threads(thread);
#pragma omp parallel for schedule(dynamic)
			for (int b = 0; b < B; b++){
				calcBlock(n_ref, 
              n_obs, 
              sigma_s, 
              idv, 
              bed_str, 
              info_s_Block[b],
						  num_s_vec[b], 
              eff_s_Block[b], 
              y_training, 
              training_indices, 
              test_indices, 
              i);
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
	return 0;
}



//' Estimate large and small effects for each block
//' 
//' @param n_ref sample size of reference panel
//' @param n_obs sample size of data
//' @param sigma_s estimate of $sigma^2_s$
//' @param idv 
//' @param bed_str filename for bed file
//' @param info_s_block_full info object for small effects
//' @param info_l_block_full info object for large effects
//' @param num_s_block
//' @param num_l_block
//' @param eff_s_block effects object for small effects per block? 
//' @param eff_l_block effects object for large effects per block?
// estimate large and small effect for each block
int DBSLMMFIT::calcBlock(int n_ref, 
                         int n_obs, 
                         double sigma_s, 
                         vector<int> idv, 
                         string bed_str, 
                         vector <INFO> info_s_block_full, //small
                         vector <INFO> info_l_block_full, //large
                         int num_s_block, 
                         int num_l_block, 
                         vector <EFF> &eff_s_block, 
                         vector <EFF> &eff_l_block,
                         arma::vec y_training,
                         arma::Col<arma::uword> training_indices, 
                         arma::Col<arma::uword> test_indices, 
                         int iter_number){
	SNPPROC cSP;
	IO cIO; 
	ifstream bed_in(bed_str.c_str(), ios::binary);
	
	// INFO small effect SNPs 
	// vector <INFO*> info_s_block(num_s_block);
	// for (int i = 0; i < num_s_block; i++)
		// info_s_block[i] = info_s_block_full[i];
	vector <INFO> info_s_block(num_s_block);
	for (int i = 0; i < num_s_block; i++)
		info_s_block[i] = info_s_block_full[i];
	// z_s
	vec z_s = zeros<vec>(num_s_block); 
	for (int i = 0; i < num_s_block; i++) 
		// z_s(i) = info_s_block[i]->z;
		z_s(i) = info_s_block[i].z;
	// small effect genotype matrix
	mat geno_s = zeros<mat>(n_ref, num_s_block); //geno_s gets filled to become the matrix Xs
	for (int i = 0; i < num_s_block; ++i) {
		vec geno = zeros<vec>(n_ref);
		double maf = 0.0; 
		// cIO.readSNPIm(info_s_block[i]->pos, n_ref, idv, bed_in, geno, maf);
		cIO.readSNPIm(info_s_block[i].pos, n_ref, idv, bed_in, geno, maf);
		cSP.nomalizeVec(geno);
		geno_s.col(i) = geno;
	}
	// pseudo EFF
	EFF eff_pseudo; 
	eff_pseudo.snp = "rs"; 
	eff_pseudo.a1 = "Y"; 
	eff_pseudo.maf = 0.0; 
	eff_pseudo.beta = 0.0; 
	
	vec beta_s = zeros<vec>(num_s_block); 
	// INFO large effect SNPs 
	if (num_l_block != 0){
		// vector <INFO*> info_l_block(num_l_block);
		vector <INFO> info_l_block(num_l_block);
		for (int i = 0; i < num_l_block; i++) 
			info_l_block[i] = info_l_block_full[i];
		// z_l
		vec z_l = zeros<vec>(num_l_block); 
		for(int i = 0; i < num_l_block; i++) 
			// z_l(i) = info_l_block[i]->z;
			z_l(i) = info_l_block[i].z;


		// large effect matrix
		mat geno_l = zeros<mat>(n_ref, num_l_block); //geno_l gets filled to become Xl matrix (in notation of supplemental materials)
		for (int i = 0; i < num_l_block; ++i) {
			vec geno = zeros<vec>(n_ref);
			double maf = 0.0; 
			// cIO.readSNPIm(info_l_block[i]->pos, n_ref, idv, bed_in, geno, maf);
			cIO.readSNPIm(info_l_block[i].pos, n_ref, idv, bed_in, geno, maf);
			cSP.nomalizeVec(geno);
			geno_l.col(i) = geno;
		}
		/* INSERT MY ASYMPTOTIC VAR CALC HERE*/
		//split geno_l and geno_s into training and test sets
		arma::mat geno_l_training = subset(geno_l, training_indices);
		arma::mat geno_l_test = subset(geno_l, test_indices);
		arma::mat geno_s_training = subset(geno_s, training_indices);
		arma::mat geno_s_test = subset(geno_s, test_indices);
		
		// calculate var(\hat\tilde y)
		arma::mat asymptotic_var = calc_asymptotic_variance(geno_l_training, 
                                                       geno_s_training, 
                                                       geno_l_test,
                                                       geno_s_test,
                                                       sigma_s,
                                                       y_training);
    // asymptotic_var should be n_test by n_test symmetric psd matrix, ie covar matrix
    //store only diagonal elements of asymptotic_var, asymptotic_var.diag()
    arma::vec avar_diag = asymptotic_var.diag();
    // define outfile
    std::string iter_number_string = to_string(iter_number);
    std::string outfile = iter_number_string + ".csv";
    //save diagonal as csv
    avar_diag.save(outfile, arma_ascii); 
    
		/* END OF FREDS ASYMPTOTIC VAR CALC CODE */
		// estimation
		vec beta_l = zeros<vec>(num_l_block); 
		estBlock(n_ref, n_obs, sigma_s, geno_s, geno_l, z_s, z_l, beta_s, beta_l);//estBlock!
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
	else{ //case of num_l_block ==0
		// estimation
		estBlock(n_ref, n_obs, sigma_s, geno_s, z_s, beta_s); // estBlock!
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
	return 0; 
}

// estimate only small effect for each block
int DBSLMMFIT::calcBlock(int n_ref, 
                         int n_obs, 
                         double sigma_s, 
                         vector<int> idv, 
                         string bed_str, 
              					 vector <INFO> info_s_block_full, 
            						 int num_s_block, 
            						 vector <EFF> &eff_s_block,
            						 arma::vec y_training,
            						 arma::Col<arma::uword> training_indices, 
            						 arma::Col<arma::uword> test_indices, 
            						 int iter_number
            						 ){
	SNPPROC cSP; // declare new SNPPROC object, cSP. Below, we'll need to populate cSP.
	IO cIO; //declare IO object, cIO
	ifstream bed_in(bed_str.c_str(), ios::binary);//ios::binary means "open in binary mode". bed_str is an argument to the function, presumably something like the file path??
	//https://www.cplusplus.com/reference/fstream/ifstream/. ifstream is a class in std, ie, std::ifstream. Used to read from files
	//https://www.cplusplus.com/doc/tutorial/files/
	//https://www.cplusplus.com/reference/string/string/c_str/. Get C string equivalent
	//Returns a pointer to an array that contains a null-terminated sequence of characters (i.e., a C-string) representing the current value of the string object.
	// INFO small effect SNPs 
	// vector <INFO*> info_s_block(num_s_block);
	// for (int i = 0; i < num_s_block; i++)
		// info_s_block[i] = info_s_block_full[i];
	vector <INFO> info_s_block(num_s_block);
	for (int i = 0; i < num_s_block; i++)
		info_s_block[i] = info_s_block_full[i];//copy info_s_block_full to info_s_block
	// z_s
	vec z_s = zeros<vec>(num_s_block); //init. z_s
	for (int i = 0; i < num_s_block; i++) 
		// z_s(i) = info_s_block[i]->z;
		z_s(i) = info_s_block[i].z;//populate z_s with the z values from the INFO object
	// small effect genotype matrix
	mat geno_s = zeros<mat>(n_ref, num_s_block);// Is num_s_block the number of blocks with small effects only? Or something else???
	// num_s_block must be the number of small effect SNPs in the block. n_ref and num_s_block are dimensions for the matrix
	for (int i = 0; i < num_s_block; ++i) {
		vec geno = zeros<vec>(n_ref); //init. geno as arma::vec of length n_ref
		double maf = 0.0; //initialize maf.
		// cIO.readSNPIm(info_s_block[i]->pos, n_ref, idv, bed_in, geno, maf);
		cIO.readSNPIm(info_s_block[i].pos, n_ref, idv, bed_in, geno, maf);// this line populates geno
		cSP.nomalizeVec(geno); // then, normalize geno
		geno_s.col(i) = geno; //write geno to a column of geno_s matrix
	}
	/* INSERT ASYMPTOTIC VARIANCE CALCS HERE */
	//split geno_l and geno_s into training and test sets
	arma::mat geno_s_training = subset(geno_s, training_indices);
	arma::mat geno_s_test = subset(geno_s, test_indices);
	
	// calculate var(\hat\tilde y)
/*	arma::mat asymptotic_var = calc_asymptotic_variance(geno_l_training, 
                                                     geno_s_training, 
                                                     geno_l_test,
                                                     geno_s_test,
                                                     sigma_s,
                                                     y_training);

	// asymptotic_var should be n_test by n_test symmetric psd matrix, ie covar matrix
	//store only diagonal elements of asymptotic_var, asymptotic_var.diag()
	arma::vec avar_diag = asymptotic_var.diag();
	// define outfile
	std::string iter_number_string = to_string(iter_number);
	std::string outfile = iter_number_string + ".csv";
	//save diagonal as csv
	avar_diag.save(outfile, arma_ascii); 
NEED METHODS FOR calculating asymptotic var when a block has no large effects */ 


	
	
	/* END ASYMPTOTIC VAR CALCS */
	// estimation
	vec beta_s = zeros<vec>(num_s_block); //num_s_block, ie, the number of small effect SNPs in the block, is the length of beta_s, ie, since every small effect SNP will be represented by one entry in beta_s
	estBlock(n_ref, n_obs, sigma_s, geno_s, z_s, beta_s);
	
	// output small effect
	for(int i = 0; i < num_s_block; i++) {
		EFF eff_s; //declare EFF object. Overwrites previous eff_s object as we loop over i. That is, eff_s will be populated but it's overwritten for every value of i.
		// eff_s.snp = info_s_block[i]->snp;
		// eff_s.a1 = info_s_block[i]->a1;
		// eff_s.maf = info_s_block[i]->maf;
		eff_s.snp = info_s_block[i].snp;
		eff_s.a1 = info_s_block[i].a1;
		eff_s.maf = info_s_block[i].maf;
		eff_s.beta = beta_s(i); 
		eff_s_block[i] = eff_s;//store eff_s for future use as entry in eff_s_block
	}
	return 0; 
}

// solve the equation Ax=b, x is a variables
vec DBSLMMFIT::PCGv(mat A, 
                    vec b, 
                    size_t maxiter, 
                    const double tol){
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

mat DBSLMMFIT::PCGm(mat A, mat B, size_t maxiter, const double tol){//like PCGv but with matrix B
	
	size_t n_iter = B.n_cols;
	mat x = zeros<mat>(A.n_rows, n_iter);
	for (size_t i = 0; i < n_iter; i++){
		x.col(i) = PCGv(A, B.col(i), maxiter, tol);
	}// end for loop
	return(x);
}// end function


//' Calculate coefficient estimates for one block of SNPs, with both large and small effects
//' 
//' @param n_ref sample size for reference panel
//' @param n_obs sample size for data
//' @param sigma_s estimate of sigma_s^2
//' @param geno_s genotypes matrix for small effect SNPs in this block
//' @param geno_l genotypes matrix for large effect SNPs in this block
//' @param z_s z for small effect SNPs in this block
//' @param z_l z for large effect SNPs in this block
//' @param beta_s coefficient estimates for small effect SNPs in this block
//' @param beta_l coefficient estimates for large effect SNPs in this block 

int DBSLMMFIT::estBlock(int n_ref, 
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
	
	return 0; 
}

int DBSLMMFIT::estBlock(int n_ref, 
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
	
	return 0; 
}
