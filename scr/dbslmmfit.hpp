/*
Deterministic version Bayesian Sparse Linear Mixed Model (DBSLMM)
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

#ifndef __DSBSLMMFIT_H__
#define __DSBSLMMFIT_H__

#include <iostream>
#include <vector>
#include <string>
#include <armadillo>
#include <numeric>
#include "dtpr.hpp"
#include "calc_asymptotic_variance.hpp"
#include "subset_to_test_and_training.hpp"


using namespace std;
using namespace arma;

class DBSLMMFIT {
public:
	// estimate large and small effect
	int est(int n_ref, 
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
         vector <INFO> test_info_l 
	);
	// estimate only small effect
	int est(int n_ref, 
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
         vector <INFO> test_info_s
	);
	// estimate large and small effect for each block
	arma::vec calcBlock(int n_ref, 
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
                     vector <INFO> test_info_l_block_full
	);
	// estimate only small effect for each block
	arma::vec calcBlock(int n_ref, 
                     int n_obs, 
                     double sigma_s, 
                     vector<int> idv, //indicator for missingness in reference data - typically this is a vector of all zeroes to indicate no missingness, since we don't care about reference data's trait values
                     string bed_str, //for reference data
                     vector <INFO> info_s_block_full, 
                     int num_s_block, 
                     vector <EFF> &eff_s_block,
                     vector <int> test_indicator,
                     string genotypes_str, 
                     vector <INFO> test_info_s_block_full
	);
	// solve x=Ab
	vec PCGv(mat A, vec b, size_t maxiter, const double tol); 
	// solve x=AB
	mat PCGm(mat A, mat B, size_t maxiter, const double tol);
	// small and large effect
	arma::field <arma::mat> estBlock(int n_ref, int n_obs, double sigma_s, mat geno_s, mat geno_l, vec z_s, vec z_l, vec &beta_s, vec &beta_l);
	// only small effect
	arma::field <arma::mat> estBlock(int n_ref, int n_obs, double sigma_s, mat geno_s, vec z_s, vec &beta_s);
};
#endif
