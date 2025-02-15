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

#include <vector>
#include <string>
#include <iostream>
#include <ctime>
#include <math.h>

#include "dtpr.hpp"
#include "dbslmmfit.hpp"
#include "dbslmm.hpp"
#include "subset_to_test_and_training.hpp"

using namespace std;

DBSLMM::DBSLMM(void) :
	version("0.3"), date("05/01/2021"), year("2021")
{}

void DBSLMM::printHeader(void)
{
	cout << endl;
	cout << "*************************************************************"<< endl;
	cout << "  Deterministic Bayesian Sparse Linear Mixed Model (DBSLMM)  " << endl;
	cout << "  Version " << version << ", " << date << "                  " << endl;
	cout << "  Visit http://www.xzlab.org/software.html For Update        " << endl;
	cout << "  (C) " << year << " Sheng Yang, Xiang Zhou                  " << endl;
	cout << "  GNU General Public License                                 " << endl;
	cout << "  For Help, Type ./dbslmm -h                                 " << endl;
	cout << "*************************************************************" << endl;
	cout << endl;

	return;
}

void DBSLMM::printHelp(void) {
	cout << " FILE I/O RELATED OPTIONS" << endl;
	cout << " -s        [filename]  " << " specify input the summary data for the small effect SNPs." << endl;
	cout << " -l        [filename]  " << " specify input the summary data for the large effect SNPs." << endl;
	cout << " -r        [filename]  " << " specify input the bfile of reference data." << endl;
	cout << " -n        [num]       " << " specify input the sample size of the summary data." << endl;
	cout << " -mafMax   [num]       " << " specify input the maximium of the difference between reference panel and summary data." << endl;
	cout << " -nsnp     [num]  " << " specify input the number of snp." << endl;
	cout << " -b        [num]       " << " specify input the block information." << endl;
	cout << " -h        [num]       " << " specify input the heritability." << endl;
	cout << " -t        [filename]  " << " specify input thread." << endl;
	cout << " -eff      [filename]  " << " specify output the estimate effect SNPs." << endl;
	return;
}

void DBSLMM::Assign(int argc, char ** argv, PARAM &cPar) {
	
	string str;
	for (int i = 0; i < argc; i++) {

		if (strcmp(argv[i], "--smallEff") == 0 || strcmp(argv[i], "-s") == 0) {

			if (argv[i + 1] == NULL || argv[i + 1][0] == '-') { continue; }
			++i;
			str.clear();
			str.assign(argv[i]);
			cPar.s = str;
		}
		else if (strcmp(argv[i], "--largeEff") == 0 || strcmp(argv[i], "-l") == 0) {

			if (argv[i + 1] == NULL || argv[i + 1][0] == '-') { continue; }
			++i;
			str.clear();
			str.assign(argv[i]);
			cPar.l = str;
		}
		else if (strcmp(argv[i], "--reference") == 0 || strcmp(argv[i], "-r") == 0) {

			if (argv[i + 1] == NULL || argv[i + 1][0] == '-') { continue; }
			++i;
			str.clear();
			str.assign(argv[i]);
			cPar.r = str;
		}
		else if (strcmp(argv[i], "--N") == 0 || strcmp(argv[i], "-n") == 0) {

			if (argv[i + 1] == NULL || argv[i + 1][0] == '-') { continue; }
			++i;
			str.clear();
			str.assign(argv[i]);
			cPar.n = atoi(str.c_str());
		}
		else if (strcmp(argv[i], "--mafMax") == 0 || strcmp(argv[i], "-mafMax") == 0) {

			if (argv[i + 1] == NULL || argv[i + 1][0] == '-') { continue; }
			++i;
			str.clear();
			str.assign(argv[i]);
			cPar.mafMax = atof(str.c_str());
		}
		else if (strcmp(argv[i], "--numSNP") == 0 || strcmp(argv[i], "-nsnp") == 0) {

			if (argv[i + 1] == NULL || argv[i + 1][0] == '-') { continue; }
			++i;
			str.clear();
			str.assign(argv[i]);
			cPar.nsnp = atoi(str.c_str());
		}
		else if (strcmp(argv[i], "--block") == 0 || strcmp(argv[i], "-b") == 0) {

			if (argv[i + 1] == NULL || argv[i + 1][0] == '-') { continue; }
			++i;
			str.clear();
			str.assign(argv[i]);
			cPar.b = str;
		}
		else if (strcmp(argv[i], "--Heritability") == 0 || strcmp(argv[i], "-h") == 0) {

			if (argv[i + 1] == NULL || argv[i + 1][0] == '-') { continue; }
			++i;
			str.clear();
			str.assign(argv[i]);
			cPar.h = atof(str.c_str());
		}
		else if (strcmp(argv[i], "--Thread") == 0 || strcmp(argv[i], "-t") == 0) {

			if (argv[i + 1] == NULL || argv[i + 1][0] == '-') { continue; }
			++i;
			str.clear();
			str.assign(argv[i]);
			cPar.t = atoi(str.c_str());
		}
		else if (strcmp(argv[i], "--EFF") == 0 || strcmp(argv[i], "-eff") == 0) {

			if (argv[i + 1] == NULL || argv[i + 1][0] == '-') { continue; }
			++i;
			str.clear();
			str.assign(argv[i]);
			cPar.eff = str;
		}
		else if (strcmp(argv[i], "--test_indicator_file") == 0 || strcmp(argv[i], "-test_indicator_file") == 0) {
		  
		  if (argv[i + 1] == NULL || argv[i + 1][0] == '-') { continue; }
		  ++i;
		  str.clear();
		  str.assign(argv[i]);
		  cPar.test_indicator_file = str.c_str();
		}
		else if (strcmp(argv[i], "--dat_str") == 0 || strcmp(argv[i], "-dat_str") == 0) {
		  
		  if (argv[i + 1] == NULL || argv[i + 1][0] == '-') { continue; }
		  ++i;
		  str.clear();
		  str.assign(argv[i]);
		  cPar.dat_str = str.c_str();
		}

		
	}
	return;
}

void DBSLMM::BatchRun(PARAM &cPar) {

	SNPPROC cSP;
	IO cIO;
	DBSLMMFIT cDBSF;

	// input check
	 cout << "Options: " << endl;
	 cout << "-s:      " << cPar.s << endl;
	 cout << "-l:      " << cPar.l << endl;
	 cout << "-r:      " << cPar.r << endl;
	 cout << "-nsnp:   " << cPar.nsnp << endl;
	 cout << "-n:      " << cPar.n << endl;
	 cout << "-mafMax: " << cPar.mafMax << endl;
	 cout << "-b:      " << cPar.b << endl;
	 cout << "-h:      " << cPar.h << endl;
	 cout << "-t:      " << cPar.t << endl;
	 cout << "-eff:    " << cPar.eff << endl;
	cout << "-test_indicator_file:  " << cPar.test_indicator_file <<endl;

	string ref_fam_str = cPar.r + ".fam"; //next line below declares the ifstream objects, including an ifstream object for reading the fam file! 
	ifstream seffstream(cPar.s.c_str()), leffstream(cPar.l.c_str()), reffstream(ref_fam_str.c_str()), beffstream(cPar.b.c_str());
	if (cPar.s.size() == 0) {
		cerr << "ERROR: -s is no parameter!" << endl;
		exit(1);
	}
	if (!beffstream) {
		cerr << "ERROR: " << cPar.b << " dose not exist!" << endl;
		exit(1);
	}
	if (!seffstream) {
		cerr << "ERROR: " << cPar.s << " dose not exist!" << endl;
		exit(1);
	}
	if (!reffstream) {
		cerr << "ERROR: " << cPar.r << " dose not exist!" << endl;
		exit(1);
	}
	if (cPar.b.size() == 0) {
		cerr << "ERROR: -b is no parameter!" << endl;
		exit(1);
	}
	if (cPar.r.size() == 0) {
		cerr << "ERROR: " << cPar.r << " dose not exist!" << endl;
		exit(1);
	}
	if (cPar.h > 1 || cPar.h < 0){
		cerr << "ERROR: -h is not correct (0, 1)!" << endl;
		exit(1);
	}
	if (cPar.t > 100 || cPar.t < 1){
		cerr << "ERROR: -t is not correct (1, 100)!" << endl;
		exit(1);
	}
	
	// get sample size of reference panel 
	char separate[] = "\t";
	cout << "Reading reference PLINK FAM file from [" << cPar.r << ".fam]" << endl;
	int n_ref = cIO.getRow(ref_fam_str); //n_ref is sample size of reference panel!
	cout << n_ref << " individuals to be included from reference FAM file." << endl;
	
	// get SNP of reference panel
	cout << "Reading reference PLINK BIM file from [" << cPar.r << ".bim]" << endl;
	map <string, ALLELE> ref_bim;
	bool constr = true; 
	if (abs(cPar.mafMax-1.0) < 1e-10){
		constr = false; 
	}
	cIO.readBim(n_ref, cPar.r, separate, ref_bim, constr); //read BIM for ref data
	int num_snp_ref = ref_bim.size(); //num_snp_ref is number of SNPs in the reference data
	cout << num_snp_ref << " SNPs to be included from reference BIM file." << endl;

	// input block file
	vector <BLOCK> block_dat; 
	cIO.readBlock(cPar.b, separate, block_dat); //readBlock is defined in scr/dtpr.cpp. It reads the files that contain the blocking information 
  //cPar.b is the --block option, ie, the file containing the blocking information, ie, the bed file for the reference data. 
  //--block ${BLOCK}.bed, ie, file path to a bed file in, eg., block_data/EUR/
  
	// input small effect summary data
	cout << "Reading summary data of small effect SNPs from [" << cPar.s << "]" << endl;
	vector <SUMM> summ_s;
	int n_s = cIO.readSumm(cPar.s, separate, summ_s); //readSumm is defined in scr/dtpr.cpp
	// What is n_s?? clearly, an integer, but is it the number of small effect snps? or is it a number of subjects, 
	// like the number of subjects used to get small effect snps?
	vector <POS> inter_s; // what does POS mean here? I get that it's the class for elements of inter_s, but what exactly does it mean?
	// see scr/dtpr.hpp for definition of POS class
	bool badsnp_s[n_s] = {false}; 
	cSP.matchRef(summ_s, ref_bim, inter_s, cPar.mafMax, badsnp_s); //matchRef is defined in scr/dtpr.cpp
	cout << "After filtering, " << inter_s.size() << " small effect SNPs are selected." << endl;
	cout << "We started with " <<summ_s.size() << " small effect SNPs" << endl;
	//read the test bim rs ids into a string vector
	std::vector <std::string> base_nums = readTestBim(cPar.dat_str + ".bim");
	//
	vector<POS> test_inter_s = makePosObjectForTestBim(base_nums, inter_s);
  cout << "test_inter_s[1].snp: " << test_inter_s[1].snp << endl;
  cout << "test_inter_s[1].pos: " << test_inter_s[1].pos << endl;
  cout << "test_inter_s[1].ps: " << test_inter_s[1].ps << endl;  
	vector <INFO> info_s; 
	vector <INFO> test_info_s; 
	
	int num_block_s = cSP.addBlock(inter_s, block_dat, info_s); //addBlock is defined in scr/dtpr.cpp & populates info_s
	int test_num_block_s = cSP.addBlock(test_inter_s, block_dat, test_info_s); //addBlock is defined in scr/dtpr.cpp & populates info_s
	cout << "test_num_block_s: " << test_num_block_s << endl;
	cout << "test_info_s[1].pos: " << test_info_s[1].pos << endl;

	 	// output samll effect badsnps 
	string badsnps_str = cPar.eff + ".badsnps"; 
	ofstream badsnpsFout(badsnps_str.c_str());
	for (size_t i = 0; i < summ_s.size(); ++i) {
		if (badsnp_s[i] == false)
			badsnpsFout << summ_s[i].snp << " " << 0 << endl;
	}
	clearVector(summ_s);

	// large effect
	vector <POS> inter_l;
	vector <INFO> info_l; 
	vector <INFO> test_info_l;
	if (leffstream) {
		// input large effect summary data
		cout << "Reading summary data of large effect SNPs from [" << cPar.l << "]" << endl;
		vector <SUMM> summ_l;
		int n_l = cIO.readSumm(cPar.l, separate, summ_l);
		// vector <POS> inter_l;
		bool badsnp_l[n_l] = {false};
		cSP.matchRef(summ_l, ref_bim, inter_l, cPar.mafMax, badsnp_l);
		vector<POS> test_inter_l = makePosObjectForTestBim(base_nums, inter_l);

		if (inter_l.size() != 0){
			int num_block_l = cSP.addBlock(inter_l, block_dat, info_l); 
		  int test_num_block_l = cSP.addBlock(test_inter_l, block_dat, test_info_l); 
		  
			cout << "After filtering, " << inter_l.size() << " large effect SNPs are selected." << endl;
		} else {
			cout << "After filtering, no large effect SNP is selected." << endl;
		}
		// output large effect badsnps 
		for (size_t i = 0; i < summ_l.size(); ++i) {
			if (badsnp_l[i] == false){
				badsnpsFout << summ_l[i].snp << " " << 1 << endl;
			}
		}
		clearVector(summ_l);
	}
	//read file containing test indices
	arma::uvec test_indicator_pre = read_indices_file(cPar.test_indicator_file) ; //read file containing training set indices
	vector<int> test_indicator =  conv_to<vector<int> >::from(test_indicator_pre);

	// output stream
	string eff_str = cPar.eff + ".txt"; 
	ofstream effFout(eff_str.c_str());
	if (inter_l.size() != 0){
		// fit model
		vector <EFF> eff_s, eff_l; 
		vector<int> idv(n_ref); // idv is an integer vector of length n_ref
		for (int i = 0; i < n_ref; i++) idv[i] = 1; 
		string bed_str = cPar.r + ".bed";
		double t_fitting = cIO.getWalltime();
		double sigma_s = cPar.h / (double)cPar.nsnp;
		cout << "Fitting model..." << endl;
		cDBSF.est(n_ref, 
            cPar.n, //number of subjects with nonmissing trait value
            sigma_s, 
            num_block_s, 
            idv, //vector of length n_ref
            bed_str, 
            info_s, 
            info_l, 
            cPar.t, 
            eff_s, 
            eff_l, 
            test_indicator, 
            cPar.dat_str, //no file ending, like ".bed" or ".bim"
            test_info_s,
            test_info_l); 
		double time_fitting = cIO.getWalltime() - t_fitting;
		cout << "Fitting time: " << time_fitting << " seconds." << endl;

		// output effect 
		for (size_t i = 0; i < eff_l.size(); ++i) {
			double beta_l_noscl = eff_l[i].beta / sqrt(2 * eff_l[i].maf * (1-eff_l[i].maf));//this is like "no scaling" for beta
			if (eff_l[i].snp != "rs" && isinf(beta_l_noscl) == false)
				effFout << eff_l[i].snp << " " << eff_l[i].a1 << " " << eff_l[i].beta << " " << beta_l_noscl << " " << 1 << endl; 
		}
		
		for (size_t i = 0; i < eff_s.size(); ++i) { 
			double beta_s_noscl = eff_s[i].beta / sqrt(2 * eff_s[i].maf * (1-eff_s[i].maf));
			if(eff_s[i].snp.size() != 0 && isinf(beta_s_noscl) == false)
				effFout << eff_s[i].snp << " " << eff_s[i].a1 << " " << eff_s[i].beta << " " << beta_s_noscl << " " << 0 << endl; 
		}
		effFout.close();
	}
	if (inter_l.size() == 0 || !leffstream){
		// fit model
		vector <EFF> eff_s; 
		vector<int> idv(n_ref);
		for (int i = 0; i < n_ref; i++) idv[i] = 1; //set all entries of idv to 1
		string bed_str = cPar.r + ".bed";
		double t_fitting = cIO.getWalltime();
		double sigma_s = cPar.h / (double)cPar.nsnp;
		cout << "Fitting model..." << endl;
		cDBSF.est(n_ref, 
            cPar.n, 
            sigma_s, 
            num_block_s, 
            idv, 
            bed_str, 
            info_s, 
            cPar.t, 
            eff_s, 
            test_indicator, 
            cPar.dat_str, 
            test_info_s); 
		double time_fitting = cIO.getWalltime() - t_fitting;
		cout << "Fitting time: " << time_fitting << " seconds." << endl;

		// output effect 
		for (size_t i = 0; i < eff_s.size(); ++i) { 
			double beta_s_noscl = eff_s[i].beta / sqrt(2 * eff_s[i].maf * (1-eff_s[i].maf));
			if(eff_s[i].snp.size() != 0 && isinf(beta_s_noscl) == false)
				effFout << eff_s[i].snp << " " << eff_s[i].a1 << " " << eff_s[i].beta << " " << beta_s_noscl << " " << 0 << endl; 
		}
	}
	return ;
}
