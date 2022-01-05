#include <armadillo>
#include <math.h>       /* floor */
#include <algorithm> /* std::sort, std::set_difference */
#include <ifstream>
#include <string>
#include "subset_to_test_and_training.hpp"


using namespace arma;
using namespace std;

//' Fisher-Yates shuffle algorithm for a permutation of integers
//' 
//' @param size number of samples to draw
//' @param max_size number of possible objects
//' @param gen for pseudorandomness
//' @references https://ideone.com/KvvCda

std::vector<int> fisher_yates_shuffle(std::size_t size, 
                                      std::size_t max_size, 
                                      std::mt19937& gen)
{
  //assert(size < max_size);
  std::vector<int> res(size);
  
  for(std::size_t i = 0; i != max_size; ++i) {
    std::uniform_int_distribution<> dis(0, i);
    std::size_t j = dis(gen);
    if (j < res.size()) {
      if (i < res.size()) {
        res[i] = res[j];
      }
      res[j] = i;
    }
  }
  return res;
}



//' (Pseudo-)Randomly sample indices, eg., to determine test set membership
//' 
//' @param n_obs total number of subjects (test plus training)
//' @param test_proportion proportion of subjects to place in test set
//' @param seed a positive integer seed for the pseudo-RNG
//' @return  integers to indicate test set membership

arma::Col<arma::uword> get_test_indices(int n_obs, 
                                        double test_proportion, 
                                        unsigned int seed){
  // calculate number of subjects to put into test set
  int n_test = floor(test_proportion * n_obs);
  // pseudo-random stuff (Mersenne twister)
  std::random_device rd;
  std::mt19937 gen(rd());
  gen.seed(seed);
  //randomly sample without replacement from the integers 0,1,...,n_obs - 1 and return n_test of them.
  std::vector<int> sampled = fisher_yates_shuffle(n_test, n_obs, gen);
  arma::Col<arma::uword> result = arma::conv_to< arma::Col<arma::uword> >::from(sampled);
  return(result);
}

//' Subset a matrix's rows by indices 
//' 
//' @param mat a matrix, eg., of genotypes, for the entire cohort, with one subject per row
//' @param test_indices vector with subject indices to go into test set
//' @return matrix of genotypes for the subsetted collection of subjects

arma::mat subset(arma::mat matrix, arma::Col<arma::uword> indices){
  arma::mat result = matrix.rows(indices);
  return(result);
}

//' Subset a vector by indices
//' 
//' @param vector a vector, arma::vec
//' @param indices arma::vec of indices to indicate which entries to extract 
//' @return vector of values for the subsetted collection of subjects

arma::vec subset(arma::vec vector, arma::Col<arma::uword> indices){
  arma::vec result = vector.elem(indices);
  return(result);
}

//' Construct an integer vector from start to end, for integers start and end
//' 
//' @param start smallest and first integer value
//' @param end largest and last integer value
//' @return integer vector, start, start + 1, ..., end

std::vector<int> make_integer_vector(int start, int end){
  std::vector<int> myVec;
  for( int i = start; i <= end; i++ ) //spacing between value is 1.
    myVec.push_back( i );
  return myVec;
}

//' Get complementary indices for, eg, test data or training data
//' 
//' @param indices indices for subjects to be placed into test data set
//' @param sample_size total combined sample size, training and test together
//' @return arma::Col integer vector containing the complement of test_indices to indicate membership in training data set

arma::Col<arma::uword> get_complementary_indices(arma::Col<arma::uword> indices, int sample_size){
  //convert to std::vector
  std::vector<int> test_std = arma::conv_to<std::vector<int> >::from(test_indices);
  std::sort(test_std.begin(), test_std.end()); //essentially overwrites test_std with the sorted vector, smallest to largest
  std::vector<int> all_indices = make_integer_vector(0, sample_size - 1);
  std::sort(all_indices.begin(), all_indices.end());
  std::vector<int> v(sample_size);
  std::vector<int>::iterator it;
  it = std::set_difference(all_indices.begin(), 
                      all_indices.end(), 
                      test_std.begin(), 
                      test_std.end(), 
                      v.begin() );
  v.resize(it - v.begin()); //result is in v
  arma::Col<arma::uword> result = arma::conv_to<arma::Col<arma::uword> >::from(v);
  return result;
}

//' Convert std::vector <string> to arma::Col<arma::uword>
//' 
//' @param string a string vector
//' @return arma::Col<arma::uword> vector, for use as indices in subsetting armadillo matrices or vectors

arma::Col<arma::uword> convert_string_to_Col(std::vector<std::string> string){
  arma::Col<arma::uword> result = conv_to< arma::Col<arma::uword> >::from(string);
  return (result);
} 

//' Read first - and only - column of a text file containing one integer per line
//' 
//' @param filepath file path for the indices file
//' @references https://www.youtube.com/watch?v=EaHFhms_Shw, https://stackoverflow.com/questions/2706076/new-how-would-i-read-a-file-that-has-3-columns-and-each-column-contains-100-n

arma::Col <arma::uword> read_indices_file(const char filepath){
  ifstream infile;
  infile.open(filepath, ios::in); //read mode
  if(infile.fail()) // checks to see if file opened 
  { 
    cout << "error - file didn't open" << endl; 
    return 1; // no point continuing if the file didn't open...
  } 
  string line;
  std::vector<string> result;
  while(string::getline(infile, line)) 
  { 
    result.push_back(line);
  } 
  infile.close(); 
  arma::Col<arma::uword> out = convert_string_to_Col(result);  
  return (out); 
}
