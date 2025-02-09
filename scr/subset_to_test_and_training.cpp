#include <armadillo>
#include <math.h>       /* floor */
#include <algorithm> /* std::sort, std::set_difference */
#include <fstream> //std::ifstream
#include <string>
#include <boost/lexical_cast.hpp>
#include <iterator>
#include <regex>


#include "subset_to_test_and_training.hpp"


using namespace arma;
using namespace std;

//' Make a vector of zeros and ones with ones at positions indicated
//' 
//' @param ones_positions an arma::uvec specifying the positions where the 1's go
//' @param length the total length of the outputted vector
//' @return an unsigned integer vector of zeros and ones
//' @reference https://gallery.rcpp.org/articles/armadillo-subsetting/

vector<int> make_ones_and_zeroes_vec(arma::uvec ones_positions, unsigned int length){
  cout << "starting make_ones_and_zeroes_vec"<<endl;
  cout <<"result length is: " << length << endl;
  arma::vec result;
  result.zeros(length); //fill vector with all zeros
  //construct a vector for replacing zeroes with ones
  arma::vec ones_vector;
  ones_vector.ones(ones_positions.n_elem);
  cout << "ones_vector has length: " << ones_positions.n_elem << endl;
  // replace zeroes with ones
  //result.elem(ones_positions) = ones_vector;
  result.elem(ones_positions) = ones_vector;
  //convert to vector<int>
  vector<int> out = conv_to<vector<int> >::from(result);
  return out;
}


//' Subset a matrix's rows by indices 
//' 
//' @param mat a matrix, eg., of genotypes, for the entire cohort, with one subject per row
//' @param test_indices vector with subject indices to go into test set
//' @return matrix of genotypes for the subsetted collection of subjects

arma::mat subset(arma::mat matrix, arma::uvec indices){
  arma::mat result = matrix.rows(indices);
  return(result);
}

//' Subset a vector by indices
//' 
//' @param vector a vector, arma::vec
//' @param indices arma::vec of indices to indicate which entries to extract 
//' @return vector of values for the subsetted collection of subjects

arma::vec subset(arma::vec vector, arma::uvec indices){
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
//' @return integer vector containing the complement of test_indices to indicate membership in training data set

arma::uvec get_complementary_indices(arma::uvec indices, int sample_size){
  //convert to std::vector
  std::vector<int> test_std = arma::conv_to<std::vector<int> >::from(indices);
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
  arma::uvec result = arma::conv_to<arma::uvec >::from(v);
  return result;
}

//' Convert std::vector <string> to indices
//' 
//' @param string a string vector
//' @return arma::uvec vector, for use as indices in subsetting armadillo matrices or vectors

arma::uvec convert_string_to_indices(std::vector <std::string> in_string){
  vector<int> vectorOfIntegers;
  castContainer(in_string, vectorOfIntegers);
  arma::uvec result = conv_to< arma::uvec >::from(vectorOfIntegers);
  return (result);
} 



//' Split a string vector with a regex
//' 
//' @references https://stackoverflow.com/questions/14265581/parse-split-a-string-in-c-using-string-delimiter-standard-c

std::vector<std::string> split(const std::string str, const std::string regex_str)
{
  std::regex regexz(regex_str);
  std::vector<std::string> list(std::sregex_token_iterator(str.begin(), str.end(), regexz, -1),
                                std::sregex_token_iterator());
  return list;
}

//' Read first column of a text file containing one integer per line with space delimiters
//' 
//' @param filepath file path for the indices file
//' @references https://www.youtube.com/watch?v=EaHFhms_Shw, https://stackoverflow.com/questions/2706076/new-how-would-i-read-a-file-that-has-3-columns-and-each-column-contains-100-n

arma::uvec read_indices_file(const string filepath){
  ifstream infile;
  infile.open(filepath.c_str(), ios::in); //read mode
/*  if(infile.fail()) // checks to see if file opened 
  { 
    cout << "error - file didn't open" << endl; 
    return 1; // no point continuing if the file didn't open...
  }*/ 
  string line;
  std::vector<std::string> result;
  while(getline(infile, line)) 
  { 
    std::vector<std::string> l0 = split(line, " "); //split line with space delimiter 
    result.push_back(l0[0]); // append only the first entry in line
  } 
  infile.close(); 
  arma::uvec out = convert_string_to_indices(result);  
  return (out); 
}


//' a single string to single integer function
//' @references https://www.py4u.net/discuss/90965

template<typename C1, typename C2>
void castContainer(const C1& source, C2& destination)
{
  typedef typename C1::value_type source_type;
  typedef typename C2::value_type destination_type;
  destination.resize(source.size());
  std::transform(source.begin(), source.end(), destination.begin(), boost::lexical_cast<destination_type, source_type>);
}

template<typename T, typename T2>
std::vector<T>& operator<<(std::vector<T>& v, T2 t)
{
  v.push_back(T(t));
  return v;
}
