# -----------------------------------------------------------------
#   Makefile for dbslmm 
# ---------------------------------------------------------------------

# Set the file type
OUTPUTD = dbslmm

# Set the path of library
ARMALIB = /net/mulan/home/yasheng/Cpp/arma/lib

# Put C++ complier 
CXX = g++
#CXX = clang++-11
CPP_FILES = main_dbslmm.cpp dtpr.cpp dbslmm.cpp dbslmmfit.cpp calc_asymptotic_variance.cpp subset_to_test_and_training.cpp helpers.cpp 
HPP_FILES = dtpr.hpp dbslmm.hpp dbslmmfit.hpp calc_asymptotic_variance.hpp subset_to_test_and_training.hpp helpers.hpp



# Set complier flags 
CXXFLAG = -static -fopenmp -O3 -std=c++11 -lm -llapacke -llapack -lblas -Wall
all: $(OUTPUTD)
$(OUTPUTD): $(CPP_FILES) $(HPP_FILES)
	$(CXX) -v $^ -o $(OUTPUTD) $(CXXFLAG) -L $(ARMALIB)

clean:
	rm -f *.o  dbslmm valid
