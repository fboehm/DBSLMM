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
CPP_FILES = main_dbslmm.cpp dtpr.cpp dbslmm.cpp dbslmmfit.cpp calc_asymptotic_variance.cpp subset_to_test_and_training.cpp 
HPP_FILES = dtpr.hpp dbslmm.hpp dbslmmfit.hpp calc_asymptotic_variance.hpp subset_to_test_and_training.hpp 



# Set complier flags 
CXXFLAG = -static -fopenmp -O3 -std=c++11 -lm -llapacke -llapack -lblas -Wall
all: $(OUTPUTD)
$(OUTPUTD): main_dbslmm.o dtpr.o dbslmm.o dbslmmfit.o calc_asymptotic_variance.o subset_to_test_and_training.o
	$(CXX) -v main_dbslmm.o dtpr.o dbslmm.o dbslmmfit.o -o $(OUTPUTD) $(CXXFLAG) -L $(ARMALIB)
main_dbslmm.o: main_dbslmm.cpp
	$(CXX) -c main_dbslmm.cpp
dbslmm.o: dbslmm.cpp dbslmm.hpp subset_to_test_and_training.hpp
	$(CXX) -c $^  $(CXXFLAG)
dbslmmfit.o: dbslmmfit.cpp dbslmmfit.hpp calc_asymptotic_variance.hpp subset_to_test_and_training.hpp
	$(CXX) -c $^ $(CXXFLAG)
dtpr.o: dtpr.cpp dtpr.hpp 
	$(CXX) -c dtpr.cpp $(CXXFLAG)
calc_asymptotic_variance.o: calc_asymptotic_variance.cpp calc_asymptotic_variance.hpp 
	$(CXX) -c calc_asymptotic_variance.cpp $(CXXFLAG)
subset_to_test_and_training.o: subset_to_test_and_training.cpp subset_to_test_and_training.hpp 
	$(CXX) -c subset_to_test_and_training.cpp $(CXXFLAG)

clean:
	rm -f *.o  dbslmm valid
