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


# Set complier flags 
CXXFLAG = -static -fopenmp -O3 -std=c++11 -lm -llapacke -llapack -lblas -Wall
all: $(OUTPUTD)
$(OUTPUTD): main_dbslmm.cpp dtpr.cpp dbslmm.cpp dbslmmfit.cpp 
	$(CXX) main_dbslmm.cpp dtpr.cpp dbslmm.cpp dbslmmfit.cpp -o $(OUTPUTD) $(CXXFLAG) -L $(ARMALIB)
main_dbslmm.o: main_dbslmm.cpp
	$(CXX) -c main_dbslmm.cpp
dbslmm.o: dbslmm.cpp dbslmm.hpp
	$(CXX) -c dbslmm.cpp  $(CXXFLAG)
dbslmmfit.o: dbslmmfit.cpp dbslmmfit.hpp
	$(CXX) -c dbslmmfit.cpp $(CXXFLAG)
dtpr.o: dtpr.cpp dtpr.hpp 
	$(CXX) -c dtpr.cpp $(CXXFLAG)


clean:
	rm -f *.o  dbslmm valid
