# Makefile
# Gurobi is installed in.
# /opt/uw/gurobi/8.1.1/

PLATFORM = linux64
GUROBI_PATH = /opt/uw/gurobi/8.1.1
INC = $(GUROBI_PATH)/include/
CPP = g++
STD = -std=c++11 -D_GLIBCXX_USE_CXX11_ABI=0
CARGS = -m64 -O3
GUROBI_CPPLIB = -L $(GUROBI_PATH)/lib -lgurobi_c++ -lgurobi81

# The following flags are only for debuging purposes
 SAN = -fsanitize=address -fsanitize=undefined -fno-sanitize-recover
 STL = -D_GLIBCXX_DEBUG -D_GLIBCXX_DEBUG_PEDANTIC -D_FORTIFY_SOURCE=2

% : %.cpp
	$(CPP) $(CARGS) $(STD) -o $@ $< -I$(INC) $(GUROBI_CPPLIB) -lm -lemon -fopenmp -lstdc++fs $(SAN) $(STL)
