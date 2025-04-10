# Makefile
# Gurobi is installed in.
# /opt/gurobi1002/

PLATFORM = linux64
GUROBI_PATH = /opt/gurobi1002/linux64
LEMON_PATH = /u/gazevedo/ResearchDocuments/waterloo-bepe/BDS/cpu149/lemon
INC = $(GUROBI_PATH)/include/
CPP = g++
STD = -std=c++11 -D_GLIBCXX_USE_CXX11_ABI=0
CARGS = -m64 -O3
GUROBI_CPPLIB = -L $(GUROBI_PATH)/lib -lgurobi_c++ -lgurobi100

# The following flags are only for debuging purposes
# SAN = -fsanitize=address -fsanitize=undefined -fno-sanitize-recover
# STL = -D_GLIBCXX_DEBUG -D_GLIBCXX_DEBUG_PEDANTIC -D_FORTIFY_SOURCE=2

% : %.cpp
	$(CPP) $(CARGS) $(STD) -o $@ $< -I$(INC) $(GUROBI_CPPLIB) -lm -lemon -fopenmp -lstdc++fs $(SAN) $(STL)
