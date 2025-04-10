# Makefile
# Gurobi is installed in.
# /Library/gurobi1001/macos_universal2/include/


PLATFORM = macos_universal2
PATH	 = /Library/gurobi1001/macos_universal2
INC      = $(PATH)/include/
CPP      = g++
STD      = -std=c++17
OPT		 = -O3 -static
CARGS    = -m64 -g
CPPLIB   = -L $(PATH)/lib -lgurobi_c++ -lgurobi100
LDFLAGS  = -L/usr/local/opt/llvm/lib
CPPFLAGS = -I/usr/local/opt/llvm/include

% : %.cpp
	$(CPP) $(CARGS) $(OPT) -o $@ $< -I$(INC) $(CPPLIB) -fopenmp -I/usr/local/include -L/usr/local/lib -lomp -lm -lemon

# clean:
# 	rm -rf *.o *_c *_c++ *.class *.log *.rlp *.lp *.bas *.ilp *.mps *.prm *.dSYM; \
# 	if [ -d $(GRBAPP) ]; then \
# 		cd $(GRBAPP); \
# 		find . ! \( -name "gurobi*.csproj" -o -name . \) -exec rm -rf {} +; \
# 	fi

