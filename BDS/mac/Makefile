# Makefile
# Gurobi is installed in.
# /Library/gurobi1001/macos_universal2/include/


PLATFORM = macos_universal2
PATH	 = /Library/gurobi1001/macos_universal2
INC      = $(PATH)/include/
CPP      = g++
STD      = -std=c++17
CARGS    = -m64 -g
CPPLIB   = -L $(PATH)/lib -lgurobi_c++ -lgurobi100


% : %.cpp
	$(CPP) $(CARGS) -o $@ $< -I$(INC) $(CPPLIB) -lm -lemon

# clean:
# 	rm -rf *.o *_c *_c++ *.class *.log *.rlp *.lp *.bas *.ilp *.mps *.prm *.dSYM; \
# 	if [ -d $(GRBAPP) ]; then \
# 		cd $(GRBAPP); \
# 		find . ! \( -name "gurobi*.csproj" -o -name . \) -exec rm -rf {} +; \
# 	fi

