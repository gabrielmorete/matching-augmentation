CPP = g++
OPT	= -Ofast	
STD = -std=c++17#-std=c++11 -D_GLIBCXX_USE_CXX11_ABI=0
# SAN = -fsanitize=address -fsanitize=undefined -fno-sanitize-recover
# STL = -D_GLIBCXX_DEBUG -D_GLIBCXX_DEBUG_PEDANTIC -D_FORTIFY_SOURCE=2

% : %.cpp
	$(CPP) $(STD) $(OPT) $(SAN) $(STL) -o $@ $<
