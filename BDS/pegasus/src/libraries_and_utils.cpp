#ifndef LIBRARIES_UTILS
#define LIBRARIES_UTILS

#include <iostream>
#include <fstream>
#include <experimental/filesystem>
#include <vector>
#include <array>
#include <cassert>
#include <string>
#include <set>
#include <algorithm>

using namespace std;

// Safe handling doubles
const double EPS = 1e-3;
int sign(double x) { return (x > EPS) - (x < -EPS); }

#define dbg(x)  cout << #x << " = " << x << endl


#endif