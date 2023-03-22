#ifndef MAIN_DEF
#define MAIN_DEF

#include <iostream>
#include <fstream>
#include <experimental/filesystem>
#include <vector>
#include <array>
#include <cassert>
#include <string>
#include <set>
#include <algorithm>
#include <cmath>
#include "gurobi_c++.h"
#include "lemon.h"
#include <omp.h>

using namespace std;

#define dbg(x)  cout << #x << " = " << x << endl

// Safe handling doubles
const double EPS = 1e-4;
/*
	-1 if x < -EPS
	0 if  -EPS <= x <= EPS
	1 if EPS < x
*/ 
int sign(double x) { return (x > EPS) - (x < -EPS); } 


// Print extra information
bool __verbose_mode = 0;


/*
	I will define a global enviroment. 
	Each graph will have a single model. Different matchings
	will be just different cost functions.
*/
GRBEnv env = GRBEnv(true);

/*
	This function receives a LP solution and returns the restrictions of
	all st-cuts with capacity < 2.
*/
vector<GRBLinExpr> FindMinCuts(double *sol, GRBVar *vars, int n, int m, ListGraph &G);

/*
	This function builds a fractional cutLP model.
*/
void BuildFractional(GRBModel &frac_model, GRBVar *frac_vars, ListGraph &G);

/*
	This function returns a optimum fractional solution to MAP.
	If no solution is found, it returns a all -1 edge map.
*/
void FractionalSolution(ListGraph::EdgeMap<double> &FracSol, 
	GRBModel &frac_model, 
	GRBVar *frac_vars, 
	ListGraph &G);

/*
	This function builds a integral cutLP model.
*/
void BuildIntegral(GRBModel &int_model, GRBVar *int_vars, ListGraph &G);

/*
	This function returns a optimum integer solution to MAP.
	If no solution is found, it returns a all -1 edge map.
*/
void IntegerSolution(ListGraph::EdgeMap<int> &IntSol, 
	GRBModel &int_model, 
	GRBVar *int_vars, 
	ListGraph &G);

/*
	Wrapper function that call the solvers.
*/
void SolveMapInstance(
	ListGraph::EdgeMap<int> &cost,
	ListGraph::EdgeMap<double> &FracSol,
	ListGraph::EdgeMap<int> &IntSol,
	ListGraph::EdgeMap<bool> &BDSSol,
	GRBModel &frac_model,
	GRBVar *frac_vars,
	GRBModel &int_model,
	GRBVar *int_vars,
	ListGraph &G,
	BDSAlgorithm &G);


/*
	Callback function class.
*/
class MinimumCut: public GRBCallback {
	public:
		ListGraph *G;
		GRBVar* vars;
		int n, m;
		MinimumCut(GRBVar* _vars, int _n, int _m, ListGraph &_G);
	protected:
		void callback();
};

#endif