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

using namespace std;

#define dbg(x)  cout << #x << " = " << x << endl


// Safe handling doubles
const double EPS = 1e-3;
int sign(double x) { return (x > EPS) - (x < -EPS); }




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
vector<GRBLinExpr> FindMinCuts(double *sol, GRBVar *vars, int n, int m);

/*
	This function builds a fractional cutLP model.
*/
void BuildFractional(GRBModel &frac_model, GRBVar *frac_vars);

/*
	This function returns a optimum fractional solution to MAP.
	If no solution is found, it returns a all -1 edge map.
*/
void FractionalSolution(ListGraph::EdgeMap<double> &FracSol, GRBModel &frac_model, GRBVar *frac_vars);

/*
	This function builds a integral cutLP model.
*/
void BuildIntegral(GRBModel &int_model, GRBVar *int_vars);

/*
	This function returns a optimum integer solution to MAP.
	If no solution is found, it returns a all -1 edge map.
*/
void IntegerSolution(ListGraph::EdgeMap<int> &IntSol, GRBModel &int_model, GRBVar *int_vars);

/*
	Wrapper function that call the solvers.
*/
void SolveMapInstance(
	ListGraph::EdgeMap<double> &FracSol,
	ListGraph::EdgeMap<int> &IntSol,
	ListGraph::EdgeMap<bool> &BDSSol,
	GRBModel &frac_model,
	GRBVar *frac_vars,
	GRBModel &int_model,
	GRBVar *int_vars);


/*
	Callback function class.
*/
class MinimumCut: public GRBCallback {
	public:
		GRBVar* vars;
		int n, m;
		MinimumCut(GRBVar* xvars, int xn, int xm);
	protected:
		void callback();
};


#endif
// /*
// /*
// 	This function returns a optimum integer solution to MAP.
// 	If no solution is found, it returns a all -1 edge map.
// */
// void IntegerSolution(ListGraph::EdgeMap<int> &IntSol){

// 	for (ListGraph::EdgeIt e(G); e != INVALID; ++e) 
// 		IntSol[e] = -1;	

// 	try {
// 		int n = countNodes(G);
// 		int m = countEdges(G);

// 		// Create an environment
// 		GRBEnv env = GRBEnv(true);
// 		// env.set("LogFile", "MAPInteger.log"); // Output may be large
// 		env.set(GRB_IntParam_OutputFlag, 0);
// 		env.start();

// 		// Create an empty model
// 		GRBModel model = GRBModel(env);
// 		model.set(GRB_IntParam_LazyConstraints, 1); // Allow callback constraints
		
// 		GRBVar vars[m];
// 		int node_u[m], node_v[m];

// 		// Setting the correct UB and OBJ.	
// 		for (ListGraph::EdgeIt e(G); e != INVALID; ++e){
// 			int u = G.id(G.u(e));
// 			int v = G.id(G.v(e));
// 			vars[G.id(e)] = model.addVar(0.0, 1.0, cost[e], GRB_BINARY, "x_" + to_string(u) + "_" + to_string(v));
// 			node_u[G.id(e)] = u;
// 			node_v[G.id(e)] = v;
// 		}

// 		// Add \delta(v) >= 2, constraints
// 		GRBLinExpr deg2[n];
// 		for (ListGraph::EdgeIt e(G); e != INVALID; ++e){
// 			int id = G.id(e);
// 			int u = G.id(G.u(e));
// 			int v = G.id(G.v(e));

// 			deg2[u] += vars[id];
// 			deg2[v] += vars[id];

// 		}

// 		for (int v = 0; v < n; v++)
// 			model.addConstr(deg2[v] >= 2, "deg2_" + to_string(v));

// 		// Set callback function
//     	MinimumCut cb = MinimumCut(vars, n, m);
//     	model.setCallback(&cb);
		
// 		// Optimize model
// 		model.optimize();

// 		// Found optimal solution
// 		if (model.get(GRB_IntAttr_SolCount) > 0){
// 			// cout << "Obj: " << model.get(GRB_DoubleAttr_ObjVal) << endl;

// 			double *sol = model.get(GRB_DoubleAttr_X, vars, m);

// 			for (ListGraph::EdgeIt e(G); e != INVALID; ++e){
// 				int id = G.id(e);
// 				int u = G.id(G.u(e));
// 				int v = G.id(G.v(e));

// 				if (sol[id] > 0.5)
// 					IntSol[e] = 1;
// 				else
// 					IntSol[e] = 0;
// 				// cout<<u + 1<<' '<<v + 1<<' '<<abs(sol[u][v])<<endl;

// 				// Sanity check
// 				assert((node_u[id] == u) and (node_v[id]) == v);
// 			}

// 			delete[] sol;
// 		}

// 	} catch(GRBException e) {
// 		cout << "Error code = " << e.getErrorCode() << endl;
// 		cout << e.getMessage() << endl;
// 	} catch(...) {
// 		cout << "Exception during optimization" << endl;
// 	}
// }

// /*
// 	This function returns a optimum fractional solution to MAP.
// 	If no solution is found, it returns a all -1 edge map.
// */
// void FractionalSolution(ListGraph::EdgeMap<double> &FracSol){
// 	for (ListGraph::EdgeIt e(G); e != INVALID; ++e) 
// 		FracSol[e] = -1;	

// 	try {
// 		int n = countNodes(G);
// 		int m = countEdges(G);

// 		// Create an environment
// 		GRBEnv env = GRBEnv(true);
// 		// env.set("LogFile", "MAPFractional.log");  // Output may be large
// 		env.set(GRB_IntParam_OutputFlag, 0);
// 		env.start();

// 		// Create an empty model
// 		GRBModel model = GRBModel(env);
// 		model.set(GRB_IntParam_Method, 0); // Forcing Primal Simplex Method
// 		// Important, since fractional solution must be an Extreme Point
		
		
// 		GRBVar vars[m];
// 		int node_u[m], node_v[m];

// 		// Setting the correct UB and OBJ.	
// 		for (ListGraph::EdgeIt e(G); e != INVALID; ++e){
// 			int u = G.id(G.u(e));
// 			int v = G.id(G.v(e));
// 			vars[G.id(e)] = model.addVar(0.0, 1.0, cost[e], GRB_CONTINUOUS, "x_" + to_string(u) + "_" + to_string(v));
// 			node_u[G.id(e)] = u;
// 			node_v[G.id(e)] = v;
// 		}

// 		// Add \delta(v) >= 2, constraints
// 		GRBLinExpr deg2[n];
// 		for (ListGraph::EdgeIt e(G); e != INVALID; ++e){
// 			int id = G.id(e);
// 			int u = G.id(G.u(e));
// 			int v = G.id(G.v(e));

// 			deg2[u] += vars[id];
// 			deg2[v] += vars[id];

// 		}

// 		for (int v = 0; v < n; v++)
// 			model.addConstr(deg2[v] >= 2, "deg2_" + to_string(v));

// 		// Optimize model
// 		model.optimize();

// 		bool found_feasible = 0;
// 		while (model.get(GRB_IntAttr_SolCount) > 0 and !found_feasible){
		
// 			double *sol = model.get(GRB_DoubleAttr_X, vars, m);
// 			pair<double, vector<Edge> > min_cut = FindMinCut(sol, n, m);

// 			// If min_cut.fist < 2, need to add constraint
// 			if (sign(min_cut.first - 2.0) < 0) { // Min cut < 2
// 				GRBLinExpr expr = 0;
				
// 				// All all edges of the global min cut as a contraint.
// 				for (Edge e : min_cut.second){
// 					int id = G.id(e);
// 					int u = G.id(G.u(e));
// 					int v = G.id(G.v(e));				
					
// 					expr += vars[id];
// 					assert((node_u[id] == u) and (node_v[id]) == v); // Sanity check
// 				}

// 				model.addConstr(expr >= 2);
// 				model.optimize();
// 			}
// 			else{
// 				// Found a feasible opt
// 				// cout << "Obj: " << model.get(GRB_DoubleAttr_ObjVal) << endl;

// 				for (ListGraph::EdgeIt e(G); e != INVALID; ++e){
// 					int id = G.id(e);
// 					int u = G.id(G.u(e));
// 					int v = G.id(G.v(e));

// 					FracSol[e] = sol[id]; 

// 					// cout<<u + 1<<' '<<v + 1<<' '<<abs(sol[u][v])<<endl;

// 					assert((node_u[id] == u) and (node_v[id]) == v); // Sanity check
// 				}
			
// 				found_feasible = 1;
// 			}

// 			delete[] sol;
// 		}


// 	} catch(GRBException e) {
// 		cout << "Error code = " << e.getErrorCode() << endl;
// 		cout << e.getMessage() << endl;
// 	} catch(...) {
// 		cout << "Exception during optimization" << endl;
// 	}
// }

// pair<double, vector<Edge> > FindMinCut(double *sol, int n, int m){
// 	ListGraph::EdgeMap<double> capacity(G);

// 	// Build EdgeMap of capacities
// 	for (ListGraph::EdgeIt e(G); e != INVALID; ++e)
// 		capacity[e] = sol[G.id(e)]; // Using current solution as capacity

// 	// Build Gomory-Hu Tree
// 	GomoryHu<ListGraph, ListGraph::EdgeMap<double> > GMH(G,capacity);
// 	GMH.run();	                

// 	ListGraph::Node min_v = G.nodeFromId(0);
// 	ListGraph::Node min_u = G.nodeFromId(1);
// 	for (ListGraph::NodeIt v(G); v != INVALID; ++v)
// 		for (ListGraph::NodeIt u(G); u != INVALID; ++u){
// 			if (GMH.minCutValue(v, u) < GMH.minCutValue(min_v, min_u)){
// 				min_v = v;
// 				min_u = u;
// 			}
// 	}

// 	vector<Edge> min_cut;
// 	for (GomoryHu<ListGraph, ListGraph::EdgeMap<double> >::MinCutEdgeIt e(GMH, min_v, min_u); e != INVALID; ++e)
// 		min_cut.push_back(e);
	
// 	return make_pair(GMH.minCutValue(min_v, min_u), min_cut);
// }

