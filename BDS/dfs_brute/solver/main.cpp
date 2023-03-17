/*
	Algorithm to brute a counter example to BDS algorithm
	
	Author : Gabriel Morete
*/
#include "src/main.h"
#include "src/lemon.h"
#include "src/bds.cpp"
#include "src/nauty_reader.cpp"
// #include "src/stdio_reader.cpp"

using namespace std;

/*
	Conventions
		- Graphs are 0-indexed and simple

	We will build a polyhedron for each graph, one for integral
	and other for fractional formulation. For each matching, we 
	will update the cost function and reoptimize.

	For fractional sulution, we impose primal simplex method.	
*/



/*
	This function receives a LP solution and retuns the restrictions of
	all st-cuts with capacity < 2.
*/
vector<GRBLinExpr> FindMinCuts(double *sol, GRBVar *vars, int n, int m){
	vector<GRBLinExpr> restrictions;
	ListGraph::EdgeMap<double> capacity(G);

	// Build EdgeMap of capacities
	for (ListGraph::EdgeIt e(G); e != INVALID; ++e)
		capacity[e] = sol[G.id(e)]; // Using current solution as capacity

	// Build Gomory-Hu Tree
	GomoryHu<ListGraph, ListGraph::EdgeMap<double> > GMH(G,capacity);
	GMH.run();	                

	// Find all cuts
	for (int i = 0; i < n; i++){
		ListGraph::Node v = G.nodeFromId(i);
		for (int j = 0; j < i; j++){
			ListGraph::Node u = G.nodeFromId(j);

			if (sign(GMH.minCutValue(v, u) - 2) < 0){
				GRBLinExpr cut = 0;
				for (GomoryHu<ListGraph, ListGraph::EdgeMap<double> >::MinCutEdgeIt e(GMH, v, u); e != INVALID; ++e){
					ListGraph::Edge f = e;
					cut += vars[G.id(f)];
				}
				restrictions.push_back(cut);
			}			
		}
	}
	
	return restrictions;
}


/*
	This function builds a fractional cutLP model.
*/
void BuildFractional(GRBModel &frac_model, GRBVar *frac_vars){
	int n = countNodes(G);
	int m = countEdges(G);

	frac_model.set(GRB_IntParam_Method, 0); // Forcing Primal Simplex Method
	// Important, since fractional solution must be an Extreme Point
	
	
	// Creating variables	
	for (ListGraph::EdgeIt e(G); e != INVALID; ++e){
		int u = G.id(G.u(e));
		int v = G.id(G.v(e));
		frac_vars[G.id(e)] = frac_model.addVar(0.0, 1.0, 0.0, GRB_CONTINUOUS, "x_" + to_string(u) + "_" + to_string(v));
	}

	// Add \delta(v) >= 2, constraints
	GRBLinExpr deg2[n];
	for (ListGraph::EdgeIt e(G); e != INVALID; ++e){
		int id = G.id(e);
		int u = G.id(G.u(e));
		int v = G.id(G.v(e));

		deg2[u] += frac_vars[id];
		deg2[v] += frac_vars[id];

	}

	for (int v = 0; v < n; v++)
		frac_model.addConstr(deg2[v] >= 2, "deg2_" + to_string(v));
}

/*
	This function returns a optimum fractional solution to MAP.
	If no solution is found, it returns a all -1 edge map.
*/
void FractionalSolution(ListGraph::EdgeMap<double> &FracSol, GRBModel &frac_model, GRBVar *frac_vars){
	for (ListGraph::EdgeIt e(G); e != INVALID; ++e) 
		FracSol[e] = -1;	

	try {
		int n = countNodes(G);
		int m = countEdges(G);

		GRBLinExpr obj = 0;
		for (ListGraph::EdgeIt e(G); e != INVALID; ++e)
			obj += cost[e] * frac_vars[G.id(e)];

		// Each matching just changes the objective function
		frac_model.setObjective(obj);

		// Optimize model
		frac_model.optimize();

		bool found_feasible = 0;
		while (frac_model.get(GRB_IntAttr_SolCount) > 0 and !found_feasible){
		
			double *sol = frac_model.get(GRB_DoubleAttr_X, frac_vars, m);

			vector<GRBLinExpr> res = FindMinCuts(sol, frac_vars, n, m);

			// If min_cut.fist < 2, need to add constraint
			if (!res.empty()) { // Min cut < 2
				for (GRBLinExpr expr : res)
					frac_model.addConstr(expr >= 2);
				frac_model.optimize();
			}
			else{
				// Found a feasible opt
				// cout << "Obj: " << frac_model.get(GRB_DoubleAttr_ObjVal) << endl;

				for (ListGraph::EdgeIt e(G); e != INVALID; ++e){
					int id = G.id(e);
					int u = G.id(G.u(e));
					int v = G.id(G.v(e));

					FracSol[e] = sol[id]; 

					// cout<<u + 1<<' '<<v + 1<<' '<<abs(sol[u][v])<<endl;
				}
			
				found_feasible = 1;
			}

			delete[] sol;
		}


	} catch(GRBException e) {
		cout << "Error code = " << e.getErrorCode() << endl;
		cout << e.getMessage() << endl;
	} catch(...) {
		cout << "Exception during optimization" << endl;
	}
}



/*
	Wrapper function that call the solvers.
*/
void SolveMapInstance(
	ListGraph::EdgeMap<double> &FracSol,
	ListGraph::EdgeMap<bool> &BDSSol,
	GRBModel &frac_model,
	GRBVar *frac_vars){

	FractionalSolution(FracSol, frac_model, frac_vars);

	if (sign(FracSol[G.edgeFromId(0)]) == -1)
		return;

	BDSAlgorithm(FracSol, BDSSol);

	
}


signed main(int argc, char *argv[]){
	// Start a global gurobi enviroment
	env.set(GRB_IntParam_OutputFlag, 0);
	env.start();

	bool stdio = 0;
	int start = 0;

	for (int i = 1; i < argc; i++){
		string s = argv[i];
		if (s == "-verbose")
			__verbose_mode = 1;
		else if (s == "-stdio")
			stdio = 1;
		else if (s == "-start"){
			s = argv[i + 1];
			start = stoi(s);
			i++;
		} else {
			cout<<"Usage: -stdio -verbose -start n"<<endl;
			return 0;
		}
	}

	RunNautyInput(start);
}

