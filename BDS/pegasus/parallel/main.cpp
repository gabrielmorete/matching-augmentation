/*
	Matching Augmentation problem
		- Linear Relaxation
		- Integer Formulation
		- BDS Approximation Algorithm 

	
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
vector<GRBLinExpr> FindMinCuts(double *sol, GRBVar *vars, int n, int m, ListGraph &G){
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
	Constructor of min cut class
*/
MinimumCut::MinimumCut(GRBVar* xvars, int xn, int xm, ListGraph &_G){
	vars = xvars;
	n = xn;
	m = xm;
	G = &_G;
}

/*
	Separator function for the MIP. 
	If the solution is not a 2ECSS it adds a cut
	separating one 2ECC.
*/
void MinimumCut::callback(){
	try {
		if (where == GRB_CB_MIPSOL){
			/*
				Solver found an integral optimal solution for the
				current formulation, must check if if there is a
				bridge or a cut.
			*/

			double *x = getSolution(vars, m);

			ListGraph::EdgeMap<bool> in_sol(*G);
			
			for (ListGraph::EdgeIt e(*G); e != INVALID; ++e){
				int id = (*G).id(e);
				int u = (*G).id((*G).u(e));
				int v = (*G).id((*G).v(e));

				if (x[id] > 0.5)
					in_sol[e] = 1;
			}

			ListGraph::NodeMap<bool> ones((*G), 1); // Must be spanning
			// H is a spanning subgraph with all edges in the solution
			SubGraph<ListGraph> H((*G), ones, in_sol);

			ListGraph::NodeMap<int> ebcc((*G), -1);
			biEdgeConnectedComponents(H, ebcc);

			int max_cmp = 0;
			for (ListGraph::NodeIt v((*G)); v != INVALID; ++v)
				max_cmp = max(max_cmp, ebcc[v]);

			if (max_cmp > 0){ // Not 2ECSS, must add a cut 

				// The cut will be all edges crossing the cut of the ebcc with id 0
				GRBLinExpr expr = 0;

				for (ListGraph::EdgeIt e((*G)); e != INVALID; ++e){
					ListGraph::Node u = (*G).u(e);
					ListGraph::Node v = (*G).v(e);
					
					if (ebcc[u] == 0 and ebcc[v] != 0)
						expr += vars[(*G).id(e)];
			
					if (ebcc[v] == 0 and ebcc[u] != 0)
						expr += vars[(*G).id(e)];
				}

				addLazy(expr >= 2);
			}

			delete[] x;
		}	
	} 
	catch (GRBException e){
		cout << "Error number: " << e.getErrorCode() << endl;
		cout << e.getMessage() << endl;
	} 
	catch (...){
		cout << "Error during callback" << endl;
	}
}

/*
	This function builds a fractional cutLP model.
*/
void BuildFractional(GRBModel &frac_model, GRBVar *frac_vars, ListGraph &G){
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
void FractionalSolution(ListGraph::EdgeMap<int> &cost,
	ListGraph::EdgeMap<double> &FracSol, 
	GRBModel &frac_model, 
	GRBVar *frac_vars,
	ListGraph &G){
	
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

			vector<GRBLinExpr> res = FindMinCuts(sol, frac_vars, n, m, G);

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
	This function builds a integral cutLP model.
*/
void BuildIntegral(GRBModel &int_model, GRBVar *int_vars, ListGraph &G){
	int n = countNodes(G);
	int m = countEdges(G);

	int_model.set(GRB_IntParam_LazyConstraints, 1); // Allow callback constraints

	for (ListGraph::EdgeIt e(G); e != INVALID; ++e){
		int u = G.id(G.u(e));
		int v = G.id(G.v(e));
		int_vars[G.id(e)] = int_model.addVar(0.0, 1.0, 0.0, GRB_BINARY, "x_" + to_string(u) + "_" + to_string(v));
	}

	// Add \delta(v) >= 2, constraints
	GRBLinExpr deg2[n];
	for (ListGraph::EdgeIt e(G); e != INVALID; ++e){
		int id = G.id(e);
		int u = G.id(G.u(e));
		int v = G.id(G.v(e));

		deg2[u] += int_vars[id];
		deg2[v] += int_vars[id];

	}

	for (int v = 0; v < n; v++)
		int_model.addConstr(deg2[v] >= 2, "deg2_" + to_string(v));
}


/*
	This function returns a optimum integer solution to MAP.
	If no solution is found, it returns a all -1 edge map.
*/
void IntegerSolution(ListGraph::EdgeMap<int> &cost,
	ListGraph::EdgeMap<int> &IntSol, 
	ListGraph::EdgeMap<bool> &BDSSol, 
	GRBModel &int_model, 
	GRBVar *int_vars,
	ListGraph &G){

	for (ListGraph::EdgeIt e(G); e != INVALID; ++e) 
		IntSol[e] = -1;	

	try {
		int n = countNodes(G);
		int m = countEdges(G);

		GRBLinExpr obj = 0;
		for (ListGraph::EdgeIt e(G); e != INVALID; ++e)
			obj += cost[e] * int_vars[G.id(e)];

		int_model.setObjective(obj);

		// Setting bds sol as a starting feasible solution
		for (ListGraph::EdgeIt e(G); e != INVALID; ++e)
			int_vars[G.id(e)].set(GRB_DoubleAttr_Start, BDSSol[e]);

		// Optimize model
		int_model.optimize();

		// Found optimal solution
		if (int_model.get(GRB_IntAttr_SolCount) > 0){
			// cout << "Obj: " << int_model.get(GRB_DoubleAttr_ObjVal) << endl;

			double *sol = int_model.get(GRB_DoubleAttr_X, int_vars, m);

			for (ListGraph::EdgeIt e(G); e != INVALID; ++e){
				int id = G.id(e);
				int u = G.id(G.u(e));
				int v = G.id(G.v(e));

				if (sol[id] > 0.5)
					IntSol[e] = 1;
				else
					IntSol[e] = 0;
				// cout<<u + 1<<' '<<v + 1<<' '<<abs(sol[u][v])<<endl;
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
	ListGraph::EdgeMap<int> &cost,
	ListGraph::EdgeMap<double> &FracSol,
	ListGraph::EdgeMap<int> &IntSol,
	ListGraph::EdgeMap<bool> &BDSSol,
	GRBModel &frac_model,
	GRBVar *frac_vars,
	GRBModel &int_model,
	GRBVar *int_vars,
	ListGraph &G){

	FractionalSolution(cost, FracSol, frac_model, frac_vars, G);

	if (sign(FracSol[G.edgeFromId(0)]) == -1)
		return;

	BDSAlgorithm(cost, FracSol, BDSSol, G);

	// If fractional solution is integral, no need to solve a MIP
	bool is_integral = 1;
	for (ListGraph::EdgeIt e(G); e != INVALID; ++e)
		if ((sign(FracSol[e]) != 0) and (sign(FracSol[e] - 1.0) != 0)) // not 0 nor 1
			is_integral = 0;

	if (is_integral){
		for (ListGraph::EdgeIt e(G); e != INVALID; ++e)
			IntSol[e] = FracSol[e];
		return;
	}

	double frac_cost = 0, BDS_cost = 0;

	for (ListGraph::EdgeIt e(G); e != INVALID; ++e){
		frac_cost += cost[e] * FracSol[e];
		BDS_cost += cost[e] * ((int) BDSSol[e]);
	}

	// ceil(frac_cost) is a lower bound on the integral cost
	if (sign(ceil(frac_cost) - BDS_cost) == 0){ // BDS sol is opt
		for (ListGraph::EdgeIt e(G); e != INVALID; ++e)
			IntSol[e] = BDSSol[e];

		return;
	}

	IntegerSolution(cost, IntSol, BDSSol, int_model, int_vars, G);
}


signed main(int argc, char *argv[]){
	// Start a global gurobi enviroment
	env.set(GRB_IntParam_OutputFlag, 0);
	env.start();

	bool stdio = 0;
	int start = 0;
	int n_threads = 1;

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
		} else if (s == "-threads"){
			s = argv[i + 1];
			n_threads = stoi(s);
			i++;
		} else {
			cout<<"Usage: -stdio -verbose -start n -threads t"<<endl;
			return 0;
		}
	}

	// if (stdio)
	// 	RunStdioInput();
	// else
		RunNautyInput(start, n_threads);
}

