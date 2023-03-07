/*
	Matching Augmentation problem
		- Linear Relaxation
		- Integer Solution
		- BDS Approximation Algorithm 

	
	Author : Gabriel Morete
*/

#include "gurobi_c++.h"
#include "src/libraries_and_utils.h"
#include "src/lemon.h"
#include "src/bds.cpp"
#include "src/nauty_reader.cpp"
#include "src/stdio_reader.cpp"

using namespace std;

/*
	Conventions
		- Graphs are 0-indexed and simple;
*/

/*
	This function receives a LP solution and retuns the edges
	of a global minimum cut and its value.
*/
pair<double, vector<Edge> > FindMinCut(double *sol, int n, int m){
	ListGraph::EdgeMap<double> capacity(G);

	// Build EdgeMap of capacities
	for (ListGraph::EdgeIt e(G); e != INVALID; ++e)
		capacity[e] = sol[G.id(e)]; // Using current solution as capacity

	// Build Gomory-Hu Tree
	GomoryHu<ListGraph, ListGraph::EdgeMap<double> > GMH(G,capacity);
	GMH.run();	                

	ListGraph::Node min_v = G.nodeFromId(0);
	ListGraph::Node min_u = G.nodeFromId(1);
	for (ListGraph::NodeIt v(G); v != INVALID; ++v)
		for (ListGraph::NodeIt u(G); u != INVALID; ++u){
			if (GMH.minCutValue(v, u) < GMH.minCutValue(min_v, min_u)){
				min_v = v;
				min_u = u;
			}
	}

	vector<Edge> min_cut;
	for (GomoryHu<ListGraph, ListGraph::EdgeMap<double> >::MinCutEdgeIt e(GMH, min_v, min_u); e != INVALID; ++e)
		min_cut.push_back(e);
	
	return make_pair(GMH.minCutValue(min_v, min_u), min_cut);
}

/*
	This is the separator function for the MIP. 
	If the solution is not a 2ECSS it adds a cut
	separating one 2ECC.
*/
class MinimumCut: public GRBCallback {
	public:
		GRBVar* vars;
		int n, m;
		MinimumCut(GRBVar* xvars, int xn, int xm){
			vars = xvars;
			n = xn;
			m = xm;
		}
	protected:
		void callback(){
			try {
				if (where == GRB_CB_MIPSOL){
					// Solver found an integral optimal solution for the
					// current formulation, must check if if there is a
					// bridge or a cut.

					double *x = getSolution(vars, m);

					ListGraph::EdgeMap<bool> in_sol(G);
					
					for (ListGraph::EdgeIt e(G); e != INVALID; ++e){
						int id = G.id(e);
						int u = G.id(G.u(e));
						int v = G.id(G.v(e));

						if (x[id] > 0.5)
							in_sol[e] = 1;
					}

					ListGraph::NodeMap<bool> ones(G, 1); // Must be spanning
					// H is a spanning subgraph with all edges in the solution
					SubGraph<ListGraph> H(G, ones, in_sol);

					ListGraph::NodeMap<int> ebcc(G, -1);
					biEdgeConnectedComponents(H, ebcc);

					int ncmp = 0;
					for (ListGraph::NodeIt v(G); v != INVALID; ++v)
						ncmp = max(ncmp, ebcc[v]);

					if (ncmp > 0){ // Not 2ECSS, must add a cut 
	
						// The cut will be all edges crossing the cut of the ebcc with id 0
						GRBLinExpr expr = 0;
	
						for (ListGraph::EdgeIt e(G); e != INVALID; ++e){
							ListGraph::Node u = G.u(e);
							ListGraph::Node v = G.v(e);
							
							if (ebcc[u] == 0 and ebcc[v] != 0)
								expr += vars[G.id(e)];
					
							if (ebcc[v] == 0 and ebcc[u] != 0)
								expr += vars[G.id(e)];
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
};

/*
	This function returns a optimum integer solution to MAP.
	If no solution is found, it returns a all -1 edge map.
*/
void IntegerSolution(ListGraph::EdgeMap<int> &IntSol){

	for (ListGraph::EdgeIt e(G); e != INVALID; ++e) 
		IntSol[e] = -1;	

	try {
		int n = countNodes(G);
		int m = countEdges(G);

		// Create an environment
		GRBEnv env = GRBEnv(true);
		// env.set("LogFile", "MAPInteger.log"); // Output may be large
		env.set(GRB_IntParam_OutputFlag, 0);
		env.start();

		// Create an empty model
		GRBModel model = GRBModel(env);
		model.set(GRB_IntParam_LazyConstraints, 1); // Allow callback constraints
		
		GRBVar vars[m];
		int node_u[m], node_v[m];

		// Setting the correct UB and OBJ.	
		for (ListGraph::EdgeIt e(G); e != INVALID; ++e){
			int u = G.id(G.u(e));
			int v = G.id(G.v(e));
			vars[G.id(e)] = model.addVar(0.0, 1.0, cost[e], GRB_BINARY, "x_" + to_string(u) + "_" + to_string(v));
			node_u[G.id(e)] = u;
			node_v[G.id(e)] = v;
		}

		// Add \delta(v) >= 2, constraints
		GRBLinExpr deg2[n];
		for (ListGraph::EdgeIt e(G); e != INVALID; ++e){
			int id = G.id(e);
			int u = G.id(G.u(e));
			int v = G.id(G.v(e));

			deg2[u] += vars[id];
			deg2[v] += vars[id];

		}

		for (int v = 0; v < n; v++)
			model.addConstr(deg2[v] >= 2, "deg2_" + to_string(v));

		// Set callback function
    	MinimumCut cb = MinimumCut(vars, n, m);
    	model.setCallback(&cb);
		
		// Optimize model
		model.optimize();

		// Found optimal solution
		if (model.get(GRB_IntAttr_SolCount) > 0){
			// cout << "Obj: " << model.get(GRB_DoubleAttr_ObjVal) << endl;

			double *sol = model.get(GRB_DoubleAttr_X, vars, m);

			for (ListGraph::EdgeIt e(G); e != INVALID; ++e){
				int id = G.id(e);
				int u = G.id(G.u(e));
				int v = G.id(G.v(e));

				if (sol[id] > 0.5)
					IntSol[e] = 1;
				else
					IntSol[e] = 0;
				// cout<<u + 1<<' '<<v + 1<<' '<<abs(sol[u][v])<<endl;

				// Sanity check
				assert((node_u[id] == u) and (node_v[id]) == v);
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
	This function returns a optimum fractional solution to MAP.
	If no solution is found, it returns a all -1 edge map.
*/
void FractionalSolution(ListGraph::EdgeMap<double> &FracSol){
	for (ListGraph::EdgeIt e(G); e != INVALID; ++e) 
		FracSol[e] = -1;	

	try {
		int n = countNodes(G);
		int m = countEdges(G);

		// Create an environment
		GRBEnv env = GRBEnv(true);
		// env.set("LogFile", "MAPFractional.log");  // Output may be large
		env.set(GRB_IntParam_OutputFlag, 0);
		env.start();

		// Create an empty model
		GRBModel model = GRBModel(env);
		model.set(GRB_IntParam_Method, 0); // Forcing Primal Simplex Method
		// Important, since fractional solution must be an Extreme Point
		
		
		GRBVar vars[m];
		int node_u[m], node_v[m];

		// Setting the correct UB and OBJ.	
		for (ListGraph::EdgeIt e(G); e != INVALID; ++e){
			int u = G.id(G.u(e));
			int v = G.id(G.v(e));
			vars[G.id(e)] = model.addVar(0.0, 1.0, cost[e], GRB_CONTINUOUS, "x_" + to_string(u) + "_" + to_string(v));
			node_u[G.id(e)] = u;
			node_v[G.id(e)] = v;
		}

		// Add \delta(v) >= 2, constraints
		GRBLinExpr deg2[n];
		for (ListGraph::EdgeIt e(G); e != INVALID; ++e){
			int id = G.id(e);
			int u = G.id(G.u(e));
			int v = G.id(G.v(e));

			deg2[u] += vars[id];
			deg2[v] += vars[id];

		}

		for (int v = 0; v < n; v++)
			model.addConstr(deg2[v] >= 2, "deg2_" + to_string(v));

		// Optimize model
		model.optimize();

		bool found_feasible = 0;
		while (model.get(GRB_IntAttr_SolCount) > 0 and !found_feasible){
		
			double *sol = model.get(GRB_DoubleAttr_X, vars, m);
			pair<double, vector<Edge> > min_cut = FindMinCut(sol, n, m);

			// If min_cut.fist < 2, need to add constraint
			if (sign(min_cut.first - 2.0) < 0) { // Min cut < 2
				GRBLinExpr expr = 0;
				
				// All all edges of the global min cut as a contraint.
				for (Edge e : min_cut.second){
					int id = G.id(e);
					int u = G.id(G.u(e));
					int v = G.id(G.v(e));				
					
					expr += vars[id];
					assert((node_u[id] == u) and (node_v[id]) == v); // Sanity check
				}

				model.addConstr(expr >= 2);
				model.optimize();
			}
			else{
				// Found a feasible opt
				// cout << "Obj: " << model.get(GRB_DoubleAttr_ObjVal) << endl;

				for (ListGraph::EdgeIt e(G); e != INVALID; ++e){
					int id = G.id(e);
					int u = G.id(G.u(e));
					int v = G.id(G.v(e));

					FracSol[e] = sol[id]; 

					// cout<<u + 1<<' '<<v + 1<<' '<<abs(sol[u][v])<<endl;

					assert((node_u[id] == u) and (node_v[id]) == v); // Sanity check
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

void SolveModel(
	ListGraph::EdgeMap<double> &FracSol,
	ListGraph::EdgeMap<int> &IntSol){

	for (ListGraph::EdgeIt e(G); e != INVALID; ++e) 
		FracSol[e] = IntSol[e] = -1;	

	try {
		int n = countNodes(G);
		int m = countEdges(G);

		// Create an environment
		GRBEnv env = GRBEnv(true);
		// env.set("LogFile", "MAPFractional.log");  // Output may be large
		env.set(GRB_IntParam_OutputFlag, 0);
		env.start();

		// Create an empty model
		GRBModel model = GRBModel(env);
		model.set(GRB_IntParam_Method, 0); // Forcing Primal Simplex Method
		// Important, since fractional solution must be an Extreme Point
		
		
		GRBVar vars[m];
		int node_u[m], node_v[m];

		// Setting the correct UB and OBJ.	
		for (ListGraph::EdgeIt e(G); e != INVALID; ++e){
			int u = G.id(G.u(e));
			int v = G.id(G.v(e));
			vars[G.id(e)] = model.addVar(0.0, 1.0, cost[e], GRB_CONTINUOUS, "x_" + to_string(u) + "_" + to_string(v));
			node_u[G.id(e)] = u;
			node_v[G.id(e)] = v;
		}

		// Add \delta(v) >= 2, constraints
		GRBLinExpr deg2[n];
		for (ListGraph::EdgeIt e(G); e != INVALID; ++e){
			int id = G.id(e);
			int u = G.id(G.u(e));
			int v = G.id(G.v(e));

			deg2[u] += vars[id];
			deg2[v] += vars[id];

		}

		for (int v = 0; v < n; v++)
			model.addConstr(deg2[v] >= 2, "deg2_" + to_string(v));

		// Optimize model
		model.optimize();

		bool found_feasible = 0;
		while (model.get(GRB_IntAttr_SolCount) > 0 and !found_feasible){
		
			double *sol = model.get(GRB_DoubleAttr_X, vars, m);
			pair<double, vector<Edge> > min_cut = FindMinCut(sol, n, m);

			// If min_cut.fist < 2, need to add constraint
			if (sign(min_cut.first - 2.0) < 0) { // Min cut < 2
				GRBLinExpr expr = 0;
				
				// All all edges of the global min cut as a contraint.
				for (Edge e : min_cut.second){
					int id = G.id(e);
					int u = G.id(G.u(e));
					int v = G.id(G.v(e));				
					
					expr += vars[id];
					assert((node_u[id] == u) and (node_v[id]) == v); // Sanity check
				}

				model.addConstr(expr >= 2);
				model.optimize();
			}
			else{
				// Found a feasible opt
				// cout << "Obj: " << model.get(GRB_DoubleAttr_ObjVal) << endl;

				for (ListGraph::EdgeIt e(G); e != INVALID; ++e){
					int id = G.id(e);
					int u = G.id(G.u(e));
					int v = G.id(G.v(e));

					FracSol[e] = sol[id]; 

					// cout<<u + 1<<' '<<v + 1<<' '<<abs(sol[u][v])<<endl;

					assert((node_u[id] == u) and (node_v[id]) == v); // Sanity check
				}
			
				found_feasible = 1;
			}

			delete[] sol;
		}

	// Now solve the IP
	for (int i = 0; i < m; i++)
		vars[i].set(GRB_CharAttr_VType, GRB_BINARY);

	model.set(GRB_IntParam_LazyConstraints, 1); // Allow callback constraints

			// Set callback function
    	MinimumCut cb = MinimumCut(vars, n, m);
    	model.setCallback(&cb);
		
		// Optimize model
		model.optimize();

		// Found optimal solution
		if (model.get(GRB_IntAttr_SolCount) > 0){
			// cout << "Obj: " << model.get(GRB_DoubleAttr_ObjVal) << endl;

			double *sol = model.get(GRB_DoubleAttr_X, vars, m);

			for (ListGraph::EdgeIt e(G); e != INVALID; ++e){
				int id = G.id(e);
				int u = G.id(G.u(e));
				int v = G.id(G.v(e));

				if (sol[id] > 0.5)
					IntSol[e] = 1;
				else
					IntSol[e] = 0;
				// cout<<u + 1<<' '<<v + 1<<' '<<abs(sol[u][v])<<endl;

				// Sanity check
				assert((node_u[id] == u) and (node_v[id]) == v);
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
	ListGraph::EdgeMap<int> &IntSol,
	ListGraph::EdgeMap<bool> &BDSSol){

	// int tries_cnt = 0;

	// do {
		FractionalSolution(FracSol);
	// } while (tries_cnt < 3 and sign(FracSol[G.edgeFromId(0)]) == -1);

	// if (sign(FracSol[G.edgeFromId(0)]) == -1)
	// 	return;

	// tries_cnt = 0;

	// do {
		IntegerSolution(IntSol);
	// } while (tries_cnt < 3 and IntSol[G.edgeFromId(0)] == -1);

	// if (IntSol[G.edgeFromId(0)] == -1)
	// 	return;

	// SolveModel(FracSol, IntSol);

	BDSAlgorithm(FracSol, BDSSol);
}


signed main(int argc, char *argv[]){
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

	if (stdio)
		RunStdioInput();
	else
		RunNautyInput(start);
}

