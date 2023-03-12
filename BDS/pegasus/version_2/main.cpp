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
// #include "src/nauty_reader.cpp"
#include "src/stdio_reader.cpp"
#include "src/main.h"

using namespace std;

/*
	Conventions
		- Graphs are 0-indexed and simple;
*/


/*
	I will define a global enviroment and model. Different matchings
	will be just different cost functions.
*/

GRBEnv env = GRBEnv(true);


/*
	This function receives a LP solution and retuns the edges of
	all st-cuts with capacity < 2..
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
	Tis function builds a fractional model

*/
void BuildFractional(GRBModel &frac_model, GRBVar *frac_vars){
	int n = countNodes(G);
	int m = countEdges(G);

	frac_model.set(GRB_IntParam_Method, 0); // Forcing Primal Simplex Method
	// Important, since fractional solution must be an Extreme Point
	
	
	// Setting the correct UB and OBJ.	
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

void BuildIntegral(GRBModel &int_model, GRBVar *int_vars){
	int n = countNodes(G);
	int m = countEdges(G);

	int_model.set(GRB_IntParam_LazyConstraints, 1); // Allow callback constraints
	// Set callback function
	MinimumCut cb = MinimumCut(int_vars, n, m);
	int_model.setCallback(&cb);

	// Setting the correct UB and OBJ.	
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
void IntegerSolution(ListGraph::EdgeMap<int> &IntSol, GRBModel &int_model, GRBVar *int_vars){
	for (ListGraph::EdgeIt e(G); e != INVALID; ++e) 
		IntSol[e] = -1;	

	try {
		int n = countNodes(G);
		int m = countEdges(G);

		GRBLinExpr obj = 0;
		for (ListGraph::EdgeIt e(G); e != INVALID; ++e)
			obj += cost[e] * int_vars[G.id(e)];

		int_model.setObjective(obj);

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

				// Sanity check
				// assert((node_u[id] == u) and (node_v[id]) == v);
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
void FractionalSolution(ListGraph::EdgeMap<double> &FracSol, GRBModel &frac_model, GRBVar *frac_vars){
	for (ListGraph::EdgeIt e(G); e != INVALID; ++e) 
		FracSol[e] = -1;	

	try {
		int n = countNodes(G);
		int m = countEdges(G);

		GRBLinExpr obj = 0;
		for (ListGraph::EdgeIt e(G); e != INVALID; ++e)
			obj += cost[e] * frac_vars[G.id(e)];

		frac_model.setObjective(obj);


		// Optimize model
		frac_model.optimize();

		bool found_feasible = 0;
		while (frac_model.get(GRB_IntAttr_SolCount) > 0 and !found_feasible){
		
			double *sol = frac_model.get(GRB_DoubleAttr_X, vars, m);
			// pair<double, vector<Edge> > min_cut = FindMinCut(sol, n, m);

			vector<GRBLinExpr> res = FindMinCuts(sol, vars, n, m);

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

					// assert((node_u[id] == u) and (node_v[id]) == v); // Sanity check
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



// void SolveModel(
// 	ListGraph::EdgeMap<double> &FracSol,
// 	ListGraph::EdgeMap<int> &IntSol){

// 	for (ListGraph::EdgeIt e(G); e != INVALID; ++e) 
// 		FracSol[e] = IntSol[e] = -1;	

// 	try {
// 		int n = countNodes(G);
// 		int m = countEdges(G);

// 		// Create an environment
// 		GRBEnv env = GRBEnv(true);

// 		// env.set("LogFile", "MAPSolver.log");  // Output may be large
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
// 			// pair<double, vector<Edge> > min_cut = FindMinCut(sol, n, m);

// 			vector<GRBLinExpr> res = FindMinCuts(sol, vars, n, m);

// 			// If min_cut.fist < 2, need to add constraint
// 			if (!res.empty()) { // Min cut < 2
// 				for (GRBLinExpr expr : res)
// 					model.addConstr(expr >= 2);
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

// 		if (sign(FracSol[G.edgeFromId(0)]) == -1) // No solution found
// 			return;

// 		/*
// 			Found a fractional extreme point solution.
// 			Now we will solve the integral version.

// 			We will warmstart the integer version by using the
// 			same model, just change the variable type
// 		*/


// 		// If relaxed formulation is integral, no need to resolve
// 		bool is_integral = 1;
// 		for (ListGraph::EdgeIt e(G); e != INVALID; ++e)
// 			if ((sign(FracSol[e]) != 0) and (sign(FracSol[e] - 1.0) != 0))
// 				is_integral = 0;

// 		if (is_integral){
// 			for (ListGraph::EdgeIt e(G); e != INVALID; ++e)
// 				IntSol[e] = FracSol[e];
// 			return;
// 		}



// 		for (int i = 0; i < m; i++) // Changing variables to binary
// 			vars[i].set(GRB_CharAttr_VType, GRB_BINARY);

// 		model.set(GRB_IntParam_LazyConstraints, 1); // Allow callback constraints

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


/*
	Wrapper function that call the solvers.
*/
void SolveMapInstance(
	ListGraph::EdgeMap<double> &FracSol,
	ListGraph::EdgeMap<int> &IntSol,
	ListGraph::EdgeMap<bool> &BDSSol,
	GRBModel &frac_model,
	GRBModel &int_model){

	FractionalSolution(FracSol, frac_model, frac_vars);

	if (sign(FracSol[G.edgeFromId(0)]) == -1)
		return;

	bool is_integral = 1;
	for (ListGraph::EdgeIt e(G); e != INVALID; ++e)
		if ((sign(FracSol[e]) != 0) and (sign(FracSol[e] - 1.0) != 0))
			is_integral = 0;

	if (is_integral){
		for (ListGraph::EdgeIt e(G); e != INVALID; ++e)
			IntSol[e] = FracSol[e];
		return;
	}

	IntegerSolution(IntSol, int_model, int_vars);

	if (IntSol[G.edgeFromId(0)] == -1)
		return;

	if (sign(FracSol[G.edgeFromId(0)]) == -1)
		return;

	BDSAlgorithm(FracSol, BDSSol);
}


bool __found_feasible;
int __cur_graph_id, __best_IP_graph_id, __best_IP_matching_id, __best_BDS_graph_id, __best_BDS_matching_id;
double __best_IP, __best_BDS;
ofstream g_out, log_out;

/*
	This function calls the LP, IP and BDS algorithms to
	solve the MAP problem and compares their outputs.
	If the output satisfies the requerements, it writes a
	file.
*/
void SolveCurrentMatching(int matching_id){
	ListGraph::EdgeMap<bool> BDSSol(G);
	ListGraph::EdgeMap<int> IntSol(G);
	ListGraph::EdgeMap<double> FracSol(G);

	SolveMapInstance(FracSol, IntSol, BDSSol);

	if (sign(FracSol[G.edgeFromId(0)]) == -1 or IntSol[G.edgeFromId(0)] == -1){
		log_out << "Exception on example g" << __cur_graph_id << " matching id " << matching_id << endl;
		return;
	}


	int cost_Int = 0;
	int cost_BDS = 0;
	double cost_Frac = 0;

	for (ListGraph::EdgeIt e(G); e != INVALID; ++e){
		int u = G.id(G.u(e));
		int v = G.id(G.v(e));

		cost_Int +=  IntSol[e] * cost[e];
		cost_Frac +=  FracSol[e] * cost[e];
		cost_BDS +=  (int)BDSSol[e] * cost[e];

		// cout<<u + 1<<' '<<v + 1<<' '<<FracSol[e]<<' '<<IntSol[e]<<' '<<BDSSol[e]<<endl;
	}

	/* 
		Found a feasible example, print to file
			- IP gap must be at least 6/5
			- BDS gap must be better than 4/3
	*/
	if (sign(5.0 * cost_Int - 6.0 * cost_Frac) >= 0 or sign(3.0 * cost_BDS - 4.0 * cost_Frac) > 0){
		if (__found_feasible == 0){ // First matching found for this graph
			// create file "g"+cnt
			g_out.open(to_string(countNodes(G)) + "/g" + to_string(__cur_graph_id));
			
			g_out << countNodes(G) <<' ' << countEdges(G) << endl << endl;
			
			for (ListGraph::EdgeIt e(G); e != INVALID; ++e)
				g_out << G.id(G.u(e)) << ' ' << G.id(G.v(e)) << endl;
			
			g_out << "----------" << endl << endl;
		}

		__found_feasible = 1;

		g_out << matching_id << ": ";

		bool first = 1;
		for (ListGraph::EdgeIt e(G); e != INVALID; ++e)
			if (cost[e] == 0){
				if (!first)
					g_out << ", ";

				g_out << G.id(G.u(e)) << " " << G.id(G.v(e));
				first = 0;
			}

		g_out<<endl;	

		g_out << "Frc: " << cost_Frac << " | ";
		for (ListGraph::EdgeIt e(G); e != INVALID; ++e)
			g_out << FracSol[e] << ' ';
		g_out << endl;

		g_out << "Int: " << cost_Int << " | ";
		for (ListGraph::EdgeIt e(G); e != INVALID; ++e)
			g_out << IntSol[e] << ' ';
		g_out << endl;

		g_out << "BDS: " << cost_BDS << " | ";
		for (ListGraph::EdgeIt e(G); e != INVALID; ++e)
			g_out << BDSSol[e] << ' ';
		g_out << endl << endl;

		// Generate entry in the log file
		log_out << "Found feasible example g" << __cur_graph_id << " matching id " << matching_id << endl;
		log_out << "Int/Frc = " << (double) cost_Int/cost_Frac << " BDS/Frc = " << (double) cost_BDS/cost_Frac << endl;
		log_out << endl;
	}

	if (sign((double)cost_Int/cost_Frac - __best_IP) > 0){
		__best_IP = (double)cost_Int/cost_Frac;
		__best_IP_graph_id = __cur_graph_id;
		__best_IP_matching_id = matching_id;
	}


	if (sign((double)cost_BDS/cost_Frac - __best_BDS) > 0){
		__best_BDS = (double)cost_BDS/cost_Frac;
		__best_BDS_graph_id = __cur_graph_id;
		__best_BDS_matching_id = matching_id;
	}

}


/*
	Backtracking algorithm to find all matchings of G.
	Matched edges are marked with cost 0 on the global
	EdgeMap cost. The running time is exponential
*/
void FindAllMatchings(int e_id, int &n, int &m, int &n_matched, int &total_matchings, 
	ListGraph::NodeMap<bool> &matched
	GRBModel &frac_model,
	GRBVars *frac_vars,
	GRBModel &int_model,
	GRBVars *int_vars){

	if (e_id >= m){
		SolveCurrentMatching(total_matchings, frac_model, frac_vars, int_model, int_vars);
		return;
	}

	if (n_matched >= n - 1){ // matching cant increase, prune
		SolveCurrentMatching(total_matchings, frac_model, frac_vars, int_model, int_vars);
		return;
	}


	// Case 1 : won't add edge e_id to the matching
	FindAllMatchings(e_id + 1, n, m, n_matched, total_matchings, matched, frac_model, frac_vars, int_model, int_vars); 

	// Case 2 : if possible, will add e_id to the matching
	ListGraph::Edge e = G.edgeFromId(e_id);
	if ((matched[G.u(e)] == 0) and (matched[G.v(e)] == 0)){ // May add e_id
		
		matched[G.u(e)] = 1;
		matched[G.v(e)] = 1;
		n_matched += 2;
		cost[e] = 0;
		total_matchings++;

		FindAllMatchings(e_id + 1, n, m, n_matched, total_matchings, matched,  frac_model, frac_vars, int_model, int_vars);

		matched[G.u(e)] = 0;
		matched[G.v(e)] = 0;
		n_matched -= 2;
		cost[e] = 1;
	}
}


/*
	Wrapper function for the matching backtrackig algorithm.
*/
void SolveAllMatchings(){
	int total_matchings = 1, n = countNodes(G), m = countEdges(G), n_matched = 0;

	// Initialize all edges to be heavy
	for (ListGraph::EdgeIt e(G); e != INVALID; ++e)
		cost[e] = 1;

	GRBModel frac_model(env);
	GRBVars frac_vars[m];
	BuildFractional(frac_model, frac_vars);

	GRBModel int_model(env);
	GRBVars int_vars[m];
	BuildIntegral(int_model, int_vars);


	ListGraph::NodeMap<bool> matched(G);

	FindAllMatchings(0, n, m, n_matched, total_matchings, matched, frac_model, frac_vars, int_model, int_vars);

	if (__found_feasible == 1)
		g_out << "Number of matchings : " << total_matchings << endl;
}



/*
	This functions receiv nauty's geng output from stdin(may modify this),
	build a LEMON graph and the log files, and calls the function
	that iterates through all matchings.
*/
void RunNautyInput(int start){
	
	// Start a global gurobi enviroment
	env.set(GRB_IntParam_OutputFlag, 0);
	env.start();


	__best_IP = __best_BDS = 1;
	__best_IP_graph_id = __best_IP_matching_id = __best_BDS_graph_id = __best_BDS_matching_id = 1;
	ofstream log_progress;

	int cnt = 0;
	while (readNautyGraph(G, cin)){
		cnt++;

		int n = countNodes(G);
		int m = countEdges(G);

		if (cnt < start) continue;


		if (cnt == 1 and start == 0){ // Create folder to log files, create log stream
			std::experimental::filesystem::create_directory("./" + to_string(n));
			log_out.open(to_string(n) + "/log"); // overwrite existing log
		}
		else if (start == cnt)
			log_out.open(to_string(n) + "/log", ios::app); // open in append mode

		// Next loop makes shure that lemon graph is consistent with the algorithm input
		int nvtx = n - 1, ok = 1;
		for (ListGraph::NodeIt v(G); v != INVALID; ++v)
			if (G.id(v) != nvtx--){
				cout << "Found inconsistency regarding vertex indexing -- Skip graph number " << cnt << endl;
				ok = 0;
				break;
			}
		// Next loop makes shure that lemon graph edge indexing is consistent
		set<int> q;
		for (ListGraph::EdgeIt e(G); e != INVALID; ++e)
			q.insert(G.id(e));

		assert(q.size() == m);
		assert(*q.rbegin() == m - 1);
		/* 
			Since the input data is massive, we will one write a file if
			there is some feasible solution.
		*/
		__found_feasible = 0;
		__cur_graph_id = cnt;

		SolveAllMatchings();

		if (__found_feasible)
			g_out.close();

		log_progress.open(to_string(n) + "/log_progress");
		log_progress << "Last fully processed graph g" << cnt << endl;
		log_progress << "Best IP/Frac: " << __best_IP << " g" << __best_IP_graph_id << " matching " << __best_IP_matching_id << endl;
		log_progress << "Best BDS/Frac: " << __best_BDS << " g" << __best_BDS_graph_id << " matching " << __best_BDS_matching_id << endl;
		log_progress.close();	
	}

	log_out.close();
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

