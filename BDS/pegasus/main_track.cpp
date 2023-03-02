/*
	Matching Augmentation problem pack
		- Linear Relaxation
		- Integer Solution
		- BDS Approximation Algorithm 

	
	Author : Gabriel Morete
*/


#include <iostream>
#include <fstream>
#include <experimental/filesystem>
#include <vector>
#include <array>
#include <cassert>
#include <string>
#include <set>
#include <algorithm>
#include "gurobi_c++.h"
#include <lemon/list_graph.h>
#include <lemon/gomory_hu.h>
#include <lemon/adaptors.h>
#include <lemon/connectivity.h>
#include <lemon/nauty_reader.h>


using namespace std;
using namespace lemon;

typedef ListGraph::Node Node;
typedef ListGraph::Edge Edge;
typedef ListGraph::NodeIt NodeIt;
typedef ListGraph::EdgeIt EdgeIt;
// typedef ListGraph::NodeMap<int> NodeMap<int>;
// typedef ListGraph::NodeMap<double> NodeMap<double>;
// typedef ListGraph::EdgeMap<int> EdgeMap<int>;
// typedef ListGraph::EdgeMap<double> EdgeMap<double>;


// Safe handling doubles
const double EPS = 1e-3;
int sign(double x) { return (x > EPS) - (x < -EPS); }


/*
	Conventions
		- Graphs are 0-indexed and simple;
*/


ListGraph G; // Declare global Graph
ListGraph::EdgeMap<int> cost(G); // Cost of the edges




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
		env.set("LogFile", "MAPInteger.log");
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
		env.set("LogFile", "MAPFractional.log");
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
		}


	} catch(GRBException e) {
		cout << "Error code = " << e.getErrorCode() << endl;
		cout << e.getMessage() << endl;
	} catch(...) {
		cout << "Exception during optimization" << endl;
	}
}


/*
	DFS Step of BDS Algorithm. The next tree edge is chosen by the following criteria.
		- Matching edge
		- Heavy edge that maximizes x^*_e
	
	O(n^2)
*/
void BDSDFS(ListGraph::Node v, 
	ListGraph::NodeMap<ListGraph::Node> &parent, 
	ListGraph::EdgeMap<double> &FracSol, 
	ListGraph::EdgeMap<bool> &BDSSol,
	int &clk,
	ListGraph::NodeMap<int> &in,
	ListGraph::NodeMap<int> &out){

	in[v] = clk++;
	
	bool found_next = 0;
	do {
		ListGraph::Arc next_arc;

		found_next = 0;

		// Here we iterate through the arcs leaving node v.
		for (ListGraph::OutArcIt a(G, v); a != INVALID; ++a){
			ListGraph::Node u = G.target(a);

			if (parent[u] == u and G.id(u) != 0){ // Unvisited non root node

				if ((!found_next) or (cost[a] == 0)){
					found_next = 1;
					next_arc = a;
				}
				
				else if (cost[next_arc] == 1 and FracSol[next_arc] < FracSol[a])
						next_arc = a;
			}
		}

		if (found_next){
			BDSSol[next_arc] = 1;
			
			ListGraph::Node u = G.target(next_arc);

			parent[u] = v;

			BDSDFS(u, parent, FracSol, BDSSol, clk, in, out);
		}

	} while (found_next);

	out[v] = clk++;
}

/*
	Retuns true if v is a proper descendent of u.
*/
bool StrictDec(Node u, Node v, ListGraph::NodeMap<int> &in, ListGraph::NodeMap<int> &out){ // is v strict dec of u
	return (in[u] < in[v]) and (out[v] < out[u]);
}

/*
	Retuns true if v is a descendent of u.
*/
bool Dec(Node u, Node v, ListGraph::NodeMap<int> &in, ListGraph::NodeMap<int> &out){ // is v strict dec of u
	return (in[u] <= in[v]) and (out[v] <= out[u]);
}

/*
	Dynamic Programming algorithm to solving the up-link only 
	agumentation problem.

	O(n^2|L|)
*/
int UpLinkDP(Node v, ListGraph::NodeMap<int> &memo, 
	ListGraph::NodeMap<Edge> &dp_edge, 
	ListGraph::NodeMap<ListGraph::Node> &parent, 
	ListGraph::NodeMap<int> &in, 
	ListGraph::NodeMap<int> &out,
	ListGraph::EdgeMap<bool> &BDSSol,
	SubGraph<ListGraph> &T){
	
	if (memo[v] != -1)
		return memo[v];

	memo[v] = countEdges(G) + 1; // infinity

	for (EdgeIt e(G); e != INVALID; ++e)
		if (BDSSol[e] == 0){ // A backedge
			
			Node u = G.u(e);
			Node w = G.v(e);
			if (Dec(w, u, in, out)) // If u is a descendent of w
				swap(w, u);

			// u is the ancestor node, w is the descendent

			if (StrictDec(u, v, in, out) and Dec(v, w, in, out)){ // feasible link

				// cout<<"considering edge "<<G.id(u) + 1<<' '<<G.id(w) + 1<<" to cover node "<<G.id(v) + 1<<endl;

				int subtree_cost = 0;

				Node h = w, lst = w, p = parent[v];
				while (h != p){
					for (SubGraph<ListGraph>::OutArcIt a(T, h); a != INVALID; ++a){
						ListGraph::Node y = T.target(a);

						if (y != lst and y != parent[h]){
							// cout<<" Transition "<<G.id(v) + 1<<' '<<G.id(y) + 1<<endl;
							subtree_cost += UpLinkDP(y, memo, dp_edge, parent, in, out, BDSSol, T);
						}
					}	

					lst = h;
					h = parent[h];
				}

				if (subtree_cost + cost[e] < memo[v]){
					// cout<<"Foud New Solution to "<<G.id(v) + 1<<' '<<G.id(u) + 1<<' '<<G.id(w) + 1<<endl;
					memo[v] = subtree_cost + cost[e];
					dp_edge[v] = e;
				}
			}
		}

	return memo[v];
}

/*
	Algorithm to recover the uplink DP solution.

	O(n + m)
*/
void RecoverUpLinkSol(Node v,
	ListGraph::NodeMap<Edge> &dp_edge, 
	ListGraph::NodeMap<ListGraph::Node> &parent, 
	ListGraph::NodeMap<int> &in, 
	ListGraph::NodeMap<int> &out,
	ListGraph::EdgeMap<bool> &BDSSol,
	SubGraph<ListGraph> &T){


	Edge e = dp_edge[v];
	BDSSol[e] = 1;

	Node u = G.u(e);
	Node w = G.v(e);
	if (Dec(w, u, in, out)) // If u is a decendent of w
		swap(w, u);

	// u is the ancestor node, w is the descendent

	// cout<<"Recovering Solution Node "<<G.id(v) + 1<<' '<<" Edge "<<G.id(u) + 1<<' '<<G.id(w) + 1<<endl;

	Node p = parent[v], lst = w;
	while (w != p){
		for (SubGraph<ListGraph>::OutArcIt a(T, w); a != INVALID; ++a){
			ListGraph::Node y = T.target(a);

			if (y != lst and y != parent[w])
				RecoverUpLinkSol(y, dp_edge, parent, in, out, BDSSol, T);
		}
			
		lst = w;
		w = parent[w];
	}
}

void UpLinkAugmentation(
	SubGraph<ListGraph> &T,
	ListGraph::EdgeMap<bool> &BDSSol,
	ListGraph::NodeMap<ListGraph::Node> &parent, 
	ListGraph::NodeMap<int> &in,
	ListGraph::NodeMap<int> &out){

	ListGraph::NodeMap<int> memo(G, -1); 
	ListGraph::NodeMap<Edge> dp_edge(G); 

	for (SubGraph<ListGraph>::OutArcIt a(T, G.nodeFromId(0)); a != INVALID; ++a){
		UpLinkDP(T.target(a), memo, dp_edge, parent, in, out, BDSSol, T);
		RecoverUpLinkSol(T.target(a), dp_edge, parent, in, out, BDSSol, T);
	}
}

/*
	Implementation of the BDS MAP algorithm.
	Must receive a extreme point solution of the cut LP.

	O(n^2|L|)
*/
void BDSAlgorithm(ListGraph::EdgeMap<double> &FracSol, ListGraph::EdgeMap<bool> &BDSSol){
	int n = countNodes(G);

	// Step 1, find a DFS Tree
	ListGraph::NodeMap<ListGraph::Node> parent(G);
	ListGraph::NodeMap<int> in(G), out(G); 

	for (ListGraph::NodeIt v(G); v != INVALID; ++v)
		parent[v] = v;

	int cnt = 0;
	BDSDFS(G.nodeFromId(0), parent, FracSol, BDSSol, cnt, in, out);

	// cout<<"BDS tree"<<endl;
	// for (NodeIt v(G); v != INVALID; ++v){
	// 	cout<<parent[v] + 1<<" <-- "<<G.id(v) + 1<<endl;
	// }




	ListGraph::NodeMap<bool> ones(G, 1); // Subgraph must be spanning

	ListGraph::EdgeMap<bool> tree_edges(G); // Graph adaptor changes with edge map
	for (ListGraph::EdgeIt e(G); e != INVALID; ++e)
		tree_edges[e] = BDSSol[e];

	// Build DFS Tree
	SubGraph<ListGraph> T(G, ones, tree_edges);
	
	assert(biEdgeConnected(G) == 1);
	assert(connected(T) == 1);


	// Step 2, uplink only augmentation
	UpLinkAugmentation(T, BDSSol, parent, in, out);

	// Sanity check, checks is BDS returned a feasible solution
	SubGraph<ListGraph> H(G, ones, BDSSol);
	assert(biEdgeConnected(H) == 1);
}


/*
	Nauty Reader
*/


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
	ListGraph::EdgeMap<int> IntSol(G);
	IntegerSolution(IntSol);

	ListGraph::EdgeMap<double> FracSol(G);
	FractionalSolution(FracSol);

	ListGraph::EdgeMap<bool> BDSSol(G);
	BDSAlgorithm(FracSol, BDSSol);

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
			- IP gap must be at least 7/5
			- BDS gap must be better than 4/3
	*/
	if (sign(5.0 * cost_Int - 7.0 * cost_Frac) >= 0 or sign(3.0 * cost_BDS - 4.0 * cost_Frac) > 0){
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
void FindAllMatchings(int e_id, int &n, int &m, int &n_matched, int &total_matchings, ListGraph::NodeMap<bool> &matched){
	if (e_id >= m){
		SolveCurrentMatching(total_matchings);
		return;
	}

	if (n_matched >= n - 1){ // matching cant increase, prune
		SolveCurrentMatching(total_matchings);
		return;
	}


	// Case 1 : won't add edge e_id to the matching
	FindAllMatchings(e_id + 1, n, m, n_matched, total_matchings, matched); 

	// Case 2 : if possible, will add e_id to the matching
	ListGraph::Edge e = G.edgeFromId(e_id);
	if ((matched[G.u(e)] == 0) and (matched[G.v(e)] == 0)){ // May add e_id
		
		matched[G.u(e)] = 1;
		matched[G.v(e)] = 1;
		n_matched += 2;
		cost[e] = 0;
		total_matchings++;

		FindAllMatchings(e_id + 1, n, m, n_matched, total_matchings, matched);

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
	// Initialize all edges to be heavy
	for (ListGraph::EdgeIt e(G); e != INVALID; ++e)
		cost[e] = 1;

	ListGraph::NodeMap<bool> matched(G);

	int total_matchings = 1, n = countNodes(G), m = countEdges(G), n_matched = 0;
	FindAllMatchings(0, n, m, n_matched, total_matchings, matched);

	if (__found_feasible == 1)
		g_out << "Number of matchings : " << total_matchings << endl;
}



/*
	This functions receiv nauty's geng output from stdin(may modify this),
	build a LEMON graph and the log files, and calls the function
	that iterates through all matchings.
*/
void RunNautyInput(){
	__best_IP = __best_BDS = 1;
	__best_IP_graph_id = __best_IP_matching_id = __best_BDS_graph_id = __best_BDS_matching_id = 1;
	ofstream log_progress;

	int cnt = 1;
	while (readNautyGraph(G, cin)){
		int n = countNodes(G);
		int m = countEdges(G);

		if (cnt == 1){ // Create folder to log files, create log stream
			std::experimental::filesystem::create_directory("./" + to_string(n));
			log_out.open(to_string(n) + "/log");
		}	

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
		cnt++;

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


signed main(){
	ListGraph G;

	RunNautyInput();

}
