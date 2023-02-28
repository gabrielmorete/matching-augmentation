/*
	Matching Augmentation problem pack
		- Linear Relaxation
		- Integer Solution
		- BDS Approximation Algorithm 

	
	Author : Gabriel Morete
*/


#include <iostream>
#include <fstream>
#include <filesystem>
#include <vector>
#include <array>
#include <cassert>
#include <string>
#include "/Library/gurobi1001/macos_universal2/include/gurobi_c++.h"
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


// I will also save the input information
int _n, _m;
vector< array<int, 3> > _edges;


/*
	This function reads the Graph from Stdio. Graph is 1-indexed
	The input format will be
		n m       number of nodes, edges
		a_1 b_1 c_1     edge bertween a_1, b_1 with cost c_1
		...
		a_m b_m c_m 
*/
void ReadStdioInput(){
	int n, m;
	cin>>n>>m;
	
	assert(n >= 3);
	assert(m >= n);

	_n = n;
	_m = m;

	for (int i = 0; i < n; i++){
		ListGraph::Node v = G.addNode();
		if (G.id(v) != i)
			cout<<"Error : vertex don't match id"<<endl;
		assert(G.id(v) == i);
	}

	for (int i = 0; i < m; i++){
		int a, b, c;
		cin>>a>>b>>c;

		a--; // 0-indexed
		b--; 

		// _edges.push_back({a, b, c});

		ListGraph::Edge e = G.addEdge(G.nodeFromId(a), G.nodeFromId(b));
		cost[e] = c;
	}
}

/*
	This function receives a LP solution and retuns the edges
	of a global minimum cut and its value.
*/
pair<double, vector<Edge> > FindMinCut(double **sol){
	int n = countNodes(G);
	ListGraph::EdgeMap<double> capacity(G);

	// Build EdgeMap of capacities
	for (ListGraph::EdgeIt e(G); e != INVALID; ++e){
		int u = G.id(G.u(e));
		int v = G.id(G.v(e));

		capacity[e] = sol[v][u]; // Using current solution as capacity
	}

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
	If there is a cut with value < 2 we will add this cut to the constraints.
*/
class MinimumCut: public GRBCallback {
	public:
		GRBVar** vars;
		int n;
		MinimumCut(GRBVar** xvars, int xn){
			vars = xvars;
			n = xn;
		}
	protected:
		void callback(){
			try {
				if (where == GRB_CB_MIPSOL){
					// Solver found an integral optimal solution for the
					// current formulation, must check if if there is a
					// minimum cut with value < 2
					// Since the solution is integral, one could just
					// search for bridges.

					double *x[n];
					for (int i = 0; i < n; i++)
						x[i] = getSolution(vars[i], n);

					ListGraph::EdgeMap<double> capacity(G);
					
					for (ListGraph::EdgeIt e(G); e != INVALID; ++e){
						int u = G.id(G.u(e));
						int v = G.id(G.v(e));

						capacity[e] = x[v][u]; // Using LP Value as capacities
					}


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

					// Vertices min_u, min_v have the minimum st-cut

					// cout<<"Min cut val"<< GMH.minCutValue(min_v, min_u) << endl;


					if (sign(GMH.minCutValue(min_v, min_u) - 2.0) < 0) { // Min cut < 2
						GRBLinExpr expr = 0;
						
						// All all edges of the global min cut as a contraint.
						for (GomoryHu<ListGraph, ListGraph::EdgeMap<double> >::MinCutEdgeIt e(GMH, min_v, min_u); e != INVALID; ++e){
							int u = G.id(G.u(e));
							int v = G.id(G.v(e));				
							
							expr += vars[u][v];
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
	This is the separator function for the MIP. 
	If the solution is not a 2ECSS it adds a cut
*/
// class MinimumCut: public GRBCallback {
// 	public:
// 		GRBVar** vars;
// 		int n;
// 		MinimumCut(GRBVar** xvars, int xn){
// 			vars = xvars;
// 			n = xn;
// 		}
// 	protected:
// 		void callback(){
// 			try {
// 				if (where == GRB_CB_MIPSOL){
// 					// Solver found an integral optimal solution for the
// 					// current formulation, must check if if there is a
// 					// bridge or a cut.

// 					double *x[n];
// 					for (int i = 0; i < n; i++)
// 						x[i] = getSolution(vars[i], n);

// 					ListGraph::EdgeMap<bool> in_sol(G);
					
// 					for (ListGraph::EdgeIt e(G); e != INVALID; ++e){
// 						int u = G.id(G.u(e));
// 						int v = G.id(G.v(e));

// 						if (x[v][u] > 0.5)
// 							in_sol[e] = 1;
// 					}

// 					ListGraph::NodeMap<bool> ones(G, 1);
// 					ListGraph H = SubGraph(G, ones, in_sol);
// 					// H is a spanning subgraph with all edges in the solution

// 					if (biEdgeConnected(H) == 0){ // Not 2ECSS, must add a cut 
// 						ListGraph::NodeMap<int> ebcc(G);
// 						biEdgeConnectedComponents(G, ebcc);

// 						GRBLinExpr expr = 0;
	
// 						for (ListGraph::EdgeIt e(G); e != INVALID; ++e){
// 							ListGraph::Node u = G.u(e);
// 							ListGraph::Node v = G.v(e);
							
// 							if (ebcc[u] == 0 and ebcc[v] != 0)
// 								expr += vars[G.id(u)][G.id(v)];
					
// 							if (ebcc[v] == 0 and ebcc[u] != 0)
// 								expr += vars[G.id(v)][G.id(u)];
// 						}

// 						addLazy(expr >= 2);
// 					}
// 			} 
// 			catch (GRBException e){
// 				cout << "Error number: " << e.getErrorCode() << endl;
// 				cout << e.getMessage() << endl;
// 			} 
// 			catch (...){
// 				cout << "Error during callback" << endl;
// 			}
// 		}
// };



/*
	This function returns a optimum integer solution to MAP.
	If no solution is found, it returns a all -1 edge map.
*/
void IntegerSolution(ListGraph::EdgeMap<int> &IntSol){

	for (ListGraph::EdgeIt e(G); e != INVALID; ++e) 
		IntSol[e] = -1;	

	try {
		int n = countNodes(G);

		// Create an environment
		GRBEnv env = GRBEnv(true);
		env.set("LogFile", "MAPInteger.log");
		env.set(GRB_IntParam_OutputFlag, 0);
		env.start();

		// Create an empty model
		GRBModel model = GRBModel(env);
		model.set(GRB_IntParam_LazyConstraints, 1); // Allow callback constraints
		
		GRBVar **vars = NULL;
		vars = new GRBVar*[n];
		for (int i = 0; i < n; i++)
			vars[i] = new GRBVar[n];
		

		// Create all variables. Optimize here
		for (int i = 0; i < n; i++)
			for (int j = 0; j <= i; j++){
				vars[i][j] = model.addVar(0.0, 0.0, 0.0, GRB_BINARY, "x_" + to_string(i) + "_" + to_string(j));
				vars[j][i] = vars[i][j];
			}    

		// Setting the correct UB and OBJ.	
		for (ListGraph::EdgeIt e(G); e != INVALID; ++e){
			int u = G.id(G.u(e));
			int v = G.id(G.v(e));

			vars[u][v].set(GRB_DoubleAttr_UB, 1);
			vars[u][v].set(GRB_DoubleAttr_Obj, cost[e]);
		}

		// Add \delta(v) >= 2, constraints
		for (int i = 0; i < n; i++){
			
			GRBLinExpr expr = 0;
			for (int j = 0; j < n; j++)
				expr += vars[i][j];
			
			model.addConstr(expr >= 2, "cut2_" + to_string(i));
		}

		// Set callback function
    	MinimumCut cb = MinimumCut(vars, n);
    	model.setCallback(&cb);
		
		// Optimize model
		model.optimize();

		// Found optimal solution
		if (model.get(GRB_IntAttr_SolCount) > 0){
			// cout << "Obj: " << model.get(GRB_DoubleAttr_ObjVal) << endl;

			double **sol = new double*[n];
			for (int i = 0; i < n; i++)
				sol[i] = model.get(GRB_DoubleAttr_X, vars[i], n);

			for (ListGraph::EdgeIt e(G); e != INVALID; ++e){
				int u = G.id(G.u(e));
				int v = G.id(G.v(e));

				if (sol[u][v] > 0.5)
					IntSol[e] = 1;
				else
					IntSol[e] = 0;

				// cout<<u + 1<<' '<<v + 1<<' '<<abs(sol[u][v])<<endl;
			}

			for (int i = 0; i < n; i++)
				delete[] sol[i];
			delete[] sol;
		}

		for (int i = 0; i < n; i++)
			delete[] vars[i];
		delete[] vars;

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

		// Create an environment
		GRBEnv env = GRBEnv(true);
		env.set("LogFile", "MAPFractional.log");
		env.set(GRB_IntParam_OutputFlag, 0);
		env.start();

		// Create an empty model
		GRBModel model = GRBModel(env);
		model.set(GRB_IntParam_Method, 0); // Forcing Primal Simplex Method
		// Important, since fractional solution must be an Extreme Point
		
		GRBVar **vars = NULL;
		vars = new GRBVar*[n];
		for (int i = 0; i < n; i++)
			vars[i] = new GRBVar[n];
		

		// Create all variables. Maybe not needed
		for (int i = 0; i < n; i++)
			for (int j = 0; j <= i; j++){
				vars[i][j] = model.addVar(0.0, 0.0, 0.0, GRB_CONTINUOUS, "x_" + to_string(i) + "_" + to_string(j));
				vars[j][i] = vars[i][j];
			}    

		// Setting the correct UB and OBJ.	
		for (ListGraph::EdgeIt e(G); e != INVALID; ++e){
			int u = G.id(G.u(e));
			int v = G.id(G.v(e));

			vars[u][v].set(GRB_DoubleAttr_UB, 1);
			vars[u][v].set(GRB_DoubleAttr_Obj, cost[e]);
		}

		// Add \delta(v) >= 2, constraints
		for (int i = 0; i < n; i++){
			
			GRBLinExpr expr = 0;
			for (int j = 0; j < n; j++)
				expr += vars[i][j];
			
			model.addConstr(expr >= 2, "cut2_" + to_string(i));
		}
		
		// Optimize model
		model.optimize();

		bool found_feasible = 0;
		while (model.get(GRB_IntAttr_SolCount) > 0 and !found_feasible){
		
			double **sol = new double*[n];
			for (int i = 0; i < n; i++)
				sol[i] = model.get(GRB_DoubleAttr_X, vars[i], n);

			pair<double, vector<Edge> > min_cut = FindMinCut(sol);

			// If min_cut.fist < 2, need to add constraint
			if (sign(min_cut.first - 2.0) < 0) { // Min cut < 2
				GRBLinExpr expr = 0;
				
				// All all edges of the global min cut as a contraint.
				for (Edge e : min_cut.second){
					int u = G.id(G.u(e));
					int v = G.id(G.v(e));				
					
					expr += vars[u][v];
				}

				model.addConstr(expr >= 2);
				model.optimize();
			}
			else{
				// Found a feasible opt
				// cout << "Obj: " << model.get(GRB_DoubleAttr_ObjVal) << endl;

				double **sol = new double*[n];
				for (int i = 0; i < n; i++)
					sol[i] = model.get(GRB_DoubleAttr_X, vars[i], n);

				for (ListGraph::EdgeIt e(G); e != INVALID; ++e){
					int u = G.id(G.u(e));
					int v = G.id(G.v(e));

					FracSol[e] = sol[u][v];
					// cout<<u<<' '<<v<<' '<<sol[u][v]
				}

				found_feasible = 1;
			}


			for (int i = 0; i < n; i++)
				delete[] sol[i];
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
	DFS Step of BDS Algorithm. The next tree edge is chosen by the following criteria.
		- Matching edge
		- Heavy edge that maximizes x^*_e
	
	O(n^2)
*/
vector< vector<Node> > tree_adj;
void BDSDFS(ListGraph::Node v, 
	ListGraph::NodeMap<int> &parent, 
	ListGraph::EdgeMap<double> &FracSol, 
	ListGraph::EdgeMap<int> &BDSSol,
	int &clk,
	ListGraph::NodeMap<int> &in,
	ListGraph::NodeMap<int> &out){

	in[v] = clk++;
	int next_arc_id = -1;

	do {
		next_arc_id = -1;

		// Here we iterate through the arcs leaving node v.
		for (ListGraph::OutArcIt e(G, v); e != INVALID; ++e){
			ListGraph::Node u = G.target(e);

			if (parent[u] == -1){

				if ((next_arc_id == -1) or (cost[e] == 0))
					next_arc_id = G.id(e);
				
				else {
					ListGraph::Arc f = G.arcFromId(next_arc_id);
					if (cost[f] == 1 and FracSol[f] < FracSol[e])
						next_arc_id = G.id(e);
				}

			}
		}

		if (next_arc_id != -1){
			ListGraph::Arc f = G.arcFromId(next_arc_id);
			BDSSol[f] = 1;
			
			ListGraph::Node u = G.target(f);

			parent[u] = G.id(v);
			tree_adj[G.id(v)].push_back(u);

			BDSDFS(u, parent, FracSol, BDSSol, clk, in, out);
		}

	} while (next_arc_id != -1);

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
	ListGraph::NodeMap<int> &parent, 
	ListGraph::NodeMap<int> &in, 
	ListGraph::NodeMap<int> &out,
	ListGraph::EdgeMap<int> &BDSSol){
	
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

				Node h = w, lst = w, p = G.nodeFromId(parent[v]);
				while (h != p){
					for (Node y : tree_adj[G.id(h)])
						if (y != lst){
							// cout<<" Transition "<<G.id(v) + 1<<' '<<G.id(y) + 1<<endl;
							subtree_cost += UpLinkDP(y, memo, dp_edge, parent, in, out, BDSSol);
						}

					lst = h;
					h = G.nodeFromId(parent[h]);
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
	ListGraph::NodeMap<int> &parent, 
	ListGraph::NodeMap<int> &in, 
	ListGraph::NodeMap<int> &out,
	ListGraph::EdgeMap<int> &BDSSol){

	Edge e = dp_edge[v];
	BDSSol[e] = 1;

	Node u = G.u(e);
	Node w = G.v(e);
	if (Dec(w, u, in, out)) // If u is a decendent of w
		swap(w, u);

	// u is the ancestor node, w is the descendent

	// cout<<"Recovering Solution Node "<<G.id(v) + 1<<' '<<" Edge "<<G.id(u) + 1<<' '<<G.id(w) + 1<<endl;

	Node p = G.nodeFromId(parent[v]), lst = w;
	while (w != p){
		for (Node y : tree_adj[G.id(w)])
			if (y != lst)
				RecoverUpLinkSol(y, dp_edge, parent, in, out, BDSSol);

		lst = w;
		w = G.nodeFromId(parent[w]);
	}
}

void UpLinkAugmentation(
	ListGraph::EdgeMap<int> &BDSSol,
	ListGraph::NodeMap<int> &parent, 
	ListGraph::NodeMap<int> &in,
	ListGraph::NodeMap<int> &out){

	ListGraph::NodeMap<int> memo(G, -1); 
	ListGraph::NodeMap<Edge> dp_edge(G); 

	for (Node v : tree_adj[0]){
		UpLinkDP(v, memo, dp_edge, parent, in, out, BDSSol);
		RecoverUpLinkSol(v, dp_edge, parent, in, out, BDSSol);
	}
}

/*
	Implementation of the BDS MAP algorithm.
	Must receive a extreme point solution of the cut LP.

	O(n^2|L|)
*/
void BDSAlgorithm(ListGraph::EdgeMap<double> &FracSol, ListGraph::EdgeMap<int> &BDSSol){
	int n = countNodes(G);

	// Step 1, find a DFS Tree
	tree_adj.resize(n);
	for (int i = 0; i < n; i++)
		tree_adj[i].clear();

	ListGraph::NodeMap<int> parent(G, -1), in(G), out(G); 

	int cnt = 0;
	parent[G.nodeFromId(0)] = 0;
	BDSDFS(G.nodeFromId(0), parent, FracSol, BDSSol, cnt, in, out);

	// cout<<"BDS tree"<<endl;
	// for (NodeIt v(G); v != INVALID; ++v){
	// 	cout<<parent[v] + 1<<" <-- "<<G.id(v) + 1<<endl;
	// }


	// Step 2, uplink only augmentation
	UpLinkAugmentation(BDSSol, parent, in, out);

	// Sanity check, checks is BDS returned a feasible solution
	ListGraph::NodeMap<bool> ones(G, 1);
	ListGraph::EdgeMap<bool> BoolBDSSol(G, 0);
	for (ListGraph::EdgeIt e(G); e != INVALID; ++e)
		BoolBDSSol[e] = BDSSol[e];
	SubGraph<ListGraph> H(G, ones, BoolBDSSol);

	assert(biEdgeConnected(H) == 1);
}


/*
	Nauty Reader
*/


bool __found_feasible;
int __cur_graph_id;
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

	ListGraph::EdgeMap<int> BDSSol(G);
	BDSAlgorithm(FracSol, BDSSol);

	int cost_Int = 0;
	int cost_BDS = 0;
	double cost_Frac = 0;

	for (ListGraph::EdgeIt e(G); e != INVALID; ++e){
		int u = G.id(G.u(e));
		int v = G.id(G.v(e));

		cost_Int +=  IntSol[e] * cost[e];
		cost_Frac +=  FracSol[e] * cost[e];
		cost_BDS +=  BDSSol[e] * cost[e];

		// cout<<u + 1<<' '<<v + 1<<' '<<FracSol[e]<<' '<<IntSol[e]<<' '<<BDSSol[e]<<endl;
	}

	// Found a feasible example, print to file
	if (sign(3.0 * cost_Int - 4.0 * cost_Frac) >= 0 or sign(3.0 * cost_BDS - 4.0 * cost_Frac) >= 0){
		if (__found_feasible == 0){ // First feasible example found
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
		log_out << "Int/Frc = " << (double) cost_Int/cost_Frac << " BDS/Frc = " << cost_BDS/cost_Frac << endl;
		log_out << endl;
	}
}


/*
	Backtracking algorithm to find all matchings of G.
	Matched edges are marked with cost 0 on the global
	EdgeMap cost. The running time is exponential
*/
void FindAllMatchings(int e_id, int &m, int &total_matchings, ListGraph::NodeMap<bool> &matched){
	if (e_id >= m){
		SolveCurrentMatching(total_matchings);
		return;
	}

	// Case 1 : won't add edge e_id to the matching
	FindAllMatchings(e_id + 1, m, total_matchings, matched); 

	// Case 2 : if possible, will add e_id to the matching
	ListGraph::Edge e = G.edgeFromId(e_id);
	if ((matched[G.u(e)] == 0) and (matched[G.v(e)] == 0)){ // May add e_id
		
		matched[G.u(e)] = 1;
		matched[G.v(e)] = 1;
		cost[e] = 0;
		total_matchings++;

		FindAllMatchings(e_id + 1, m, total_matchings, matched);

		matched[G.u(e)] = 0;
		matched[G.v(e)] = 0;
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

	int total_matchings = 1, m = countEdges(G);
	FindAllMatchings(0, m, total_matchings, matched);

	if (__found_feasible == 1)
		g_out << "Number of matchings : " << total_matchings << endl;
}



/*
	This functions receiv nauty's geng output from stdin(may modify this),
	build a LEMON graph and the log files, and calls the function
	that iterates through all matchings.
*/
void RunNautyInput(){
	int cnt = 1;
	while (readNautyGraph(G, cin)){
		int n = countNodes(G);
		int m = countEdges(G);

		if (cnt == 1){ // Create folder to log files, create log stream
			filesystem::create_directory("./" + to_string(n));
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
	}

	log_out.close();
}


signed main(){
	ListGraph G;

	RunNautyInput();

}



// signed main(){
// 	ReadStdioInput();

// 	ListGraph::EdgeMap<int> IntSol(G);
// 	IntegerSolution(IntSol);

// 	ListGraph::EdgeMap<double> FracSol(G);
// 	FractionalSolution(FracSol);

// 	ListGraph::EdgeMap<int> BDSSol(G);
// 	BDSAlgorithm(FracSol, BDSSol);


// 	int cost_Int = 0;
// 	int cost_BDS = 0;
// 	double cost_Frac = 0;

// 	for (ListGraph::EdgeIt e(G); e != INVALID; ++e){
// 		int u = G.id(G.u(e));
// 		int v = G.id(G.v(e));

// 		cost_Int +=  IntSol[e] * cost[e];
// 		cost_Frac +=  FracSol[e] * cost[e];
// 		cost_BDS +=  BDSSol[e] * cost[e];

// 		cout<<u + 1<<' '<<v + 1<<' '<<FracSol[e]<<' '<<IntSol[e]<<' '<<BDSSol[e]<<endl;
// 	}

// 	cout<<"Cost Fractional "<<cost_Frac<<endl;
// 	cout<<"Cost Integral "<<cost_Int<<endl;
// 	cout<<"Cost BDS "<<cost_BDS<<endl;
// }
