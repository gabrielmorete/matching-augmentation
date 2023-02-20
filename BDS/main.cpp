#include <iostream>
#include <vector>
#include <array>
#include <cassert>
#include <string>
#include "/Library/gurobi1001/macos_universal2/include/gurobi_c++.h"
#include <lemon/list_graph.h>
#include <lemon/gomory_hu.h>

using namespace std;
using namespace lemon;  

typedef long long int ll;

const double EPS = 1e-8;
int sign(double x) { return (x > EPS) - (x < -EPS); }


/*
	Conventions
		- On lemon, graphs are 0 indexed;
		- Graph is simple
*/


ListGraph G; // Declare global Graph
ListGraph::EdgeMap<int> cost(G); // Cost of the edges


// I will also save the input information
int _n, _m;
vector< array<int, 3> > _edges;


/*
	This function reads the Graph. The input format will be
	n m       number of nodes, edges
	a b c     edge bertween a, b with cost c
*/
void ReadInput(){
	int n, m;
	cin>>n>>m;
	
	_n = n;
	_m = m;

	for (int i = 0; i < n; i++){
		ListGraph::Node v = G.addNode();
		if (G.id(v) != i)
			cout<<"Error : node don't match id"<<endl;
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
	This is the separator function. If there is a cut with value < 2
	we will add this cut to the constraints.
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
				if (where == GRB_CB_MIPSOL){ // Found a solution
					ListGraph::EdgeMap<double> capacities(G);

					double *x[n];
					for (int i = 0; i < n; i++)
						x[i] = getSolution(vars[i], n);

					ListGraph::EdgeMap<double> capacity(G);
					
					for (ListGraph::EdgeIt e(G); e != INVALID; ++e){
						int u = G.id(G.u(e));
						int v = G.id(G.v(e));

						capacity[e] = x[v][u];
					}


					GomoryHu<ListGraph, ListGraph::EdgeMap<double> > GMH(G,capacity);
					GMH.run();	                

					double min_val = 3;
					ListGraph::Node min_v = G.nodeFromId(0);
					ListGraph::Node min_u = G.nodeFromId(1);

					for (ListGraph::NodeIt v(G); v != INVALID; ++v)
						for (ListGraph::NodeIt u(G); u != INVALID; ++u){
							if (GMH.minCutValue(v, u) < GMH.minCutValue(min_v, min_u)){
								min_v = v;
								min_u = u;
							}
					}

					if (sign(GMH.minCutValue(min_v, min_u)) < 0) { // Min cut < 2
						GRBLinExpr expr = 0;
						
						for (GomoryHu<ListGraph, ListGraph::EdgeMap<double> >::MinCutEdgeIt e(GMH, min_v, min_u); e != INVALID; ++e){
							int u = G.id(G.u(e));
							int v = G.id(G.v(e));				
							
							expr += vars[u][v];
						}

						addLazy(expr >= 2);
					}
			}		

			} catch (GRBException e) {
				cout << "Error number: " << e.getErrorCode() << endl;
				cout << e.getMessage() << endl;
			} catch (...) {
				cout << "Error during callback" << endl;
			}
		
	}
};






// This function receivs a graph and returns a optmal fractional solution
void IntegerSolution(){
	try {
		
		int n = countNodes(G);

		// Create an environment
		GRBEnv env = GRBEnv(true);
		env.set("LogFile", "Integer.log");
		env.start();

		// Create an empty model
		GRBModel model = GRBModel(env);
		model.set(GRB_IntParam_LazyConstraints, 1);
		
		GRBVar **vars = NULL;
		vars = new GRBVar*[n];
		for (int i = 0; i < n; i++)
			vars[i] = new GRBVar[n];
		

		for (int i = 0; i < n; i++)
			for (int j = 0; j < i; j++){
				// vars[i][j] = model.addVar(0.0, 0.0, 0.0, GRB_BINARY, "x_" + to_string(i) + "_" + to_string(j));
				// vars[j][i] = vars[i][j];
			}    

		for (ListGraph::EdgeIt e(G); e != INVALID; ++e){
			int u = G.id(G.u(e));
			int v = G.id(G.v(e));

			// vars[u][v].set(GRB_DoubleAttr_UB, 1);
			// vars[u][v].set(GRB_DoubleAttr_Obj, cost[e]);
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

		// cout << x.get(GRB_StringAttr_VarName) << " "
		//  << x.get(GRB_DoubleAttr_X) << endl;
		// cout << y.get(GRB_StringAttr_VarName) << " "
		//  << y.get(GRB_DoubleAttr_X) << endl;
		// cout << z.get(GRB_StringAttr_VarName) << " "
		//  << z.get(GRB_DoubleAttr_X) << endl;

		cout << "Obj: " << model.get(GRB_DoubleAttr_ObjVal) << endl;

	} catch(GRBException e) {
		cout << "Error code = " << e.getErrorCode() << endl;
		cout << e.getMessage() << endl;
	} catch(...) {
		cout << "Exception during optimization" << endl;
	}

}




signed main(){
	ReadInput();
}
