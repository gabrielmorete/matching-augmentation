/*
	Matching Augmentation problem pack
		- Linear Relaxation
		- Integer Solution
		- BDS Approximation Algorithm 

	
	Author : Gabriel Morete
*/


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

// Safe handling doubles
const double EPS = 1e-8;
int sign(double x) { return (x > EPS) - (x < -EPS); }


/*
	Conventions
		- Input is 1-indexed
		- On lemon, graphs are 0-indexed;
		- Graph is simple, for now
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
	
	assert(n >= 3);
	assert(m >= n);

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
		env.start();

		// Create an empty model
		GRBModel model = GRBModel(env);
		model.set(GRB_IntParam_LazyConstraints, 1);
		
		GRBVar **vars = NULL;
		vars = new GRBVar*[n];
		for (int i = 0; i < n; i++)
			vars[i] = new GRBVar[n];
		

		// Create all variables. Maybe not needed
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

		if (model.get(GRB_IntAttr_SolCount) > 0){
			cout << "Obj: " << model.get(GRB_DoubleAttr_ObjVal) << endl;

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
		env.start();

		// Create an empty model
		GRBModel model = GRBModel(env);
		model.set(GRB_IntParam_LazyConstraints, 1);
		
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

		// Set callback function
    	MinimumCut cb = MinimumCut(vars, n);
    	model.setCallback(&cb);
		
		// Optimize model
		model.optimize();

		if (model.get(GRB_IntAttr_SolCount) > 0){
			cout << "Obj: " << model.get(GRB_DoubleAttr_ObjVal) << endl;

			double **sol = new double*[n];
			for (int i = 0; i < n; i++)
				sol[i] = model.get(GRB_DoubleAttr_X, vars[i], n);

			for (ListGraph::EdgeIt e(G); e != INVALID; ++e){
				int u = G.id(G.u(e));
				int v = G.id(G.v(e));

				FracSol[e] = sol[u][v];

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



signed main(){
	ReadInput();

	ListGraph::EdgeMap<int> IntSol(G);
	IntegerSolution(IntSol);

	ListGraph::EdgeMap<double> FracSol(G);
	FractionalSolution(FracSol);

	for (ListGraph::EdgeIt e(G); e != INVALID; ++e){
		int u = G.id(G.u(e));
		int v = G.id(G.v(e));

		cout<<u + 1<<' '<<v + 1<<' '<<IntSol[e]<<' '<<FracSol[e]<<endl;
	}

}
