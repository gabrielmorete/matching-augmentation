/*
	Program receives two lists A, B of points and tests if every 
	element of the list A can be expressed as a convex combination of the 
	elements of the list B.

	Let a \in A, and q be a rational number

		q . x >= \sum_{b in B} l_b b
		\sum_{b in B} l_b = 1
		l_b >= 0, b \in B

	To run the program, flags -verbose -coef are not mandatory.
	Use -verbose for extra information
	Use -coef to overwrite the default values of the coefficients

	Author : Gabriel Morete	
*/

#include <iostream>
#include <sstream>
#include <fstream>
#include <vector>
#include <bitset>
#include <string>
#include <cassert>
#include <iomanip>
#include <lemon/list_graph.h>
#include <lemon/gomory_hu.h>
#include <lemon/matching.h>
#include <lemon/adaptors.h>
#include <lemon/connectivity.h>
#include <lemon/nauty_reader.h>
#include "gurobi_c++.h"


using namespace std;
using namespace lemon;

void print(ListGraph &G){
	for (ListGraph::EdgeIt e(G); e != INVALID; ++e)
		cout << G.id(G.v(e)) << ' ' << G.id(G.u(e)) << ", ";
	cout << endl;
}

void print(SubGraph<ListGraph> &G){
	for (SubGraph<ListGraph>::EdgeIt e(G); e != INVALID; ++e)
		cout << G.id(G.v(e)) << ' ' << G.id(G.u(e)) << ", ";
	cout << endl;
}

void print(ListGraph &G, long long int msk){

	for (ListGraph::EdgeIt e(G); e != INVALID; ++e)
		if (msk & (1 << (G.id(e))))
			cout << G.id(G.v(e)) << ' ' << G.id(G.u(e)) << ", ";
	
	cout << endl;
}

/*
	This function reads the Graph from Stdio. Graph is 0-indexed
	The input format will be
		n m       number of nodes, edges
		a_1 b_1    edge between a_1, b_1
		...
		a_m b_m
*/
void ReadStdioGraph(ListGraph &G){
	int n, m;
	cin>>n>>m;
	
	assert(n >= 3);
	assert(m >= n);

	for (int i = 0; i < n; i++){
		ListGraph::Node v = G.addNode();
		if (G.id(v) != i)
			cout<<"Error : vertex don't match id"<<endl;
		assert(G.id(v) == i);
	}

	for (int i = 0; i < m; i++){
		int a, b;
		cin>>a>>b;

		assert(a < n);
		assert(b < n);

		G.addEdge(G.nodeFromId(a), G.nodeFromId(b));
	}

	set<int> q;
	for (ListGraph::EdgeIt e(G); e != INVALID; ++e)
		q.insert(G.id(e));

	cout << q.size() << endl;

	assert(q.size() == m);
	assert(*q.rbegin() == m - 1);
}

// Checks if the graph is 4-reg, 4-ec
bool check(ListGraph &G){
	for (ListGraph::NodeIt v(G); v != INVALID; ++v)
		if (countIncEdges(G, v) != 4){
			// cout << "Not 4-regular" << endl;
			return 0;
		}

	int m = countEdges(G);
	// Check 4-edge connected by definition

	ListGraph::NodeMap<bool> ones(G, 1);
	ListGraph::EdgeMap<bool> edges(G, 1);

	for (int i = 0; i < m; i++)
		for (int j = i + 1; j < m; j++)
			for (int k = j + 1; k < m; k++){
				edges[G.edgeFromId(i)] = 0;
				edges[G.edgeFromId(j)] = 0;
				edges[G.edgeFromId(k)] = 0;

				SubGraph<ListGraph> H(G, ones, edges);

				if (connected(H) == 0){
					return 0;
				}

				edges[G.edgeFromId(i)] = 1;
				edges[G.edgeFromId(j)] = 1;
				edges[G.edgeFromId(k)] = 1;
			}
	return true;
}


// Gurobi enviroment
GRBEnv env = GRBEnv(true);


/*
	Callback function class.
*/
class MinimumCut: public GRBCallback {
	public:
		ListGraph *G;
		GRBVar** vars;
		int n, m;
		
		// Constructor for min cut
		MinimumCut(GRBVar** _vars, int _n, int _m, ListGraph &_G){
			vars = _vars;
			n = _n;
			m = _m;
			G = &_G;
		}
	protected:
		/*
			Separator function for the MIP. 
			If the solution is not a 2ECSS it adds a cut
			separating one 2ECC.
		*/
		void callback(){
			try {
				if (where == GRB_CB_MIPSOL){
					/*
						Solver found an integral optimal solution for the
						current formulation, must check if if there is a
						bridge or a cut.
					*/

					double *x[3];

					for (int i = 0; i < 3; i++)
						x[i] = getSolution(vars[i], m);

					for (int g = 0; g < 3; g++){

						ListGraph::EdgeMap<bool> in_sol(*G);
						
						for (ListGraph::EdgeIt e(*G); e != INVALID; ++e){
							int id = (*G).id(e);
							int u = (*G).id((*G).u(e));
							int v = (*G).id((*G).v(e));

							if (x[g][id] > 0.5)
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

							// Will add all cuts crossing 2ECC
							GRBLinExpr expr[max_cmp + 1];

							for (ListGraph::EdgeIt e((*G)); e != INVALID; ++e){
								ListGraph::Node u = (*G).u(e);
								ListGraph::Node v = (*G).v(e);
								
								if (ebcc[u] != ebcc[v]){
									expr[ ebcc[u] ] += vars[g][(*G).id(e)];
									expr[ ebcc[v] ] += vars[g][(*G).id(e)];
								}
							}

							for (int i = 0; i < max_cmp; i++)
								addLazy(expr[i] >= 2);
						}

					}	

					for (int i = 0; i < 3; i++)	
						delete[] x[i];
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
	coefficient is fixed to be 2/3.  model finds the combination with the minimum number of elements
*/
vector< vector<pair<int, int>> > ConvexComb(int e, ListGraph &G, int op = 0){ //op = 0 (<=), op = 1 (==)

	try{
		GRBModel model(env);
		model.set(GRB_IntParam_LazyConstraints, 1); // Allow callback constraints

		int n = countNodes(G);
		int m = countEdges(G);

		GRBVar **x = new GRBVar*[3]; // used to minimize number of positive variables

		for (int i = 0; i < 3; i++)
			x[i] = new GRBVar[m];

		MinimumCut cb = MinimumCut(x, n, m, G);
		model.setCallback(&cb);

		for (int i = 0; i < 3; i++)
			for (int j = 0; j < m; j++)
				x[i][j] = model.addVar(0.0, 1.0, 0.0, GRB_BINARY, "X_" + to_string(i) + "_" + to_string(j));


		for (int j = 0; j < m; j++){	
			GRBLinExpr conv = 0;
			for (int i = 0; i < 3; i++)
				conv += x[i][j];

			if (j != e)
				model.addConstr(conv <= 2, "e_" + to_string(j) + "<= 2"); // each edge appers at most twice
			else
				model.addConstr(conv <= 0, "e_" + to_string(j) + "<= 2"); // removed edge
		}		
		
		model.optimize();

		if (model.get(GRB_IntAttr_SolCount) == 0){
			cout << "Conjecture is false" << endl;
			exit(0);
		}

		double *opt_sol[3];

		for (int i = 0; i < 3; i++)
			opt_sol[i] = model.get(GRB_DoubleAttr_X, x[i], m);

		vector< vector<pair<int, int>> > ans(3);

		for (int i = 0; i < 3; i++){
			long long int x = 0;
			for (int j = 0; j < m; j++)
				if (opt_sol[i][j] > 0.5)
					ans[i].push_back({G.id(G.u( G.edgeFromId(j) )), G.id(G.v(  G.edgeFromId(j) ))});
		} 

		for (int i = 0; i < 3; i++)
			delete[] opt_sol[i];

		for (int i = 0; i < 3; i++)
			delete[] x[i];
		delete [] x;

		return ans;
	
	} catch(GRBException e) {
		cout << "Error code = " << e.getErrorCode() << endl;
		cout << e.getMessage() << endl;
	} catch(...) {
		cout << "Exception during optimization" << endl;
	}

	return {};
}


signed main(int argc, char *argv[]){
	env.set(GRB_IntParam_OutputFlag, 0);
	env.start();

	ListGraph G;


	if (argc == 1){
		ReadStdioGraph(G);
		if (check(G) == 0){
			cout << "Invalid input" << endl;
			return 0;
		}


		print(G);

		int n = countNodes(G);
		int m = countEdges(G);

		for (int i = 0; i < m; i++){
			int u = min(G.id(G.u(G.edgeFromId(i))), G.id(G.v(G.edgeFromId(i))));
			int v = max(G.id(G.u(G.edgeFromId(i))), G.id(G.v(G.edgeFromId(i))));

			cout << '\t' << u << ' ' << v << ' ' << endl;
			
			auto comb = ConvexComb(i, G);

			for (auto x : comb){
				cout << "\t\t"; 
				for (auto y : x)
					cout << y.first << ' ' << y.second << ", ";
				cout << endl;
			}	
		}
	}
	else{
		while (readNautyGraph(G, cin)){
			if (!check(G))
				continue;

			int m = countEdges(G);
			for (int i = 0; i < m; i++)
				auto comb = ConvexComb(i, G);
		}
	}
}
