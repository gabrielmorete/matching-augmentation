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

void print(ListGraph &G, int msk){

	for (ListGraph::EdgeIt e(G); e != INVALID; ++e)
		if (msk & (1 << (G.id(e))))
			cout << G.id(G.v(e)) << ' ' << G.id(G.u(e)) << ", ";
	
	cout << endl;
}


bool test(SubGraph<ListGraph> &G){ // for now, only simple cycle
	if (biEdgeConnected(G) != 1)
		return false;

	map<pair<int, int>, int> frq;

	SubGraph<ListGraph>::NodeMap<int> deg(G, 0);


	for (SubGraph<ListGraph>::EdgeIt e(G); e != INVALID; ++e){
		int u = min(G.id(G.u(e)), G.id(G.v(e)));
		int v = max(G.id(G.u(e)), G.id(G.v(e)));
		
		deg[G.v(e)]++;
		deg[G.u(e)]++;

		frq[{u, v}]++;
	}

	for (SubGraph<ListGraph>::EdgeIt e(G); e != INVALID; ++e){
		int u = min(G.id(G.u(e)), G.id(G.v(e)));
		int v = max(G.id(G.u(e)), G.id(G.v(e)));
		
		if (frq[{u, v}] > 1) // doubled edge
			if (deg[G.u(e)] > 2 and deg[G.v(e)] > 2) // no head
				return false;
	}

	return true;
}

bool test2(SubGraph<ListGraph> &G){ // for now, only simple cycle
	if (biEdgeConnected(G) != 1)
		return false;

	map<pair<int, int>, int> frq;

	SubGraph<ListGraph>::NodeMap<int> deg(G, 0);
	SubGraph<ListGraph>::NodeMap<int> mark(G, 0);



	for (SubGraph<ListGraph>::EdgeIt e(G); e != INVALID; ++e){
		int u = min(G.id(G.u(e)), G.id(G.v(e)));
		int v = max(G.id(G.u(e)), G.id(G.v(e)));
		
		deg[G.v(e)]++;
		deg[G.u(e)]++;

		if (frq[{u, v}] == 1){
			mark[G.v(e)] += 1;
			mark[G.u(e)] += 1;
		}

		frq[{u, v}]++;
	}

	/*
		mark[v] == 0 -> not incidente to multiedge
		mark[v] == 1 -> head of a chain
		mark[v] == 1 -> inner vertex of a chain

	*/

	for (SubGraph<ListGraph>::NodeIt v(G); v != INVALID; ++v){
		if (mark[v] == 1){ // head of the chain
			if (deg[v] == 2) // leaf node
				continue;
			else{
				// need to find the other head to check if it is a leaf
				SubGraph<ListGraph>::Node lst = v, u = v;

				do {
					for (SubGraph<ListGraph>::OutArcIt a(G, v); a != INVALID; ++a)
						if (G.target(a) != lst){
							lst = u;
							u = G.target(a);
							break;
						}
				} while(mark[u] == 2);

				if (deg[u] > 2) // also not a leaf
					return false;
			}
		}

	}

	return true;
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


// Gurobi enviroment
GRBEnv env = GRBEnv(true);


/*
	Model to find the smallest coefficient to each point.

	Note: Rebuilding the model everytime is not efficient, but 
	for the current application is ok.
*/
double ConvexComb(double *sol, int dim, int G, vector<int> H, int op = 0){ //op = 0 (<=), op = 1 (==)
	
	int n = H.size(); // number of points
	
	try{
		GRBModel model(env);
		GRBVar lambda[n];

		model.set(GRB_IntParam_Method, 0); // Forcing Primal Simplex Method

		// one variable for each subgraph (combination coefficient)
		for (int i = 0; i < n; i++)
			lambda[i] = model.addVar(0.0, 1.0, 0, GRB_CONTINUOUS, "l_" + to_string(i) );

		GRBLinExpr conv;
		for (int i = 0; i < n; i++)
			conv += lambda[i];

		model.addConstr(conv == 1, "conv_comb");

		// Model will thy to optmize the coefficient of the combination
		GRBVar coef = model.addVar(0.0, GRB_INFINITY, 1, GRB_CONTINUOUS, "coef");

		// combination constraint
		for (int j = 0; j < dim; j++){ // Each edge
			GRBLinExpr comb = 0;
			
			for (int i = 0; i < n; i++){ // for each subgraph
				int x = 0;
				if (H[i] & (1<<j)) // edge is in the graph
					x = 1;

				comb += x * lambda[i];
			}

			int y = 0;
			if (G & (1 << j)) // edge is in G
				y = 1;

			if (op)
				model.addConstr( comb == coef * y, "coord_" + to_string(j)); 
			else	
				model.addConstr( comb <= coef * y, "coord_" + to_string(j)); 
		}

		model.optimize();
		assert(model.get(GRB_IntAttr_SolCount) > 0);

		double *opt_sol = model.get(GRB_DoubleAttr_X, lambda, n);

		for (int i = 0; i < n; i++)
			sol[i] = opt_sol[i];

		delete[] opt_sol;

		return model.get(GRB_DoubleAttr_ObjVal);
	
	} catch(GRBException e) {
		cout << "Error code = " << e.getErrorCode() << endl;
		cout << e.getMessage() << endl;
	} catch(...) {
		cout << "Exception during optimization" << endl;
	}

	return -1;
}

double ConvexComb2(double *sol, int dim, int G, vector<int> H, int op = 0){ //op = 0 (<=), op = 1 (==)
	
	int n = H.size(); // number of points
	
	double coef = 2.0/3.0;

	try{
		GRBModel model(env);
		GRBVar lambda[n];
		GRBVar y[n]; // used to minimize number of positive variables


		// model.set(GRB_IntParam_Method, 0); // Forcing Primal Simplex Method

		// one variable for each subgraph (combination coefficient)
		for (int i = 0; i < n; i++)
			lambda[i] = model.addVar(0.0, 1.0, 0, GRB_CONTINUOUS, "l_" + to_string(i) );

		for (int i = 0; i < n; i++){
			y[i] = model.addVar(0.0, 1.0, 1, GRB_BINARY, "l_" + to_string(i) );
			model.addConstr(lambda[i] <= y[i]);
		}

		GRBLinExpr conv;
		for (int i = 0; i < n; i++)
			conv += lambda[i];

		model.addConstr(conv == 1, "conv_comb");

		
		// combination constraint
		for (int j = 0; j < dim; j++){ // Each edge
			GRBLinExpr comb = 0;
			
			for (int i = 0; i < n; i++){ // for each subgraph
				int x = 0;
				if (H[i] & (1<<j)) // edge is in the graph
					x = 1;

				comb += x * lambda[i];
			}

			int y = 0;
			if (G & (1 << j)) // edge is in G
				y = 1;

			if (op)
				model.addConstr( comb == coef * y, "coord_" + to_string(j)); 
			else	
				model.addConstr( comb <= coef * y, "coord_" + to_string(j)); 
		}

		model.optimize();
		assert(model.get(GRB_IntAttr_SolCount) > 0);

		double *opt_sol = model.get(GRB_DoubleAttr_X, lambda, n);

		for (int i = 0; i < n; i++)
			sol[i] = opt_sol[i];

		delete[] opt_sol;

		return model.get(GRB_DoubleAttr_ObjVal);
	
	} catch(GRBException e) {
		cout << "Error code = " << e.getErrorCode() << endl;
		cout << e.getMessage() << endl;
	} catch(...) {
		cout << "Exception during optimization" << endl;
	}

	return -1;
}


signed main(){
	env.set(GRB_IntParam_OutputFlag, 0);
	env.start();


	ListGraph G;

	ReadStdioGraph(G);

	print(G);

	int n = countNodes(G);
	int m = countEdges(G);
	ListGraph::NodeMap<bool> ones(G, 1);	

	vector<int> base;

	assert(m < 30);
	for (int msk = 1; msk < (1 << m); msk++){
		if (__builtin_popcount(msk) < n)
			continue;

		ListGraph::EdgeMap<bool> mask(G, 0);
		for (int i = 0; i < m; i++)
			if (msk & (1<<i))
				mask[G.edgeFromId(i)] = 1;

		SubGraph<ListGraph> H(G, ones, mask);
	
		if (test(H)){
			base.push_back(msk);
		}
	}

	// Now I have every valid 2ECSS
	// Test convex comb.
	double sol[base.size()];

	int fmsk =  (1<<m) - 1;
	cout << ConvexComb(sol, m, fmsk, base) << endl;

	for (int i = 0; i < m; i++){
		int u = min(G.id(G.u(G.edgeFromId(i))), G.id(G.v(G.edgeFromId(i))));
		int v = max(G.id(G.u(G.edgeFromId(i))), G.id(G.v(G.edgeFromId(i))));

		cout << '\t' << u << ' ' << v << ' ' << ConvexComb2(sol, m, fmsk - (1<<i), base, 1) << endl;
	
		for (int i = 0; i < base.size(); i++)
			if (sol[i] > 0.01){
				cout << "\t\t" << sol[i] << ' '; 
				print(G, base[i]);
			}
	}
}