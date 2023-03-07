#include "libraries_and_utils.h"
#include "lemon.h"
#include "main.h"

/*
	Stdio reader
*/


/*
	This function reads the Graph from Stdio. Graph is 0-indexed
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

	for (int i = 0; i < n; i++){
		ListGraph::Node v = G.addNode();
		if (G.id(v) != i)
			cout<<"Error : vertex don't match id"<<endl;
		assert(G.id(v) == i);
	}

	for (int i = 0; i < m; i++){
		int a, b, c;
		cin>>a>>b>>c;

		// a--; // 0-indexed
		// b--; 

		ListGraph::Edge e = G.addEdge(G.nodeFromId(a), G.nodeFromId(b));
		cost[e] = c;
	}

	set<int> q;
	for (ListGraph::EdgeIt e(G); e != INVALID; ++e)
		q.insert(G.id(e));

	assert(q.size() == m);
	assert(*q.rbegin() == m - 1);
}

void RunStdioInput(){
	ReadStdioInput();

	ListGraph::EdgeMap<int> IntSol(G);
	ListGraph::EdgeMap<double> FracSol(G);
	ListGraph::EdgeMap<bool> BDSSol(G);

	SolveMapInstance(FracSol,IntSol,BDSSol);

	int cost_Int = 0;
	int cost_BDS = 0;
	double cost_Frac = 0;

	for (ListGraph::EdgeIt e(G); e != INVALID; ++e){
		int u = G.id(G.u(e));
		int v = G.id(G.v(e));

		cost_Int +=  IntSol[e] * cost[e];
		cost_Frac +=  FracSol[e] * cost[e];
		cost_BDS +=  BDSSol[e] * cost[e];

		cout<<u + 1<<' '<<v + 1<<' '<<FracSol[e]<<' '<<IntSol[e]<<' '<<BDSSol[e]<<endl;
	}

	cout<<"Cost Fractional "<<cost_Frac<<endl;
	cout<<"Cost Integral "<<cost_Int<<endl;
	cout<<"Cost BDS "<<cost_BDS<<endl;
}
