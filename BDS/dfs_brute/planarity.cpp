// Check if a graph is planar
// Returns the subdivision of K33 or K5
// Run with any commando line argument to
// display human-readable information

#include <iostream>
#include <algorithm>
#include <vector>
#include <map>
#include <unordered_map>
#include <set>
#include <tuple>
#include <cassert>
#include <fstream>
#include <lemon/list_graph.h>
#include <lemon/planarity.h>

using namespace std;
using namespace lemon;

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
		string s1, s2;
		cin >> a >> b;

		assert(a < n);
		assert(b < n);

		G.addEdge(G.nodeFromId(a), G.nodeFromId(b));
	}

	set<int> q;
	for (ListGraph::EdgeIt e(G); e != INVALID; ++e)
		q.insert(G.id(e));

	assert(q.size() == m);
	assert(*q.rbegin() == m - 1);
}

signed main(int argc, char *argv[]){
	ListGraph G;
	ReadStdioGraph(G);

	PlanarEmbedding<ListGraph> Plan(G);

	if (!Plan.run()){
		if (argc > 1)
			cout << "Non-planar" << endl;	
		
		
		ListGraph::NodeMap<int> deg(G);
		for (ListGraph::EdgeIt e(G); e != INVALID; ++e)
			if (Plan.kuratowski(e)){
				if (argc > 1)
					cout << G.id(G.v(e)) << ' ' << G.id(G.u(e)) << endl;
				deg[G.v(e)]++;
				deg[G.u(e)]++;
			}

		int mx = 0;
		for (ListGraph::NodeIt v(G); v != INVALID; ++v)
			mx = max(mx, deg[v]);	

		if (mx <= 3){
			if (argc > 1)
				cout << "K_{3,3}" << endl;
			return 1;
		}
		else{
			if (argc > 1)
				cout << "K_5" << endl;
			return 2;
		}
	}
	else{
		if (argc > 1)
			cout << "Planar" << endl;
		return 0;
	}
}
