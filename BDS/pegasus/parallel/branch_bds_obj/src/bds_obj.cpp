#include "bds_obj.h"

/*
	I need lemon to find global min cut in a eficient and stable way.
	For bds I can implement it more efficiently without lemon.
	You may find a lemon based implementation in bds_lemon.cpp
*/

/*
	Implementation of the BDS MAP algorithm.
	Must receive a extreme point solution of the cut LP.

	O(n|L|)

	Declare the structure once per graph. 
	Then just updateit for different costs and lp solutions.
*/

/*
	DFS Step of BDS Algorithm. The next tree edge is chosen by the following criteria.
		- Matching edge
		- Heavy edge that maximizes x^*_e
	
	O(n + m)
*/
void BDSAlgorithm::Dfs(int v){
	in[v] = clk++;

	int matched = 1 - cost[adj[v][0]];

	for (int i = 1 + matched; i < adj[v].size(); i++)
		assert(sign(lp[adj[v][i - 1]] - lp[adj[v][i]]) >= 0);


	for (int e : adj[v]){
		// Edges are sorted and the algorithm runs in the support
		if (sign(lp[e]) <= 0)
			break;

		int u = e_u[e] + e_v[e] - v;
		if (in[u] == 0){
			parent[u] = v;
			in_sol[e] = 1;
			tree_adj[v].push_back(u);

			Dfs(u);
		}
	}

	out[v] = clk++;
}

/*
	Retuns true if v is a proper descendent of u.
*/
inline bool BDSAlgorithm::StrictDec(int u, int v){ // is v strict dec of u
	return (in[u] < in[v]) and (out[v] < out[u]);
}

/*
	Retuns true if v is a descendent of u.
*/
inline bool BDSAlgorithm::Dec(int u, int v){ // is v strict dec of u
	return (in[u] <= in[v]) and (out[v] <= out[u]);
}


/*
	Dynamic programming on a tree, follows reverse topological sort.

	memo_edge is ists original cost + the cost of covering every danglind subtree.
	At time it is coveding {v, parent[v]}, the cost of the esdge ist he same as
	memo[v]

	O(n|L|)
*/
void BDSAlgorithm::UpLinkDP(int v){
	for (auto x : tree_adj[v])
		UpLinkDP(x);

	dp_edge[v] = cover[v][0];
	for (int e : cover[v])
		if (memo_edge[e] < memo_edge[dp_edge[v]])
			dp_edge[v] = e;

	if (parent[v] != 0)
		for (int e : cover[parent[v]])
			if (Dec(v, e_u[e]) == 0 and Dec(v, e_v[e]) == 0)
				memo_edge[e] += memo_edge[ dp_edge[v] ];
}

/*
	Algorithm to recover the uplink DP solution.

	O(n + m)
*/
void BDSAlgorithm::RecoverUpLinkSol(int v){
	int eid = dp_edge[v];	

	in_sol[eid] = 1;

	int u = e_u[eid];
	int w = e_v[eid];

	if (Dec(w, u)) // w is the lower vertex
		swap(u, w);

	int lst = w;
	while (w != parent[v]){
		for (int h : tree_adj[w])
			if (h != lst)
				RecoverUpLinkSol(h);

		lst = w;
		w = parent[w];
	}
}


void BDSAlgorithm::UpLinkAugmentation(){
	for (int v = 0; v < n; v++)
		cover[v].clear();

	for (int i = 0; i < m; i++)
		if (!in_sol[i] and sign(lp[i]) > 0){
			int u = e_u[i];
			int v = e_v[i];
			if (Dec(u, v)) // u is the lower vertex
				swap(u, v);

			while (u != v){ // link i covers {v, parent[v]}
				cover[u].push_back(i);

				if (parent[u] == u){
					dbg(u);
					dbg(lp[i]);

				}

				assert(u != parent[u]);
				u = parent[u];
			}

			memo_edge[i] = cost[i];
		}

	// Preprocess
	for (auto u : tree_adj[0]){
		UpLinkDP(u);
		// RecoverUpLinkSol(u);
	}
}

BDSAlgorithm::BDSAlgorithm(){
	
}


/*
	Builds structure for new graph
*/
BDSAlgorithm::BDSAlgorithm(ListGraph &G){
	n = countNodes(G);
	m = countEdges(G);


	adj.resize(n);
	tree_adj.resize(n);

	for (int i = 0; i < n; i++)
		adj[i].clear();

	in.resize(n);
	parent.resize(n);
	out.resize(n);
	dp_edge.resize(n);
	cover.resize(n);

	e_u.resize(m);
	e_v.resize(m);
	cost.resize(m);
	in_sol.resize(m);
	lp.resize(m);
	memo_edge.resize(m);


	for (ListGraph::EdgeIt e(G); e != INVALID; ++e){
		int eid = G.id(e);
		e_u[eid] = G.id(G.u(e));
		e_v[eid] = G.id(G.v(e));
		adj[e_v[eid]].push_back(eid);
		adj[e_u[eid]].push_back(eid);
	}
}

/*
	Update costs and lp value
*/
void BDSAlgorithm::Update(ListGraph::EdgeMap<int> &_cost, ListGraph::EdgeMap<double> &_FracSol, ListGraph &G){
	for (ListGraph::EdgeIt e(G); e != INVALID; ++e){
		int eid = G.id(e);
		lp[eid] = _FracSol[e];
		cost[eid] = _cost[e];
		in_sol[eid] = 0; 
	}

	// Matched edge is the first
	for (int v = 0; v < n; v++){
		int p = 0;

		for (int e = 0; e < adj[v].size(); e++)
			if (cost[adj[v][e]] == 0 and (sign(lp[adj[v][e]]) > 0)) // If matched edge is zero, skip
				p = e;
		
		swap(adj[v][0], adj[v][p]);	
	
		int matched = 1 - cost[adj[v][0]];

		sort(adj[v].begin() + matched, adj[v].end(),
			[this](int a, int b){ // Sort by increase lp value, skip matching edge
				return lp[a] > lp[b];
			}
		);
	}

	// for (int v = 0; v < n; v++){

	// 	for (int e : adj[v])
	// 		cout<<"(" << cost[e] << ", " << lp[e] << ") ";
	// 	cout<<endl;

	// 	int matched = 1 - cost[adj[v][0]];

	// 	for (int i = 1 + matched; i < adj[v].size(); i++)
	// 		assert(sign(lp[adj[v][i - 1]] - lp[adj[v][i]]) >= 0);
	// }
}

#warning "remove returns"

/*
	Run BDS algorithm
*/
void BDSAlgorithm::Run(ListGraph::EdgeMap<int> &_cost, 
	ListGraph::EdgeMap<bool> &BDSSol, 
	ListGraph::EdgeMap<double> &FracSol, 
	ListGraph &G){
	
	Update(_cost, FracSol, G);	

	// Step 1, find a DFS Tree
	for (int v = 0; v < n; v++){
		tree_adj[v].clear();
		parent[v] = v;
		in[v] = 0;
	}

	clk = 1;
	Dfs(0);

	if (__verbose_mode){
		cout << "BDS Tree Found" << endl;
		for (int v = 0; v < n; v++)
			cout << parent[v] << " <-- " << v << endl;	
	}

	// Step 2, uplink only augmentation
	UpLinkAugmentation();

	return; /////////////////////////////////////////////////////////////////////////////////////////////////////////////

	for (ListGraph::EdgeIt e(G); e != INVALID; ++e)
		BDSSol[e] = in_sol[G.id(e)];

	ListGraph::NodeMap<bool> ones(G, 1);

	// Sanity check, checks if BDS returned a feasible solution
	SubGraph<ListGraph> H(G, ones, BDSSol);
	assert(biEdgeConnected(H) == 1);
}


