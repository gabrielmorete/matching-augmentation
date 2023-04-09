#include "bds_obj.h"

/*
	I need lemon to find global min cut in a eficient and stable way.
	For bds I can implement it more efficiently without lemon.
	You may find a lemon based implementation in bds_lemon.cpp
*/

/*
	Implementation of the BDS MAP algorithm.
	Must receive a extreme point solution of the cut LP.

	O(n + |L|log(|L|))

	Declare the structure once per graph. 
	Then just updateit for different costs and lp solutions.
*/


BDSAlgorithm::BDSAlgorithm(){}

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
	covered.resize(n);
	link.resize(n);

	e_u.resize(m);
	e_v.resize(m);
	cost.resize(m);
	in_sol.resize(m);
	lp.resize(m);

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
	
		// Be careful when matched edge is the first but LP value is zero
		int matched = 0;
		if ((cost[ adj[v][0] ] == 0) and (sign(lp[ adj[v][0] ]) > 0) )
			matched = 1;

		sort(adj[v].begin() + matched, adj[v].end(),
			[this](int a, int b){ // Sort by increase lp value, skip matching edge
				return lp[a] > lp[b];
			}
		);
	}
}

void BDSAlgorithm::PrintAndCheck(){
	for (int v = 0; v < n; v++){
		cout<<v<<": ";
		for (int e : adj[v])
			cout<<"(" << (e_v[e] + e_u[e] - v) << ", " << cost[e] << ", " << lp[e] << ", " << in_sol[e] << ") ";
		cout<<" | " << link[v].first << ", " << link[v].second << endl;

		int matched = 1 - cost[adj[v][0]];

		for (int i = 1 + matched; i < adj[v].size(); i++)
			assert(sign(lp[adj[v][i - 1]] - lp[adj[v][i]]) >= 0);
	}
}


/*
	DFS Step of BDS Algorithm. The next tree edge is chosen by the following criteria.
		- Matching edge
		- Heavy edge that maximizes x^*_e
	
	O(n + m)
*/
void BDSAlgorithm::Dfs(int v){
	in[v] = clk++;

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
	Unit edge tree covering

	O(n + |L|)
*/
int BDSAlgorithm::UpLinkCover(int v){
	int my_val = 0;
	for (auto x : tree_adj[v]){
		my_val += UpLinkCover(x);
		
		if (link[x].first != -1)
			if ( link[v].first == -1 or in[ link[x].second ] < in[ link[v].second ] )
				link[v] = link[x];
	}

	my_val += covered[v];

	if (my_val == 0){
		assert(link[v].first != -1);

		covered[ link[v].second ] += -1;
		my_val += 1;
		in_sol[ link[v].first ] = 1;
	}

	return my_val;
}


void BDSAlgorithm::UpLinkAugmentation(int root){
	for (int v = 0; v < n; v++){
		link[v] = {-1, v};
		covered[v] = 0;
	}

	for (int i = 0; i < m; i++) // Preprocessing edges
		if (!in_sol[i] and (sign(lp[i]) > 0)){
			int u = e_u[i];
			int v = e_v[i];
			if (Dec(u, v)) // u is the lower vertex
				swap(u, v);

			if (in[link[u].second] > in[v]){
				link[u].first = i;
				link[u].second = v;
			}
		}

	for (auto u : tree_adj[root])
		UpLinkCover(u);
}

/*
	Run BDS algorithm
*/
void BDSAlgorithm::Run(ListGraph::EdgeMap<int> &_cost, 
	ListGraph::EdgeMap<bool> &BDSSol, 
	ListGraph::EdgeMap<double> &FracSol, 
	ListGraph &G){
	
	Update(_cost, FracSol, G);	

	int wrst = 0;
	vector<int> wrst_sol(m, 0);

	for (int r = 0; r < n; r++){
		// Step 1, find a DFS Tree
		for (int v = 0; v < n; v++){
			tree_adj[v].clear();
			parent[v] = v;
			in[v] = 0;
		}

		clk = 1;
		Dfs(r);

		// Sanity check, found a tree
		assert(parent[r] == r);
		for (int v = 0; v < n; v++)
			if (v != r)
				assert(parent[v] != v);

		// Step 2, uplink only augmentation
		UpLinkAugmentation(r);

		int cur = 0;
		for (int i = 0; i < m; i++)
			cur += in_sol[i];

		if (cur > wrst){
			wrst = cur;
			for (int i = 0; i < m; i++)
				wrst_sol[i] = in_sol[i];
		}

	}


	if (__verbose_mode){
		cout << "BDS Tree Found" << endl;
		for (int v = 0; v < n; v++)
			cout << parent[v] << " <-- " << v << endl;	
		PrintAndCheck();
	}

	for (ListGraph::EdgeIt e(G); e != INVALID; ++e)
		BDSSol[e] = wrst_sol[G.id(e)];


	// Sanity check, checks if edges are from the support
	for (ListGraph::EdgeIt e(G); e != INVALID; ++e)
		if (BDSSol[e] and (sign(FracSol[e]) <= 0))
			assert(0);


	// Sanity check, checks if BDS returned a feasible solution
	ListGraph::NodeMap<bool> ones(G, 1);

	SubGraph<ListGraph> H(G, ones, BDSSol);
	assert(biEdgeConnected(H) == 1);
}
