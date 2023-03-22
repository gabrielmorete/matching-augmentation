#include "main.h"
#include "lemon.h"

/*
	I need lemon to find global min cut in a eficient and stable way
	For bds I can implement it more efficiently without lemon.
*/

class BDS{
	public:
		int n, m, clk;
		vector<vector<int>> adj;
		vector<bool> cost, in_sol;
		vector<int> e_u, e_v, in, out, parent, memo;
		vector<double> lp;
		bool updated;

		BDS(ListGraph &G){
			n = countNodes(G);
			m = countEdges(G);

			adj.resize(n);
			in.resize(n);
			parent.resize(n);
			out.resize(n);
			memo.resize(n);

			cost.resize(m);
			in_sol.resize(m);
			lp.resize(m);


			for (ListGraph::EdgeIt e(G); e != INVALID; ++e){
				int eid = G.id(e);
				int e_u[eid] = G.id(G.u(e));
				int e_v[eid] = G.id(G.v(e));
				adj[v].pb(eid);
				adj[u].pb(eid);
			}
			updated = false;
		}

		void update(ListGraph::EdgeMap<int> &_cost, ListGraph::EdgeMap<double> &_FracSol, ListGraph &G){
			updated = true;

			for (ListGraph::EdgeIt e(G); e != INVALID; e++){
				int eid = G.id(e);
				lp[eid] = FracSol[e];
				cost[eid] = _cost[e];
				in_sol[eid] = 0; 
			}

			// Matched edge is the first
			for (int v = 0; v < n; v++){
				int p = 0;

				for (int e = 0; e < adj[v].size(); e++)
					if (cost[e] == 0 and sign(lp[e]) > 0) // If matched edge is zero, skip
						p = e;
				swap(adj[v][0], adj[v][e]);	
			}

			int matched = 1 - cost[adj[v][0]];

			sort(adj[v].begin() + matched, adj[v].end(),
				[](int a, int b){ // Sort by increase lp value, skip matching edge
					return lp[a] > lp[v];
				}
			)
		}

		/*
			DFS Step of BDS Algorithm. The next tree edge is chosen by the following criteria.
				- Matching edge
				- Heavy edge that maximizes x^*_e
			
			O(n + m)
		*/

		void dfs(int v){
			in[v] = clk++;

			for (int e : adj[v]){
				if (sign(lp[e]) <= 0) // edges are sorted
					break;

				int u = e_u[e] + e_v[e] - v;
				if (parent[u] == u and u != 0){
					parent[u] = v;
					in_sol[e] = 1;
					tree_adj[v].push_back(u);

					dfs(u);
				}
			}

			out[v] = clk++;
		}

		/*
			Retuns true if v is a proper descendent of u.
		*/
		inline bool StrictDec(int u, int v){ // is v strict dec of u
			return (in[u] < in[v]) and (out[v] < out[u]);
		}

		/*
			Retuns true if v is a descendent of u.
		*/
		inline bool Dec(int u, int v){ // is v strict dec of u
			return (in[u] <= in[v]) and (out[v] <= out[u]);
		}


		void run(ListGraph::EdgeMap<bool> &BDSSol; ListGraph &G){
			assert(updated);
			updated = 0;


			// Step 1, find a DFS Tree
			for (int v = 0; v < n; v++)
				parent[v] = v;

			clk = 0;
			dfs(0);

			if (__verbose_mode){
				cout << "BDS Tree Found" << endl;
				for (int v = 0; v < n; v++)
					cout << parent[v] << " <-- " << v << endl;	
			}


			// Step 2, uplink only augmentation
			UpLinkAugmentation(cost, BDSSol, FracSol, parent, in, out, T, G);


			for (ListGraph::EdgeIt e(G); e != INVALID; e++)
				BDSSol[e] = in_sol[G.id(e)];

			ListGraph::NodeMap<bool> ones(G, 1);

			// Sanity check, checks if BDS returned a feasible solution
			SubGraph<ListGraph> H(G, ones, BDSSol);
			assert(biEdgeConnected(H) == 1);

			// Sanity check, checks if edges are from the support
			for (ListGraph::EdgeIt e(G); e != INVALID; ++e)
				if (BDSSol[e] and (sign(FracSol[e]) <= 0))
					assert(0);

	}

	protected:


};


/*
	Dynamic Programming algorithm to solving the up-link only 
	agumentation problem.

	O(n^2|L|)
*/
int UpLinkDP(Node v, ListGraph::NodeMap<int> &memo, 
	ListGraph::NodeMap<Edge> &dp_edge,
	ListGraph::EdgeMap<int> &cost,
	ListGraph::NodeMap<ListGraph::Node> &parent, 
	ListGraph::NodeMap<int> &in, 
	ListGraph::NodeMap<int> &out,
	ListGraph::EdgeMap<bool> &BDSSol,
	ListGraph::EdgeMap<double> &FracSol,
	SubGraph<ListGraph> &T,
	ListGraph &G){
	
	if (memo[v] != -1)
		return memo[v];

	memo[v] = countEdges(G) + 1; // infinity

	for (EdgeIt e(G); e != INVALID; ++e)
		if ((BDSSol[e] == 0) and (sign(FracSol[e]) > 0)){ // A backedge of the support
			
			Node u = G.u(e);
			Node w = G.v(e);
			if (Dec(w, u, in, out)) // If u is a descendent of w
				std::swap(w, u);

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
							subtree_cost += UpLinkDP(y, memo, dp_edge, cost, parent, in, out, BDSSol, FracSol, T, G);
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

// void UpLinkDP(Node v, 
// 	ListGraph::NodeMap<int> &memo, 
// 	ListGraph::NodeMap<Edge> &dp_edge, 
// 	ListGraph::EdgeMap<int> &link_cost, 
// 	ListGraph::NodeMap<ListGraph::Node> &parent, 
// 	ListGraph::NodeMap<int> &in, 
// 	ListGraph::NodeMap<int> &out,
// 	ListGraph::EdgeMap<bool> &BDSSol,
// 	SubGraph<ListGraph> &T){
	
// 	for (SubGraph<ListGraph>::OutArcIt a(T, v); a != INVALID; ++a){
// 		ListGraph::Node y = T.target(a);

// 		if (y != parent[v]){
// 			// cout<<" Transition "<<G.id(v) + 1<<' '<<G.id(y) + 1<<endl;
// 			UpLinkDP(y, memo, dp_edge, link_cost, parent, in, out, BDSSol, T);
// 			}
// 	}		
	
// 	memo[v] = countEdges(G) + 1; // infinity

// 	vector<ListGraph::Edge> prop; 

// 	for (EdgeIt e(G); e != INVALID; ++e)
// 		if (BDSSol[e] == 0){ // A backedge
			
// 			Node u = G.u(e);
// 			Node w = G.v(e);
// 			if (Dec(w, u, in, out)) // If u is a descendent of w
// 				std::swap(w, u);

// 			// u is the ancestor node, w is the descendent

// 			if (StrictDec(u, v, in, out) and Dec(v, w, in, out)){ // feasible link
			
// 				if (link_cost[e] + cost[e] < memo[v]){
// 					// cout<<"Foud New Solution to "<<G.id(v) + 1<<' '<<G.id(u) + 1<<' '<<G.id(w) + 1<<endl;
// 					memo[v] = link_cost[e] + cost[e];
// 					dp_edge[v] = e;
// 				}
// 			}
// 			else if (StrictDec(u, parent[v], in, out) and Dec(parent[v], w, in, out) and !Dec(parent[v], w, in, out)){ // feasible link
// 				prop.push_back(e);
// 			}
// 		}

// 	if (parent[v] != v)	
// 		for (ListGraph::Edge e : prop)
// 			link_cost[e] += memo[v];

// 	// return memo[v];
// }

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

	Node u = T.u(e);
	Node w = T.v(e);
	if (Dec(w, u, in, out)) // If u is a decendent of w
		std::swap(w, u);

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
	ListGraph::EdgeMap<int> &cost,
	ListGraph::EdgeMap<bool> &BDSSol,
	ListGraph::EdgeMap<double> &FracSol,
	ListGraph::NodeMap<ListGraph::Node> &parent, 
	ListGraph::NodeMap<int> &in,
	ListGraph::NodeMap<int> &out,
	SubGraph<ListGraph> &T,
	ListGraph &G){

	ListGraph::NodeMap<int> memo(G, -1); 
	ListGraph::NodeMap<ListGraph::Edge> dp_edge(G); 
	// ListGraph::EdgeMap<int> link_cost(G); 


	for (SubGraph<ListGraph>::OutArcIt a(T, G.nodeFromId(0)); a != INVALID; ++a){
		UpLinkDP(T.target(a), memo, dp_edge, cost, parent, in, out, BDSSol, FracSol, T, G);
		RecoverUpLinkSol(T.target(a), dp_edge, parent, in, out, BDSSol, T);
	}
}

/*
	Implementation of the BDS MAP algorithm.
	Must receive a extreme point solution of the cut LP.

	O(n^2|L|)
*/
