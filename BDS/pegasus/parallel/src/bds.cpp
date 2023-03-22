#include "main.h"
#include "lemon.h"

/*
	I need lemon to find global min cut in a eficient and stable way
	For bds I can implement it more efficiently without lemon.
*/

class BDS{
	protected:
		int n, m, clk;
		vector<vector<int>> adj;
		vector<bool> cost, in_sol;
		vector<int> e_u, e_v, in, out, parent, memo;
		vector<double> lp;
		bool updated;

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

		vector<vector<int>> cover;

		/*
			Dynamic programming on a tree, follows reverse topological sort
		*/

		int UpLinkDP(int v){
			for (auto x : tree_adj[v])
				UpLinkDP(x);

			dp_edge[v] = cover[v].begin();
			for (int e : cover[v])
				if (memo_edge[e] < memo_edge[dp_edge[v]])
					dp_edge[v] = e;

			if (parent[v] != 0)
				for (int e : cover[parent[v]])
					if (Dec(e_u[e], v) == 0 and Dec(e_v[e], u) == 0)
						memo_edge[e] += memo_edge[ dp_edge[v] ];
		}

		/*
			Algorithm to recover the uplink DP solution.

			O(n + m)
		*/
		void RecoverUpLinkSol(int v){
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


		void UpLinkAugmentation(){
			fill(memo.begin(), memo.end(), -1);

			cover.resize(n); // do this above
			for (int v = 0; v < n; v++)
				cover[v].clear();

			for (int i = 0; i < m; i++)
				if (!in_sol[i]){
					int u = e_u[i];
					int v = e_v[i];
					if (Dec(u, v)) // u is the lower vertex
						swap(u, v);

					while (u != v){ // link i covers {v, parent[v]}
						cover[u].push_back[i];
						u = parent[u];
					}

					memo_edge[i] = cost[i];
				}


			// Preprocess
			for (auto u : tree_adj[0]){
				UpLinkDP(u);
				RecoverUpLinkSol(u);
			}
		}


		/*
			Implementation of the BDS MAP algorithm.
			Must receive a extreme point solution of the cut LP.

			O(n^2|L|)
		*/

	public:


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
			UpLinkAugmentation();


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


};

