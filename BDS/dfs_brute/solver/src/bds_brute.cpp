#include "main.h"
#include "lemon.h"

// #define all(x) x.begin(),x.end()
// #define dbg(x) cout << #x << " = " << x << endl
// #define chapa cout<<"oi meu chapa"<<endl;


const int MAXN = 20;

struct bds_brute{
	int n, m, eid;
	int is_matched[MAXN];
	int edge_cost[MAXN * MAXN], a[MAXN * MAXN], b[MAXN * MAXN];
	double lp[MAXN * MAXN];
	vector< pair<int, int> > adj[MAXN];
	string name;

	int max_cost, min_cost, clk;
	int pre[MAXN], in[MAXN], out[MAXN], in_sol[MAXN * MAXN];
	vector<int> tree_adj[MAXN];

	void DFS(int v){
		in[v] = clk++;
		int nxt;
		
		// Edges are sorted by matching + decreasing value of lp
		// Sanity check
		for (int i = 1 + is_matched[v]; i < adj[v].size(); i++)
			assert(sign(lp[adj[v][i - 1].second] - lp[adj[v][i].second]) >= 0);
		
		for (auto x : adj[v]){
			int u = x.first;
			int id = x.second;

			if (pre[u] != -1)
				continue;

			pre[u] = v;
			in_sol[id] = 1;
			tree_adj[v].push_back(u);
			DFS(u);
		}

		out[v] = clk++;
	}

	/*
		Retuns true if v is a proper descendent of u.
	*/
	bool StrictDec(int u,int v){ // is v strict dec of u
		return (in[u] < in[v]) and (out[v] < out[u]);
	}

	/*
		Retuns true if v is a descendent of u.
	*/
	bool Dec(int u, int v){ // is v strict dec of u
		return (in[u] <= in[v]) and (out[v] <= out[u]);
	}

	int memo[MAXN];

	int Augmentation(int v){
		if (memo[v] != -1)
			return memo[v];

		memo[v] = m + 1; // INF

		for (int i = 0; i < m; i++){
			if (in_sol[i])
				continue;

			int u = a[i];
			int w = b[i];
			if (Dec(w, u))
				swap(u, w); // u is the ancestor

			if (StrictDec(u, v) and Dec(v, w)){ // feasible link
				int cur_cost = edge_cost[i];

				int lst = w;
				while (w != pre[v]){
					for (int x : tree_adj[w])
						if (x != lst)
							cur_cost += Augmentation(x);

					lst = w;
					w = pre[w];
				}
		
				memo[v] = min(memo[v], cur_cost);
			}
		}

		return memo[v];
	}

	map<long long int, bool> done[MAXN]; //  may just use a massive matrix
	// set<pair<int, long long int>> done;


	void AllDFS(){
		for (int v = 0; v < n; v++){ // For every root
		
			for (int u = 0; u < n; u++){ // Clear data
				pre[u] = -1;
				memo[u] = -1;

				tree_adj[u].clear();
			}

			for (int i = 0; i < m; i++)
				in_sol[i] = 0;


			clk = 0;
			pre[v] = v;
			DFS(v); // Generate de DFS tree with root v

			long long int msk = 0;
			for (int i = 0; i < m; i++) // bitmask fo the tree
				if (in_sol[i])
					msk |= (1ll<<i);

			if (done[v][msk] == 0){ // new tree
				done[v][msk] = 1;
				
				int cur_cost = 0;

				for (int u : tree_adj[v])
					cur_cost += Augmentation(u);

				for (int i = 0; i < m; i++)
					cur_cost += in_sol[i] * edge_cost[i];

				max_cost = max(max_cost, cur_cost);
				min_cost = min(min_cost, cur_cost);
			}
		}
	}

	int values[MAXN], sz[MAXN][MAXN];
	void Backtracking(int v, int stage){
		if (v == n){
			AllDFS();
			return;
		}

		if (stage == values[v] + 1)
			return Backtracking(v + 1, 1);

		sort(adj[v].begin() + sz[v][stage - 1], adj[v].begin() + sz[v][stage]);

		do {
			Backtracking(v, stage + 1);
		} while (next_permutation(adj[v].begin() + sz[v][stage - 1], adj[v].begin() + sz[v][stage]));
	}

	void print(){
		for (int u = 0; u < n; u++){
			cout<<u<<": ";
			for (auto x : adj[u])
				cout<<"("<<x.first<<", "<<lp[x.second]<<") ";
			cout<<endl;
		}
	}

	void ReadGraph(ListGraph::EdgeMap<double> &FracSol){
		n = countNodes(G);
		m = countEdges(G);

		eid = 0;
		for (int i = 0; i < n; i++)
			adj[i].clear();

		for (int i = 0; i < m; i++){
			ListGraph::Edge ed = G.edgeFromId(i);
			int v = G.id(G.v(ed));
			int u = G.id(G.u(ed));

			adj[v].push_back(make_pair(u, eid));
			adj[u].push_back(make_pair(v, eid));

			a[eid] = v;
			b[eid] = u;
			edge_cost[eid] = cost[ed];
			dbg(FracSol[ed]);
			lp[eid] = FracSol[ed];

			eid++;
		}
	}


	pair<int, int> brute_all(ListGraph::EdgeMap<double> &FracSol){
		
		ReadGraph(FracSol);
		
		max_cost = 0; 
		min_cost = m + 1;


			for(int v = 0; v < n; v++){

		// Matching edge will be the first edge
		int p = 0;
		for (int i = 0; i < adj[v].size(); i++)
			if (edge_cost[adj[v][i].second] == 0)
				p = i;
		
		swap(adj[v][0], adj[v][p]);

		// Edges of the adj will be sorted by decreasing LP value.
		// We will only permute edges with the same value

		is_matched[v] = 1 - edge_cost[ adj[v][0].second ];

		sort(adj[v].begin() + is_matched[v], adj[v].end(),
			[this](pair<int, int> a, pair<int, int> b){
				return lp[a.second] > lp[b.second];
			}
		);

		while ((adj[v].size() > 1) and ( sign( lp[adj[v].back().second] ) <= 0 ))
			adj[v].pop_back();


		// sz[v][i] will be the amount of edges with the ith lp value
		// sz[v][0] = 0 since we will acumulate the sum

		sz[v][0] = 0;

		values[v] = 1;
		sz[v][1] = 1;

		if (is_matched[v]){ // Matching edge is different
			values[v] = 2;
			sz[v][2] = 1;
		}

		for (int i = 1 + is_matched[v]; i < adj[v].size(); i++){
			int e1 = adj[v][i].second;
			int e2 = adj[v][i - 1].second;

			if (sign(lp[e1] - lp[e2]) != 0){
				values[v]++;
				sz[v][values[v]] = 0;
			}

			sz[v][values[v]]++;
		}
		

		for (int i = 1; i <= values[v]; i++)
			sz[v][i] += sz[v][i - 1];

		// cout<<v<<" : ";
		// for (int i = 0; i < values[v]; i++)
		// 	cout<<sz[v][i]<<' ';
		// cout<<endl;
	}
	

	Backtracking(0, 1);


		return make_pair(min_cost, max_cost);
	}
};
