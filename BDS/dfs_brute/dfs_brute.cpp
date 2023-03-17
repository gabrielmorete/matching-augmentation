#include <iostream>
#include <algorithm>
#include <vector>
#include <map>
#include <unordered_map>
#include <set>
#include <tuple>
#include <cassert>
#include <fstream>

using namespace std;

#define all(x) x.begin(),x.end()
#define dbg(x) cout << #x << " = " << x << endl
#define chapa cout<<"oi meu chapa"<<endl;

const double EPS = 1e-3;
int sign(double x) { return (x > EPS) - (x < -EPS); }

const int MAXN = 20;

int n, m, eid;
int cost[MAXN * MAXN], a[MAXN * MAXN], b[MAXN * MAXN];
double lp[MAXN * MAXN];
vector< pair<int, int> > adj[MAXN];
string name;
void ReadGraph(){
	cin>>n>>m;

	eid = 0;
	for (int i = 0; i < n; i++)
		adj[i].clear();

	for (int i = 0; i < m; i++){
		int v, u, c;
		double x;

		cin>>v>>u>>c>>x;

		// cout<<"--"<<v<<' '<<u<<' '<<c<<' '<<x<<endl;

		adj[v].push_back(make_pair(u, eid));
		adj[u].push_back(make_pair(v, eid));

		a[eid] = v;
		b[eid] = u;
		cost[eid] = c;
		lp[eid] = x;

		eid++;
	}

	cin>>name;
}


int max_cost, min_cost, clk;
int pre[MAXN], in[MAXN], out[MAXN], in_sol[MAXN * MAXN];
vector<int> tree_adj[MAXN];

void DFS(int v){
	in[v] = clk++;
	int nxt;
	
	// Edges are sorted by matching + decreasing value
	for (int i = 2; i < adj[v].size(); i++)
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
			int cur_cost = cost[i];

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
	bool new_tree = 0;
	for (int v = 0; v < n; v++){
	
		for (int u = 0; u < n; u++){
			pre[u] = -1;
			memo[u] = -1;

			tree_adj[u].clear();
		}

		for (int i = 0; i < m; i++)
			in_sol[i] = 0;


		clk = 0;
		pre[v] = v;
		DFS(v);

		long long int msk = 0;
		for (int i = 0; i < m; i++)
			if (in_sol[i])
				msk |= (1ll<<i);
		
		if (done[v][msk] == 0){
			done[v][msk] = 1;
			
			new_tree = 0;

			int cur_cost = 0;

			for (int u : tree_adj[v])
				cur_cost += Augmentation(u);

			for (int i = 0; i < m; i++)
				cur_cost += in_sol[i] * cost[i];

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

	if (stage == values[v])
		return Backtracking(v + 1, 1);

	sort(adj[v].begin() + sz[v][stage - 1], adj[v].begin() + sz[v][stage]);

	do {
		Backtracking(v, stage + 1);
	} while (next_permutation(adj[v].begin() + sz[v][stage - 1], adj[v].begin() + sz[v][stage]));
}


signed main(){
	ReadGraph();
	max_cost = 0; 
	min_cost = m + 1;


	for(int v = 0; v < n; v++){

		sort(adj[v].begin(), adj[v].end(), // sort all edges
			[](pair<int, int> a, pair<int, int> b){
				return lp[a.second] > lp[b.second];
			}
		);


		// Matching edge will be the first edge
		int p = 0;
		for (int i = 0; i < adj[v].size(); i++)
			if (cost[adj[v][i].second] == 0)
				p = i;
		
		swap(adj[v][0], adj[v][p]);

		// Edges of the adj will be sorted by decreasing LP value.
		// We will only permute edges with the same value

		sort(adj[v].begin() + 1, adj[v].end(),
			[](pair<int, int> a, pair<int, int> b){
				return lp[a.second] > lp[b.second];
			}
		);

		while ((adj[v].size() > 1) and ( sign( lp[adj[v].back().second] ) <= 0 ))
			adj[v].pop_back();


		values[v] = 1;
		sz[v][0] = 1;

		for (int i = 1; i < adj[v].size(); i++){
			int e1 = adj[v][i].second;
			int e2 = adj[v][i - 1].second;

			if (sign(lp[e1] - lp[e2]) != 0)
				values[v]++;

			sz[v][values[v] - 1]++;
		}
		

		for (int i = 1; i < values[v]; i++)
			sz[v][i] += sz[v][i - 1];

		// cout<<v<<" : ";
		// for (int i = 0; i < values[v]; i++)
		// 	cout<<sz[v][i]<<' ';
		// cout<<endl;
	}
	

	// for (int u = 0; u < n; u++){
	// 	cout<<u<<": ";
	// 	for (auto x : adj[u])
	// 		cout<<"("<<x.first<<", "<<lp[x.second]<<") ";
	// 	cout<<endl;
	// }	

	Backtracking(0, 1);


	double sol = 0;
	for (int i = 0; i < m; i++)
		sol += lp[i] * cost[i];


	// dbg(min_cost);
	// dbg(max_cost);
	// dbg(sol);
	


	double fmin = min_cost/sol;
	double fmax = max_cost/sol;

	cout<<fmin<<' '<<fmax<<' '<<sol<<endl;

	if (sign(fmin - 1.4) >= 0 or sign(fmax - 1.4) > 0){
		ofstream log("log_bds_brute", ios::app); // open in append mode
		log << name << ' ' << fmin << ' ' << fmax << endl;
		log.close();
	}

}
