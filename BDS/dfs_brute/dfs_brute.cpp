#include <iostream>
#include <algorithm>
#include <vector>
#include <map>
#include <unordered_map>
#include <tuple>

using namespace std;

#define all(x) x.begin(),x.end()
#define dbg(x) cout<<#x<<" = "<<x<<endl;
#define chapa cout<<"oi meu chapa"<<endl;

const int MAXN = 10;

int n, m, eid;
int cost[MAXN * MAXN], a[MAXN * MAXN], b[MAXN * MAXN];
double lp[MAXN * MAXN];
vector< pair<int, int> > adj[MAXN];

void ReadGraph(){
	cin>>n>>m;

	eid = 0;
	for (int i = 0; i < n; i++)
		adj[i].clear();

	for (int i = 0; i < m; i++){
		int v, u, c;
		double x;

		cin>>v>>u>>c>>x;

		adj[v].push_back(make_pair(u, eid));
		adj[u].push_back(make_pair(v, eid));

		a[eid] = v;
		b[eid] = u;
		cost[eid] = c;
		lp[eid] = x;

		eid++;
	}
}


int max_cost, min_cost, clk;
int pre[MAXN], in[MAXN], out[MAXN], in_sol[MAXN * MAXN];
vector<int> tree_adj[MAXN];

void DFS(int v){
	in[v] = clk++;
	int nxt;
	
	// dbg(v);

	do{
		nxt = -1;

		for (auto x : adj[v]){
			int u = x.first;
			int id = x.second;

			if (pre[u] != -1)
				continue;

			if ((cost[id] == 0) or (nxt == -1))
				nxt = id;
			else if ((cost[nxt] == 1) and (lp[id] > lp[nxt]))
				nxt = id;

		}

		if (nxt != -1){
			dbg(lp[nxt]);
			assert(lp[nxt] > 0.01);
			int u = a[nxt] + b[nxt] - v;

			in_sol[nxt] = 1;
			tree_adj[v].push_back(u);
			pre[u] = v;

			DFS(u);
		}
	} while(nxt != -1);

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

map<pair<int, long long int>, bool> done; //  may just use a massive matrix


void AllDFS(){
	// cout<<"oi"<<endl;
	for (int v = 0; v < n; v++){
		
		for (int u = 0; u < n; u++){
			pre[u] = -1;
			memo[u] = -1;

			tree_adj[u].clear();
		}
		
		cout<<"new"<<endl;

		for (int u = 0; u < n; u++){
			cout<<u<<": ";
			for (auto x : adj[u])
				cout<<"("<<x.first<<", "<<lp[x.second]<<") ";
			cout<<endl;

		}

		dbg(v);

		clk = 0;
		pre[v] = v;
		DFS(v);
		// cout<<"ahhhh"<<endl;

		long long int msk;
		for (int i = 0; i < m; i++)
			if (in_sol[i])
				msk |= (1<<i);
		if (done[make_pair(v, msk)])
			continue;	
		done[make_pair(v, msk)] = 1;

		int cur_cost = 0;

		for (int u : tree_adj[v])
			cur_cost += Augmentation(u);

		// for (int u = 0; u < n; u++){
		// 	dbg(memo[u]);
		// }

		// exit(0);

		for (int i = 0; i < m; i++){
			cur_cost += in_sol[i] * cost[i];
			in_sol[i] = 0;
		}

		// dbg(cur_cost);

		max_cost = max(max_cost, cur_cost);
		min_cost = min(min_cost, cur_cost);
	}

	// cout<<min_cost<<' '<<max_cost<<endl;
}


int zero[MAXN];

void Backtracking(int v){
	if (v == n){
		AllDFS();
		return;
	}


	sort(adj[v].begin() + 1, adj[v].begin() + (adj[v].size() - zero[v]));

	do{
		Backtracking(v + 1);
		// chapa;
	} while (next_permutation(adj[v].begin() + 1, adj[v].begin() + (adj[v].size() - zero[v])));
}


signed main(){
	ReadGraph();
	max_cost = 0; 
	min_cost = m + 1;


	for(int v = 0; v < n; v++){

		int p = 0;
		for (int i = 0; i < adj[v].size(); i++)
			if (cost[adj[v][i].second] == 0)
				p = i;
		
		swap(adj[v][0], adj[v][p]);

		for (int i = (int)adj[v].size() - 1; i > 0; i--)
			for (int j = i; j > 0; j--)
				if (lp[ adj[v][j].second ] < 0.01){
					swap(adj[v][i], adj[v][j]);
					break;
				}

		for (int i = 0; i < adj[v].size(); i++)
			if (lp[ adj[v][i].second ] < 0.01)
				zero[v]++;		


		// matching edge is always the first	
	}


	Backtracking(0);


	double sol = 0;
	for (int i = 0; i < m; i++)
		sol += lp[i] * cost[i];

	cout<<min_cost<<' '<<max_cost<<endl;
	dbg(sol);
	dbg(min_cost/sol);
	dbg(max_cost/sol);
}
