#include "bits/stdc++.h"
using namespace std;

int n;
vector<int> deg, lft, rgt, max_deg;
vector< vector<int> > adj;

vector<array<int, 3>> edges;

void add_edge(int u, int v, int c){
	edges.push_back({u, v, c});
	deg[u]++;
	deg[v]++;
}

void remove_edge(){
	auto e = edges.back();
	deg[e[0]]--;
	deg[e[1]]--;
	edges.pop_back();	
}

bool check(){
	for (int i = 1; i <= n; i++)
		if (deg[i] <= 2)
			return false;
	return true;	
}

int cnt;

void print(){
	cout << n << ' ' << edges.size() << endl;
	for (auto x : edges){
		cout << x[0] - 1<< ' ' << x[1] - 1<< ' ' << x[2] << endl;
	}
	cout << cnt++ << endl << endl;
}


vector<int> order;
void run(int p, int q){
	// cout << p << ' ' << q << ' ' << order.size() << ' ' << adj[order[p]].size() << endl;
	
	if (p == order.size()){
		if (check())
			print();
		return;
	}

	if (q == adj[order[p]].size()){
		if (deg[order[p]] <= 2)
			return;

		run(p + 1, 0);
		return;
	}

	run(p, q + 1);

	if (deg[order[p]] <= max_deg[order[p]] and deg[adj[order[p]][q]] <= max_deg[deg[adj[order[p]][q]]]){
		add_edge(order[p], adj[order[p]][q], 1);
		run(p, q + 1);
		remove_edge();
	}
}



signed main(int argc, char *argv[]){
	if (argc != 2){
		cout << "Usage: number of vertices" << endl;
		return 1;
	}

	string s = argv[1];
	n = stoi(s);

	deg.resize(n + 1);
	adj.resize(n + 1);
	rgt.resize(n + 1);
	lft.resize(n + 1);
	max_deg.resize(n + 1);
	fill(deg.begin(), deg.end(), 0);

	// if (n % 6 != 0){
	// 	cout << "not a multiple of 6" << endl;
	// 	return 1;
	// }

	// add edges from the caterpillar
	for (int i = 1; i <= 6; i += 2){
		add_edge(i, i + 1, 0);
		add_edge(i + 1, i + 2, 1);
	}

	for (int i = 7; i <= n - 4; i += 4){
		add_edge(i, i + 1, 0);
		add_edge(i, i + 2, 1);
		add_edge(i + 2, i + 3, 0);
		add_edge(i + 3, i + 4, 1);
	}

	// final part
	add_edge(n - 3, n - 2, 0);
	add_edge(n - 3, n - 1, 1);
	add_edge(n - 1, n, 0);

	// possible good extra edges

	add_edge(2, n, 1);

	// generate marks

	for (int i = 8; i < n; i += 8)
		lft[i] = 1;

	for (int i = 12; i < n; i += 8)
		rgt[i] = 1;


	// speedup, fix first left and right

	add_edge(8, 3, 1);
	add_edge(8, 5, 1);

	add_edge(16, 7, 1);
	add_edge(16, 13, 1);



	add_edge(12, 2, 1);
	add_edge(12, 4, 1);
	add_edge(12, 6, 1);
	
	// build adjacency
	for (int i = 9; i < n; i++)
		if (lft[i] == 0 and rgt[i] == 0)
			adj[1].push_back(i); // adjacent to 1

	for (int i = 16; i < n; i += 8)
		for (int j = 7; j < i - 2; j++)
			if (lft[j] == 0 and rgt[j] == 0)
				adj[i].push_back(j);

	for (int i = 20; i < n; i += 8)
		for (int j = 7; j < i - 2; j++)
			if (lft[j] == 0 and rgt[j] == 0)
				adj[i].push_back(j);

	for (int i = 7; i < n - 3; i++)
		if (lft[i] == 0 and rgt[i] == 0)
			adj[n - 1].push_back(i);	

	for (int i = 7; i < n - 3; i++)
		if (lft[i] == 0 and rgt[i] == 0)
			adj[n].push_back(i);

	fill(max_deg.begin(), max_deg.end(), 4);


	// for (int i = 16; i < n; i += 8){
	// 	order.push_back(i);
	// }

	for (int i = 20; i < n; i += 8){
		order.push_back(i);
		// max_deg[i]= 4;
	}

	order.push_back(1);
	order.push_back(n - 1);
	order.push_back(n);

	// max_deg[1] = max_deg[2] = 4;
	// for (int i = 7; i < 4; i += 4)
	// 	max_deg[i] = 4;

	run(0, 0);

	cout << cnt << endl;

	// print();
}