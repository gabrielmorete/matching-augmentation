#include "bits/stdc++.h"
using namespace std;


int n;
vector<array<int, 3>> edges;

void add_edge(int u, int v, int c){
	edges.push_back({u, v, c});
}

void print(){
	cout << n << ' ' << edges.size() << endl;
	for (auto x : edges){
		cout << x[0] << ' ' << x[1] << ' ' << x[2] << endl;
	}
}


signed main(int argc, char *argv[]){
		if (argc != 2){
			cout << "Usage:number of vertices" << endl;
			return 1;
		}

		string s = argv[1];
		n = stoi(s);

		// if (n % 6 != 0){
		// 	cout << "not a multiple of 6" << endl;
		// 	return 1;
		// }


		// add edges from the caterpillar
		for (int i = 1; i <= 6; i += 2){
			add_edge(i, i + 1, 0);
			add_edge(i + 1, i + 2, 1);
		}

		for (int i = 7; i <= n - 3; i += 4){
			add_edge(i, i + 1, 0);
			add_edge(i, i + 2, 1);
			add_edge(i + 2, i + 3, 0);
			add_edge(i + 3, i + 4, 1);
		}

		// add(n - 3, n - 1, 1);
		// add(n - 1, n, 0);


		// // only for n = 22

		// // left side edges

		// add(8, 3, 1);
		// add(8, 5, 1);

		// add(16, 9, 1);
		// add(16, 13, 1);

		// // right side
		// add(12, 2, 1);
		// add(12, 4, 1);
		// add(12, 6, 1);

		// add(20, 10, 1);
		// add(20, 14, 1);
		// add(20, 18, 1);

		// // head
		// add(1, 15, 1)
		// add(1, 18, 1)
		// add(1, 19, 1)


		print();
}