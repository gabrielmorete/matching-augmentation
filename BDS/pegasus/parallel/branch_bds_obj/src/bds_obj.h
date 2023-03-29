#ifndef BDS_DEF
#define BDS_DEF

#include <vector>
#include <cassert>
#include <set>
#include <algorithm>
#include <cmath>
#include "lemon.h"

using namespace std;

class BDSAlgorithm{
	protected:
		int n, m, clk;
		vector<bool> cost, in_sol;
		vector<int> e_u, e_v, memo_edge; 
		vector<int> in, out, parent, dp_edge;
		vector<double> lp;

		vector<vector<int>> adj, tree_adj, cover;
		
		/*
			DFS Step of BDS Algorithm. The next tree edge is chosen by the following criteria.
				- Matching edge
				- Heavy edge that maximizes x^*_e
			
			O(n + m)
		*/

		void Dfs(int v);

		/*
			Retuns true if v is a proper descendent of u.
		*/
		inline bool StrictDec(int u, int v); // is v strict dec of u


		/*
			Retuns true if v is a descendent of u.
		*/
		inline bool Dec(int u, int v); // is v strict dec of u
		

		/*
			Dynamic programming on a tree, follows reverse topological sort.

			memo_edge is ists original cost + the cost of covering every danglind subtree.
			At time it is coveding {v, parent[v]}, the cost of the esdge ist he same as
			memo[v]

			O(n|L|)
		*/
		void UpLinkDP(int v);

		/*
			Algorithm to recover the uplink DP solution.

			O(n + m)
		*/
		void RecoverUpLinkSol(int v);


		void UpLinkAugmentation();

	public:
		BDSAlgorithm();

		/*
			Builds structure for new graph
		*/
		BDSAlgorithm(ListGraph &G);
		/*
			Update costs and lp value
		*/
		void Update(ListGraph::EdgeMap<int> &_cost, ListGraph::EdgeMap<double> &_FracSol, ListGraph &G);

		/*
			Run BDS algorithm
		*/
		void Run(ListGraph::EdgeMap<int> &_cost, 
			ListGraph::EdgeMap<bool> &BDSSol, 
			ListGraph::EdgeMap<double> &FracSol, 
			ListGraph &G);
};



#endif