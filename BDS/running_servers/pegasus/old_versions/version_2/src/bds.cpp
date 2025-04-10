#include "main.h"
#include "lemon.h"


/*
	DFS Step of BDS Algorithm. The next tree edge is chosen by the following criteria.
		- Matching edge
		- Heavy edge that maximizes x^*_e
	
	O(n^2)
*/
void BDSDFS(ListGraph::Node v, 
	ListGraph::NodeMap<ListGraph::Node> &parent, 
	ListGraph::EdgeMap<double> &FracSol, 
	ListGraph::EdgeMap<bool> &BDSSol,
	int &clk,
	ListGraph::NodeMap<int> &in,
	ListGraph::NodeMap<int> &out){

	in[v] = clk++;
	
	bool found_next = 0;
	do {
		ListGraph::Arc next_arc;

		found_next = 0;

		// Here we iterate through the arcs leaving node v.
		for (ListGraph::OutArcIt a(G, v); a != INVALID; ++a){
			ListGraph::Node u = G.target(a);

			if (sign(FracSol[a]) <= 0) // Run algorithm on the support
				continue;

			if (parent[u] == u and G.id(u) != 0){ // Unvisited non root node

				if ((!found_next) or (cost[a] == 0)){
					found_next = 1;
					next_arc = a;
				}
				
				else if (cost[next_arc] == 1 and FracSol[next_arc] < FracSol[a])
						next_arc = a;
			}
		}

		if (found_next){
			BDSSol[next_arc] = 1;
			
			ListGraph::Node u = G.target(next_arc);

			parent[u] = v;

			BDSDFS(u, parent, FracSol, BDSSol, clk, in, out);
		}

	} while (found_next);

	out[v] = clk++;
}

/*
	Retuns true if v is a proper descendent of u.
*/
bool StrictDec(Node u, Node v, ListGraph::NodeMap<int> &in, ListGraph::NodeMap<int> &out){ // is v strict dec of u
	return (in[u] < in[v]) and (out[v] < out[u]);
}

/*
	Retuns true if v is a descendent of u.
*/
bool Dec(Node u, Node v, ListGraph::NodeMap<int> &in, ListGraph::NodeMap<int> &out){ // is v strict dec of u
	return (in[u] <= in[v]) and (out[v] <= out[u]);
}

/*
	Dynamic Programming algorithm to solving the up-link only 
	agumentation problem.

	O(n^2|L|)
*/
int UpLinkDP(Node v, ListGraph::NodeMap<int> &memo, 
	ListGraph::NodeMap<Edge> &dp_edge, 
	ListGraph::NodeMap<ListGraph::Node> &parent, 
	ListGraph::NodeMap<int> &in, 
	ListGraph::NodeMap<int> &out,
	ListGraph::EdgeMap<bool> &BDSSol,
	ListGraph::EdgeMap<double> &FracSol,
	SubGraph<ListGraph> &T){
	
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
							subtree_cost += UpLinkDP(y, memo, dp_edge, parent, in, out, BDSSol, FracSol, T);
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

	Node u = G.u(e);
	Node w = G.v(e);
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
	SubGraph<ListGraph> &T,
	ListGraph::EdgeMap<bool> &BDSSol,
	ListGraph::EdgeMap<double> &FracSol,
	ListGraph::NodeMap<ListGraph::Node> &parent, 
	ListGraph::NodeMap<int> &in,
	ListGraph::NodeMap<int> &out){

	ListGraph::NodeMap<int> memo(G, -1); 
	ListGraph::NodeMap<Edge> dp_edge(G); 
	// ListGraph::EdgeMap<int> link_cost(G); 


	for (SubGraph<ListGraph>::OutArcIt a(T, G.nodeFromId(0)); a != INVALID; ++a){
		UpLinkDP(T.target(a), memo, dp_edge, parent, in, out, BDSSol, FracSol, T);
		RecoverUpLinkSol(T.target(a), dp_edge, parent, in, out, BDSSol, T);
	}
}

/*
	Implementation of the BDS MAP algorithm.
	Must receive a extreme point solution of the cut LP.

	O(n^2|L|)
*/
void BDSAlgorithm(ListGraph::EdgeMap<double> &FracSol, ListGraph::EdgeMap<bool> &BDSSol){
	int n = countNodes(G);

	// Step 1, find a DFS Tree
	ListGraph::NodeMap<ListGraph::Node> parent(G);
	ListGraph::NodeMap<int> in(G), out(G); 

	for (ListGraph::NodeIt v(G); v != INVALID; ++v)
		parent[v] = v;

	int cnt = 0;
	BDSDFS(G.nodeFromId(0), parent, FracSol, BDSSol, cnt, in, out);

	if (__verbose_mode){
		cout << "BDS Tree Found" << endl;
		for (NodeIt v(G); v != INVALID; ++v)
			cout << G.id(parent[v]) << " <-- " << G.id(v) << endl;	
	}


	ListGraph::NodeMap<bool> ones(G, 1); // Subgraph must be spanning

	ListGraph::EdgeMap<bool> tree_edges(G); // Graph adaptor changes with edge map
	for (ListGraph::EdgeIt e(G); e != INVALID; ++e)
		tree_edges[e] = BDSSol[e];

	// Build DFS Tree
	SubGraph<ListGraph> T(G, ones, tree_edges);
	
	assert(biEdgeConnected(G) == 1);
	assert(connected(T) == 1);


	// Step 2, uplink only augmentation
	UpLinkAugmentation(T, BDSSol, FracSol, parent, in, out);

	// Sanity check, checks if BDS returned a feasible solution
	SubGraph<ListGraph> H(G, ones, BDSSol);
	assert(biEdgeConnected(H) == 1);

	// Sanity check, checks if edges are from the support
	for (ListGraph::EdgeIt e(G); e != INVALID; ++e)
		if (BDSSol[e] and (sign(FracSol[e]) <= 0))
			assert(0);
}