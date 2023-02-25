#include "bits/stdc++.h"
#include <filesystem>
#include <lemon/list_graph.h>
#include <lemon/nauty_reader.h>


using namespace std;
using namespace lemon;

ListGraph G;
ListGraph::EdgeMap<int> cost(G);

ofstream g_out, log_out;

/*
	This function calls the LP, IP and BDS algorithms to
	solve the MAP problem and compares their outputs.
	If the output satisfies the requerements, it writes a
	file.
*/
void SolveCurrentMatching(int matching_id){
	g_out << matching_id << ": ";

	bool first = 1;
	for (ListGraph::EdgeIt e(G); e != INVALID; ++e)
		if (cost[e] == 0){
			if (!first)
				g_out << ", ";

			g_out << G.id(G.u(e)) << " " << G.id(G.v(e));
			first = 0;
		}

	g_out<<endl;	
}


/*
	Backtracking algorithm to find all matchings of G.
	Matched edges are marked with cost 0 on the global
	EdgeMap cost. The running time is exponential
*/
void FindAllMatchings(int e_id, int &m, int &total_matchings, ListGraph::NodeMap<bool> &matched){
	if (e_id >= m){
		SolveCurrentMatching(total_matchings);
		return;
	}

	// Case 1 : won't add edge e_id to the matching
	FindAllMatchings(e_id + 1, m, total_matchings, matched); 

	// Case 2 : if possible, will add e_id to the matching
	ListGraph::Edge e = G.edgeFromId(e_id);
	if ((matched[G.u(e)] == 0) and (matched[G.v(e)] == 0)){ // May add e_id
		
		matched[G.u(e)] = 1;
		matched[G.v(e)] = 1;
		cost[e] = 0;
		total_matchings++;

		FindAllMatchings(e_id + 1, m, total_matchings, matched);

		matched[G.u(e)] = 0;
		matched[G.v(e)] = 0;
		cost[e] = 1;
	}
}


/*
	Wrapper function for the matching backtrackig algorithm.
*/
void SolveAllMatchings(){
	// Initialize all edges to be heavy
	for (ListGraph::EdgeIt e(G); e != INVALID; ++e)
		cost[e] = 1;

	ListGraph::NodeMap<bool> matched(G);

	int total_matchings = 1, m = countEdges(G);
	FindAllMatchings(0, m, total_matchings, matched);

	g_out << endl << "Number of matchings : " << total_matchings << endl;
}



/*
	This functions receiv nauty's geng output from stdin(may modify this),
	build a LEMON graph and the log files, and calls the function
	that iterates through all matchings.
*/
void RunNautyInput(){
	int cnt = 1;
	while (readNautyGraph(G, cin)){
		int n = countNodes(G);
		int m = countEdges(G);

		if (cnt == 1){ // Create folder to log files, create log stream
			filesystem::create_directory("./" + to_string(n));
			log_out.open(to_string(n) + "/log");
		}	

		// Next loop makes shure that lemon graph is consistent with the algorithm input
		int nvtx = n - 1, ok = 1;
		for (ListGraph::NodeIt v(G); v != INVALID; ++v)
			if (G.id(v) != nvtx--){
				cout << "Found inconsistency regarding vertex indexing -- Skip graph number " << cnt << endl;
				ok = 0;
				break;
			}

		if (ok){
			// create file "g"+cnt
			g_out.open(to_string(n) + "/g" + to_string(cnt));

			g_out << n <<' ' << m << endl << endl;

			for (ListGraph::EdgeIt e(G); e != INVALID; ++e)
				g_out << G.id(G.u(e)) << ' ' << G.id(G.v(e)) << endl;

			g_out << "----------" << endl << endl;
		}	


		SolveAllMatchings();
		cnt++;

		g_out.close();
	}

	log_out.close();
}






signed main(){
	ListGraph G;

	RunNautyInput();

}