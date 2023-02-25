#include <iostream>
#include <fstream>
#include <filesystem>
#include <lemon/list_graph.h>
#include <lemon/nauty_reader.h>

using namespace std;
using namespace lemon;

// Safe handling doubles
const double EPS = 1e-3;
int sign(double x) { return (x > EPS) - (x < -EPS); }


ListGraph G;
ListGraph::EdgeMap<int> cost(G);

bool __found_feasible;
int __cur_graph_id;
ofstream g_out, log_out;

/*
	This function calls the LP, IP and BDS algorithms to
	solve the MAP problem and compares their outputs.
	If the output satisfies the requerements, it writes a
	file.
*/
void SolveCurrentMatching(int matching_id){
	ListGraph::EdgeMap<int> IntSol(G);
	// IntegerSolution(IntSol);

	ListGraph::EdgeMap<double> FracSol(G);
	// FractionalSolution(FracSol);

	ListGraph::EdgeMap<int> BDSSol(G);
	// BDSAlgorithm(FracSol, BDSSol);

	int cost_Int = 0;
	int cost_BDS = 0;
	double cost_Frac = 0;

	for (ListGraph::EdgeIt e(G); e != INVALID; ++e){
		int u = G.id(G.u(e));
		int v = G.id(G.v(e));

		cost_Int +=  IntSol[e] * cost[e];
		cost_Frac +=  FracSol[e] * cost[e];
		cost_BDS +=  BDSSol[e] * cost[e];

		// cout<<u + 1<<' '<<v + 1<<' '<<FracSol[e]<<' '<<IntSol[e]<<' '<<BDSSol[e]<<endl;
	}

	// Found a feasible example, print to file
	if (sign(3.0 * cost_Int - 4.0 * cost_Frac) >= 0 or sign(3.0 * cost_BDS - 4.0 * cost_Frac) >= 0){
		if (__found_feasible == 0){ // First feasible example found
			// create file "g"+cnt
			g_out.open(to_string(countNodes(G)) + "/g" + to_string(__cur_graph_id));
			
			g_out << countNodes(G) <<' ' << countEdges(G) << endl << endl;
			
			for (ListGraph::EdgeIt e(G); e != INVALID; ++e)
				g_out << G.id(G.u(e)) << ' ' << G.id(G.v(e)) << endl;
			
			g_out << "----------" << endl << endl;
		}

		__found_feasible = 1;

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

		g_out << "Frc :" << cost_Frac << " | ";
		for (ListGraph::EdgeIt e(G); e != INVALID; ++e)
			g_out << FracSol[e] << ' ';
		g_out << endl;

		g_out << "Int :" << cost_Int << " | ";
		for (ListGraph::EdgeIt e(G); e != INVALID; ++e)
			g_out << IntSol[e] << ' ';
		g_out << endl;

		g_out << "BDS :" << cost_BDS << " | ";
		for (ListGraph::EdgeIt e(G); e != INVALID; ++e)
			g_out << BDSSol[e] << ' ';
		g_out << endl << endl;

		// Generate entry in the log file
		log_out << "Found feasible example g" << __cur_graph_id << " matching id " << matching_id << endl;
		log_out << "Int/Frc = " << (double) cost_Int/cost_Frac << " BDS/Frc = " << cost_BDS/cost_Frac << endl;
		log_out << endl;
	}
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

	if (__found_feasible == 1)
		g_out << "Number of matchings : " << total_matchings << endl;
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

		/* 
			Since the input data is massive, we will one write a file if
			there is some feasible solution.
		*/
		__found_feasible = 0;
		__cur_graph_id = cnt;

		SolveAllMatchings();
		cnt++;

		if (__found_feasible)
			g_out.close();
	}

	log_out.close();
}


signed main(){
	ListGraph G;

	RunNautyInput();

}