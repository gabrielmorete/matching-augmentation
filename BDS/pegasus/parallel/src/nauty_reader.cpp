#include "lemon.h"
#include "main.h"

/*
	Nauty Reader
*/

bool __found_feasible;
int __cur_graph_id, __best_IP_graph_id, __best_IP_matching_id, __best_BDS_graph_id, __best_BDS_matching_id;
double __best_IP, __best_BDS;
ofstream g_out, log_out;

/*
	This function calls the LP, IP and BDS algorithms to
	solve the MAP problem and compares their outputs.
	If the output satisfies the requerements, it writes a
	file.
*/
void SolveCurrentMatching(int matching_id,
	ListGraph::EdgeMap<int> &cost,
	GRBModel &frac_model,
	GRBVar *frac_vars,
	GRBModel &int_model,
	GRBVar *int_vars,
	ListGraph &G){

	ListGraph::EdgeMap<bool> BDSSol(G);
	ListGraph::EdgeMap<int> IntSol(G);
	ListGraph::EdgeMap<double> FracSol(G);

	SolveMapInstance(cost, FracSol, IntSol, BDSSol, frac_model, frac_vars, int_model, int_vars, G);

	if (sign(FracSol[G.edgeFromId(0)]) == -1 or IntSol[G.edgeFromId(0)] == -1){
		#pragma omp critical
		{
			ofstream excep(to_string(countNodes(G))+"/exception", std::ios_base::app);
			excep << "Exception on example g" << __cur_graph_id << " matching id " << matching_id << endl;
			excep.close();
		}
		return;
	}


	int cost_Int = 0;
	int cost_BDS = 0;
	double cost_Frac = 0;

	for (ListGraph::EdgeIt e(G); e != INVALID; ++e){
		int u = G.id(G.u(e));
		int v = G.id(G.v(e));

		cost_Int +=  IntSol[e] * cost[e];
		cost_Frac +=  FracSol[e] * cost[e];
		cost_BDS +=  (int)BDSSol[e] * cost[e];
	}

	/* 
		Found a feasible example, print to file
			- IP gap must be at least 6/5
			- BDS gap must be better than 4/3
	*/
	if (sign(5.0 * cost_Int - 6.0 * cost_Frac) >= 0 or sign(3.0 * cost_BDS - 4.0 * cost_Frac) > 0){
		if (__found_feasible == 0){ // First matching found for this graph
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

		g_out << "Frc: " << cost_Frac << " | ";
		for (ListGraph::EdgeIt e(G); e != INVALID; ++e)
			g_out << FracSol[e] << ' ';
		g_out << endl;

		g_out << "Int: " << cost_Int << " | ";
		for (ListGraph::EdgeIt e(G); e != INVALID; ++e)
			g_out << IntSol[e] << ' ';
		g_out << endl;

		g_out << "BDS: " << cost_BDS << " | ";
		for (ListGraph::EdgeIt e(G); e != INVALID; ++e)
			g_out << BDSSol[e] << ' ';
		g_out << endl << endl;

		
		#pragma omp critical
		{
			// Generate entry in the log file
			log_out << "Found feasible example g" << __cur_graph_id << " matching id " << matching_id << endl;
			log_out << "Int/Frc = " << (double) cost_Int/cost_Frac << " BDS/Frc = " << (double) cost_BDS/cost_Frac << endl;
			log_out << endl;
		}
	}

	#pragma omp critical
	{
		if (sign((double)cost_Int/cost_Frac - __best_IP) > 0){
			__best_IP = (double)cost_Int/cost_Frac;
			__best_IP_graph_id = __cur_graph_id;
			__best_IP_matching_id = matching_id;
		}


		if (sign((double)cost_BDS/cost_Frac - __best_BDS) > 0){
			__best_BDS = (double)cost_BDS/cost_Frac;
			__best_BDS_graph_id = __cur_graph_id;
			__best_BDS_matching_id = matching_id;
		}
	}	
}


/*
	Backtracking algorithm to find all matchings of G.
	Matched edges are marked with cost 0 on the global
	EdgeMap cost. The running time is exponential
*/
void FindAllMatchings(int e_id, int &n, int &m, int &n_matched, int &total_matchings, 
	ListGraph::NodeMap<bool> &matched,
	ListGraph::EdgeMap<int> &cost,
	GRBModel &frac_model,
	GRBVar *frac_vars,
	GRBModel &int_model,
	GRBVar *int_vars,
	ListGraph &G){

	if (e_id >= m){
		SolveCurrentMatching(total_matchings, cost, frac_model, frac_vars, int_model, int_vars, G);
		return;
	}

	if (n_matched >= n - 1){ // matching cant increase, prune
		SolveCurrentMatching(total_matchings, cost, frac_model, frac_vars, int_model, int_vars, G);
		return;
	}


	// Case 1 : won't add edge e_id to the matching
	FindAllMatchings(e_id + 1, n, m, n_matched, total_matchings, matched, cost, frac_model, frac_vars, int_model, int_vars, G); 

	// Case 2 : if possible, will add e_id to the matching
	ListGraph::Edge e = G.edgeFromId(e_id);
	if ((matched[G.u(e)] == 0) and (matched[G.v(e)] == 0)){ // May add e_id
		
		matched[G.u(e)] = 1;
		matched[G.v(e)] = 1;
		n_matched += 2;
		cost[e] = 0;
		total_matchings++;

		FindAllMatchings(e_id + 1, n, m, n_matched, total_matchings, matched, cost, frac_model, frac_vars, int_model, int_vars, G);

		matched[G.u(e)] = 0;
		matched[G.v(e)] = 0;
		n_matched -= 2;
		cost[e] = 1;
	}
}


/*
	Wrapper function for the matching backtrackig algorithm.
*/
void SolveAllMatchings(ListGraph &G){
	int n = countNodes(G), m = countEdges(G);

	ListGraph::EdgeMap<int> cost(G); // Cost of the edges
	// Initialize all edges to be heavy
	for (ListGraph::EdgeIt e(G); e != INVALID; ++e)
		cost[e] = 1;

	GRBModel frac_model(env);
	GRBVar frac_vars[m];
	BuildFractional(frac_model, frac_vars, G);

	GRBModel int_model(env);
	GRBVar int_vars[m];
	BuildIntegral(int_model, int_vars, G);
	MinimumCut cb = MinimumCut(int_vars, n, m, G);
	int_model.setCallback(&cb);

	ListGraph::NodeMap<bool> matched(G);

	int total_matchings = 1, n_matched = 0;
	FindAllMatchings(0, n, m, n_matched, total_matchings, matched, cost, frac_model, frac_vars, int_model, int_vars, G);

	if (__found_feasible == 1)
		g_out << "Number of matchings : " << total_matchings << endl;
}


int NUM_THREADS = 1;

bool ReadGraph(int &cnt, int &my_cnt, ListGraph &G){
	bool ok = 1;
	#pragma omp critical 
	{ 
		ok = (bool)(readNautyGraph(G, cin));
		cnt += ok;
		my_cnt = cnt;
	}
	return ok;
}

void PrintLogProgress(int n, int cnt){

	#pragma omp critical
	{
		ofstream log_progress(to_string(n) + "/log_progress");
		log_progress << "Last read graph" << cnt << endl;
		log_progress << "Best IP/Frac: " << __best_IP << " g" << __best_IP_graph_id << " matching " << __best_IP_matching_id << endl;
		log_progress << "Best BDS/Frac: " << __best_BDS << " g" << __best_BDS_graph_id << " matching " << __best_BDS_matching_id << endl;
		log_progress.close();	
	}
}


/*
	This functions receiv nauty's geng output from stdin(may modify this),
	build a LEMON graph and the log files, and calls the function
	that iterates through all matchings.
*/
void RunNautyInput(int start){
	__best_IP = __best_BDS = 1;
	__best_IP_graph_id = __best_IP_matching_id = __best_BDS_graph_id = __best_BDS_matching_id = 1;

	int cnt = 0;

    #pragma omp parallel num_threads(NUM_THREADS) \
    private(G, __found_feasible, __cur_graph_id, g_out, log_out)\
    shared(cnt, __best_BDS_graph_id, __best_BDS_matching_id, __best_IP_graph_id, __best_IP_matching_id, __best_IP, __best_BDS)
	{
		ListGraph G; // Declare global Graph
		int my_cnt;
		while (ReadGraph(cnt, my_cnt, G)){	


			int n = countNodes(G);
			int m = countEdges(G);

			if (cnt < start) continue;


			if (cnt == 1 and start == 0){ // Create folder to log files, create log stream
				// std::experimental::filesystem::create_directory("./" + to_string(n));
				log_out.open(to_string(n) + "/log"); // overwrite existing log
			}
			else if (start == cnt)
				log_out.open(to_string(n) + "/log", ios::app); // open in append mode

			// Next loop makes shure that lemon graph is consistent with the algorithm input
			int nvtx = n - 1, ok = 1;
			for (ListGraph::NodeIt v(G); v != INVALID; ++v)
				if (G.id(v) != nvtx--){
					cout << "Found inconsistency regarding vertex indexing -- Skip graph number " << cnt << endl;
					ok = 0;
					break;
				}

			// Next loop makes shure that lemon graph edge indexing is consistent
			set<int> q;
			for (ListGraph::EdgeIt e(G); e != INVALID; ++e)
				q.insert(G.id(e));

			assert(q.size() == m);
			assert(*q.rbegin() == m - 1);
			/* 
				Since the input data is massive, we will one write a file if
				there is some feasible solution.
			*/
			__found_feasible = 0;
			__cur_graph_id = my_cnt;

			SolveAllMatchings(G);

			if (__found_feasible)
				g_out.close();

			PrintLogProgress(n, cnt);
		}
	}

	log_out.close();
}