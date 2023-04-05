/*
	Nauty Reader

	Read graph6 input in parallel, generate all matchings and calls solvers.
*/

#include "lemon.h"
#include "main.h"

#include <cstdlib>

/*
	gaps for each algorithm.
		- IP gap is >= __IP_dividend/__IP_divisor
		- BDS gap is > __BDS_dividend/__BDS_divisor
*/

const double __IP_dividend = 4;
const double __IP_divisor = 3;
const double __BDS_dividend = 7;
const double __BDS_divisor = 5;

bool __found_feasible;
int __cur_graph_id, __best_IP_graph_id, __best_IP_matching_id, __best_BDS_graph_id, __best_BDS_matching_id;
double __best_IP, __best_BDS;
int __cur_graph_thread[100];

#pragma omp threadprivate(__found_feasible, __cur_graph_id)

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
	ListGraph &G,
	BDSAlgorithm &BDS){


	ListGraph::EdgeMap<bool> BDSSol(G);
	ListGraph::EdgeMap<int> IntSol(G);
	ListGraph::EdgeMap<double> FracSol(G);

	SolveMapInstance(cost, FracSol, IntSol, BDSSol, frac_model, frac_vars, int_model, int_vars, G, BDS);


	if (__cur_graph_id == 978 and matching_id == 696)
		BDS.PrintAndCheck();

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

	ofstream g_out;

	/* 
		Found a feasible example, print to file
	*/
	if (sign(__IP_divisor * cost_Int - __IP_dividend * cost_Frac) >= 0 or sign(__BDS_divisor * cost_BDS - __BDS_dividend * cost_Frac) > 0){
		if (__found_feasible == 0){ // First matching found for this graph
			// create file "g"+cnt
			g_out.open(to_string(countNodes(G)) + "/g" + to_string(__cur_graph_id));
			
			g_out << countNodes(G) <<' ' << countEdges(G) << endl << endl;
			
			for (ListGraph::EdgeIt e(G); e != INVALID; ++e)
				g_out << G.id(G.u(e)) << ' ' << G.id(G.v(e)) << endl;
			
			g_out << "----------" << endl << endl;
		}
		else
			g_out.open(to_string(countNodes(G)) + "/g" + to_string(__cur_graph_id), ios::app);

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

		g_out.close();
		
		#pragma omp critical
		{
			ofstream log_out(to_string(countNodes(G)) + "/log", ios::app); // clear log file
			log_out << "Found feasible example g" << __cur_graph_id << " matching id " << matching_id << endl;
			log_out << "Int/Frc = " << (double) cost_Int/cost_Frac << " BDS/Frc = " << (double) cost_BDS/cost_Frac << endl;
			log_out << endl;
			log_out.close();
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
	ListGraph &G,
	BDSAlgorithm &BDS){

	if (e_id >= m){
		SolveCurrentMatching(total_matchings, cost, frac_model, frac_vars, int_model, int_vars, G, BDS);
		return;
	}

	if (n_matched >= n - 1){ // matching cant increase, prune
		SolveCurrentMatching(total_matchings, cost, frac_model, frac_vars, int_model, int_vars, G, BDS);
		return;
	}

	// Case 1 : won't add edge e_id to the matching
	FindAllMatchings(e_id + 1, n, m, n_matched, total_matchings, matched, cost, frac_model, frac_vars, int_model, int_vars, G, BDS); 

	// Case 2 : if possible, will add e_id to the matching
	ListGraph::Edge e = G.edgeFromId(e_id);
	if ((matched[G.u(e)] == 0) and (matched[G.v(e)] == 0)){ // May add e_id
		
		matched[G.u(e)] = 1;
		matched[G.v(e)] = 1;
		n_matched += 2;
		cost[e] = 0;
		total_matchings++;

		FindAllMatchings(e_id + 1, n, m, n_matched, total_matchings, matched, cost, frac_model, frac_vars, int_model, int_vars, G, BDS);

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

	BDSAlgorithm BDS(G);

	int total_matchings = 1, n_matched = 0;
	FindAllMatchings(0, n, m, n_matched, total_matchings, matched, cost, frac_model, frac_vars, int_model, int_vars, G, BDS);
}


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

void PrintLogProgress(int n, int cnt, int last){
	#pragma omp critical
	{ // If you interrupt the algorithm, may be empty
		ofstream log_progress(to_string(n) + "/log_progress");
		log_progress << "Last read graph " << cnt << endl; // Careful with this, I'm not using mutex
		log_progress << "Smallest unprocessed graph " << last << endl; // Careful with this, I'm not using mutex
		log_progress << "Best IP/Frac: " << __best_IP << " g" << __best_IP_graph_id << " matching " << __best_IP_matching_id << endl;
		log_progress << "Best BDS/Frac: " << __best_BDS << " g" << __best_BDS_graph_id << " matching " << __best_BDS_matching_id << endl;
		log_progress.close();	
	}
}

int ReadLogProgress(int n){
	ifstream log_progress(to_string(n) + "/log_progress");
	int start = 0;
	if (log_progress){
		string s;
		getline(log_progress, s);
		for (int i = 0; i < 4; i++)
			log_progress>>s;
		start = stoi(s);

		for (int i = 0; i < 6; i++){
			log_progress>>s;
			if (i == 2)
				__best_IP = stod(s);
			if (i == 3)
				__best_IP_graph_id = stoi(s.substr(1));
			if (i == 5)
				__best_IP_matching_id = stoi(s);
		}

		for (int i = 0; i < 6; i++){
			log_progress>>s;
			if (i == 2)
				__best_BDS = stod(s);
			if (i == 3)
				__best_BDS_graph_id = stoi(s.substr(1));
			if (i == 5)
				__best_BDS_matching_id = stoi(s);
		}
	}
	else{
		cout << "Can't open log_progress file"<<endl;
		exit(1);
	}

	cout << " -- log_progress data --" << endl;
	cout << " Start graph: g" << start << endl; // Careful with this, I'm not using mutex
	cout << " Best IP/Frac: " << __best_IP << " g" << __best_IP_graph_id << " matching " << __best_IP_matching_id << endl;
	cout << " Best BDS/Frac: " << __best_BDS << " g" << __best_BDS_graph_id << " matching " << __best_BDS_matching_id << endl;

	return start;
}


/*
	This functions receiv nauty's geng output from stdin(may modify this),
	build a LEMON graph and the log files, and calls the function
	that iterates through all matchings.
*/
void RunNautyInput(int start, int n_threads = 1){
	__best_IP = __best_BDS = 1;
	__best_IP_graph_id = __best_IP_matching_id = __best_BDS_graph_id = __best_BDS_matching_id = 1;

	if (start < 0)
			cout << " Running solver with " << "-log_start -threads " << n_threads << endl;
	else
		cout << " Running solver with " << "-start " << start << " -threads " << n_threads << endl;
	
	cout << " IP gap >= " << __IP_dividend << "/" << __IP_divisor << endl;
	cout << " BDS gap > " << __BDS_dividend << "/" << __BDS_divisor << endl;

	int cnt = 0;

	std::system("export GOMP_CPU_AFFINITY=32-65");


    #pragma omp parallel num_threads(n_threads) \
    shared(cnt, __best_BDS_graph_id, __best_BDS_matching_id, __best_IP_graph_id, __best_IP_matching_id, __best_IP, __best_BDS)
	{
		int id = omp_get_thread_num();	
		ListGraph G; // Declare global Graph
		int my_cnt;
		while (ReadGraph(cnt, my_cnt, G)){	

			__cur_graph_thread[id] = my_cnt;

			int n = countNodes(G);
			int m = countEdges(G);

			#pragma omp critical
			{	
				if (start == -1)
					start = ReadLogProgress(n);
				if (start == 0){ // Create folder to log files, create log stream
					std::experimental::filesystem::create_directory("./" + to_string(n));
					ofstream log_out(to_string(countNodes(G)) + "/log"); // clear log file
					log_out.close();
					start = -2; // Invalid option	
				}
			}

			if (cnt < start) continue;

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

			int min_id = __cur_graph_thread[0]; // not critical
			for (int i = 1; i < n_threads; i++)
				min_id = min(min_id, __cur_graph_id);

			PrintLogProgress(n, cnt, min_id);
		}
	}
}