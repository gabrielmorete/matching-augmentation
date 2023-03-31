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


void Print(ListGraph::EdgeMap<int> &cost, ListGraph::EdgeMap<double> &FracSol, ListGraph &G){
	#pragma omp critical
	{
		cout << countNodes(G) << ' ' << countEdges(G) << endl;
		for (ListGraph::EdgeIt e(G); e != INVALID; ++e){
	 		cout << G.id(G.u(e)) << " " << G.id(G.v(e)) << " " << cost[e] << ' ' << FracSol[e] << endl;
		}
	}
}

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
	ListGraph &G){

	int m = countEdges(G);

	ListGraph::EdgeMap<double> FracSol(G);
	FractionalSolution(cost, FracSol, frac_model, frac_vars, G);

	bool is_half_integral = 1,is_integral = 1;
	for (ListGraph::EdgeIt e(G); e != INVALID; ++e){
		if ((sign(FracSol[e]) != 0) and (sign(FracSol[e] - 0.5) != 0) and (sign(FracSol[e] - 1.0) != 0))
			is_half_integral = 0;
		if ((sign(FracSol[e]) != 0) and (sign(FracSol[e] - 1.0) != 0))
			is_integral = 0;
	}

	if (is_integral)
		return;

	ListGraph::NodeMap<double> cut_val(G, 0);	

	if (is_half_integral)
		cout<<"oi"<<endl;

	if (is_half_integral){
		for (ListGraph::EdgeIt e(G); e != INVALID; ++e){
			cut_val[G.u(e)] += FracSol[e];
			cut_val[G.v(e)] += FracSol[e];
		}	

		for (ListGraph::NodeIt v(G); v != INVALID; ++v){
			bool ok = 0;
			if (sign(cut_val[v] - 2.0) != 0){
				ok = 1;
			}

			if (ok){ // test if there is a change
				GRBModel frac_model_2(env);
				GRBVar frac_vars_2[m];
				ListGraph::EdgeMap<double> FracSol_2(G);

				BuildFractional2(frac_model_2, frac_vars_2, G);
				FractionalSolution(cost, FracSol_2, frac_model_2, frac_vars_2, G);


				double obj_1 = 0, obj_2 = 0;

				for (ListGraph::EdgeIt e(G); e != INVALID; ++e){
					obj_1 += FracSol[e];
					obj_2 += FracSol_2[e];
				}
				
				if (sign(FracSol_2[G.edgeFromId(0)]) < 0){
					cout << "v :" << G.id(v) << endl;
					Print(cost, FracSol, G);
					cout << "----" << endl;
					Print(cost, FracSol_2, G);
				}

				assert(0);
			}

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
	ListGraph &G){

	if (e_id >= m){
		SolveCurrentMatching(total_matchings, cost, frac_model, frac_vars, G);
		return;
	}

	if (n_matched >= n - 1){ // matching cant increase, prune
		SolveCurrentMatching(total_matchings, cost, frac_model, frac_vars, G);
		return;
	}

	// Case 1 : won't add edge e_id to the matching
	FindAllMatchings(e_id + 1, n, m, n_matched, total_matchings, matched, cost, frac_model, frac_vars, G); 

	// Case 2 : if possible, will add e_id to the matching
	ListGraph::Edge e = G.edgeFromId(e_id);
	if ((matched[G.u(e)] == 0) and (matched[G.v(e)] == 0)){ // May add e_id
		
		matched[G.u(e)] = 1;
		matched[G.v(e)] = 1;
		n_matched += 2;
		cost[e] = 0;
		total_matchings++;

		FindAllMatchings(e_id + 1, n, m, n_matched, total_matchings, matched, cost, frac_model, frac_vars, G);

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

	ListGraph::NodeMap<bool> matched(G);

	int total_matchings = 1, n_matched = 0;
	FindAllMatchings(0, n, m, n_matched, total_matchings, matched, cost, frac_model, frac_vars, G);
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
	

	int cnt = 0;

	// std::system("export GOMP_CPU_AFFINITY=32-64");


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