/*
	Library of functions that handle stdio.
*/

#include "lemon.h"
#include "main.h"
// #include "nauty_reader.cpp"

/*
	This function reads the Graph and esge costs from Stdio. Graph is 0-indexed
	The input format will be
		n m       number of nodes, edges
		a_1 b_1 c_1     edge between a_1, b_1 with cost c_1
		...
		a_m b_m c_m 
*/
void ReadStdioInput(ListGraph::EdgeMap<int> &cost, ListGraph &G){
	cout << "CCCCCCCCCCC" << endl;
	int n, m;
	cin>>n>>m;
	
	assert(n >= 3);
	assert(m >= n);

	for (int i = 0; i < n; i++){
		ListGraph::Node v = G.addNode();
		if (G.id(v) != i)
			cout<<"Error : vertex don't match id"<<endl;
		assert(G.id(v) == i);
	}

	for (int i = 0; i < m; i++){
		int a, b, c;
		cin>>a>>b>>c;

		ListGraph::Edge e = G.addEdge(G.nodeFromId(a), G.nodeFromId(b));
		cost[e] = c;
	}

	set<int> q;
	for (ListGraph::EdgeIt e(G); e != INVALID; ++e)
		q.insert(G.id(e));

	assert(q.size() == m);
	assert(*q.rbegin() == m - 1);
	cout << "BBBBBBBBB" << endl;
}


/*
	This function reads the Graph from Stdio. Graph is 0-indexed
	The input format will be
		n m       number of nodes, edges
		a_1 b_1    edge between a_1, b_1
		...
		a_m b_m
*/
void ReadStdioGraph(ListGraph &G){
	int n, m;
	cin>>n>>m;
	
	assert(n >= 3);
	assert(m >= n);

	for (int i = 0; i < n; i++){
		ListGraph::Node v = G.addNode();
		if (G.id(v) != i)
			cout<<"Error : vertex don't match id"<<endl;
		assert(G.id(v) == i);
	}

	for (int i = 0; i < m; i++){
		int a, b, c;
		cin>>a>>b;

		ListGraph::Edge e = G.addEdge(G.nodeFromId(a), G.nodeFromId(b));
	}

	set<int> q;
	for (ListGraph::EdgeIt e(G); e != INVALID; ++e)
		q.insert(G.id(e));

	assert(q.size() == m);
	assert(*q.rbegin() == m - 1);
}


void RunStdioInput(){
	ListGraph G;
	ListGraph::EdgeMap<int> cost(G);
	cout << "ZZZZZZZZZZZZ" << endl;

	if (__all_matchings){
		
		__best_IP = __best_BDS = 1;
		__best_IP_graph_id = __best_IP_matching_id = __best_BDS_graph_id = __best_BDS_matching_id = 1;


		ReadStdioGraph(G);

		assert(biEdgeConnected(G) == 1);

		int n = countNodes(G);
		cout << "Warning: overwriting files in folder " << n << endl;
		std::experimental::filesystem::create_directory("./" + to_string(n));
		ofstream log_out(to_string(countNodes(G)) + "/log"); // clear log file
		log_out.close();

		SolveAllMatchings(G);

		PrintLogProgress(n, 0, 0);
	}
	else {
		cout << "FFFFFFFFFFFFFFFFF" << endl;
		ReadStdioInput(cost, G);

		assert(biEdgeConnected(G) == 1);

		int n = countNodes(G);
		int m = countEdges(G);

		ListGraph::EdgeMap<int> IntSol(G);
		ListGraph::EdgeMap<double> FracSol(G);
		ListGraph::EdgeMap<bool> BDSSol(G);

		GRBModel frac_model(env);
		GRBVar frac_vars[m];
		BuildFractional(frac_model, frac_vars, G);

		GRBModel int_model(env);
		GRBVar int_vars[m];
		BuildIntegral(int_model, int_vars, G);
		MinimumCut cb = MinimumCut(int_vars, n, m, G);
		int_model.setCallback(&cb);

		BDSAlgorithm BDS(G);

		cout << "AHHHHHHHHHHHHHH" << endl;

		SolveMapInstance(cost, FracSol, IntSol, BDSSol, frac_model, frac_vars, int_model, int_vars, G, BDS);
		int cost_Int = 0;
		int cost_BDS = 0;
		double cost_Frac = 0;

		for (ListGraph::EdgeIt e(G); e != INVALID; ++e){
			int u = G.id(G.u(e));
			int v = G.id(G.v(e));

			cost_Int +=  IntSol[e] * cost[e];
			cost_Frac +=  FracSol[e] * cost[e];
			cost_BDS +=  BDSSol[e] * cost[e];

			// cout << u << ' ' << v<< ' ' << FracSol[e] << ' ' << IntSol[e] << ' ' << BDSSol[e] << endl;
		}

		int id;
		cin >> id;

		// create file "g"+cnt
		ofstream g_out("g" + to_string(id));
		
		g_out << countNodes(G) <<' ' << countEdges(G) << endl << endl;
		
		for (ListGraph::EdgeIt e(G); e != INVALID; ++e)
			g_out << G.id(G.u(e)) << ' ' << G.id(G.v(e)) << endl;
		
		g_out << "----------" << endl << endl;
		
		g_out << 0 << ": ";

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

		// cout << "Cost Fractional " << cost_Frac << endl;
		// cout << "Cost Integral " << cost_Int << endl;
		// cout << "Cost BDS " << cost_BDS << endl;
	}
}
