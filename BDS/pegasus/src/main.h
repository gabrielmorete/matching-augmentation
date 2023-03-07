#ifndef MAIN_DEF
#define MAIN_DEF

/*
	This function receives a LP solution and retuns the edges
	of a global minimum cut and its value.
*/
pair<double, vector<Edge> > FindMinCut(double *sol, int n, int m);

/*
	This function returns a optimum integer solution to MAP.
	If no solution is found, it returns a all -1 edge map.
*/
void IntegerSolution(ListGraph::EdgeMap<int> &IntSol);

/*
	This function returns a optimum fractional solution to MAP.
	If no solution is found, it returns a all -1 edge map.
*/
void FractionalSolution(ListGraph::EdgeMap<double> &FracSol);



#endif