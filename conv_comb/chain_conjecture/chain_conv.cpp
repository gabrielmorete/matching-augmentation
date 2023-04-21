/*
	Program receives two lists A, B of points and tests if every 
	element of the list A can be expressed as a convex combination of the 
	elements of the list B.

	Let a \in A, and q be a rational number

		q . x >= \sum_{b in B} l_b b
		\sum_{b in B} l_b = 1
		l_b >= 0, b \in B

	To run the program, flags -verbose -coef are not mandatory.
	Use -verbose for extra information
	Use -coef to overwrite the default values of the coefficients

	Author : Gabriel Morete	
*/

#include <iostream>
#include <sstream>
#include <fstream>
#include <vector>
#include <string>
#include <cassert>
#include <iomanip>
#include <lemon.h>
using namespace std;


signed main(){
	ListGraph G;

	ListGraph::Node v = G.addNode();
	ListGraph::Node u = G.addNode();
	G.addEdge(v, u)

	cout<<countEdges(G)<<endl;

	G.addEdge(v, u)

	cout<<countEdges(G)<<endl;
	

}