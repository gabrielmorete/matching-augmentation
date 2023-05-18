/*
	Author : Gabriel Morete	
*/

#include <iostream>
#include <sstream>
#include <fstream>
#include <vector>
#include <string>
#include <cassert>
#include <iomanip>
#include <cmath>
using namespace std;

// Working with doubles
const double EPS = 1e-4;
int sign(double x) { return (x > EPS) - (x < -EPS); }


class ExtremePoint{
	public:
		vector<double> point;

		ExtremePoint () {}

	    friend istream &operator >> (istream &is, ExtremePoint &e); 
	    friend ostream &operator << (ostream &os, ExtremePoint &e);

	    int getDim () {return point.size();}
	    
	    double operator[] (int i){
	    	return point[i];
	    }

	    void clear(){
	    	point.clear();
	    }

	    void append(double val){
	    	point.push_back(val);	
	    }
};

istream &operator >>(istream &is, ExtremePoint &p){
	string in;

	if (!getline(is, in)) // use istream standard error
		return is;

    p.clear();

    stringstream ss_in(in);

    /*
		the input is now in a stream, going to read 
		each element (point coordinate) separately
    */

    while (ss_in >> in){
    	// in contains a coordinate

    	double val;

    	if (in.find("/") != in.string::npos){ // fractional coordinate
    		int pos = in.find("/");
    		val = stod(in.substr(0, pos)); // splitting the string in two parts
    		val /= stod(in.substr(pos + 1));
    	}
    	else // integer coordinate
    		val = stod(in);

    	p.append(val);
    }

    return is;
}


ostream &operator << (ostream &os, ExtremePoint &p) 
{
    os << "[ ";
    for (int i = 0; i < p.getDim(); i++) {
		os << p[i] << " ";
    }
    os << "]";
    return os;
}

signed main(int argc, char const *argv[]){

	if (argc < 3){
		cerr << "Usage : frac_points_file edge_list" << endl;
		exit(1);
	}

	fstream edge_file(argv[2]);
	if (!edge_file){
		cerr << "Can't open edges file" << endl;
		exit(1);
	}

	vector<pair<int, int>> edges;
	int a, b;
	while (edge_file >> a >> b)
		edges.push_back({a, b});

	int m = edges.size();
	assert(m > 0);

	fstream frac_file(argv[1]);
	if (!frac_file){
		cerr << "Can't open fractial points file" << endl;
		exit(1);
	}
	
	ExtremePoint fx;
	int mxo = 0;
	while (frac_file >> fx){
		assert(fx.getDim() == m + 1);

		vector<double> deg(m + 1, 0); // m >= n

		int co = 0;

		for (int e = 0; e < m; e++){
			deg[ edges[e].first ] += fx[e + 1];
			deg[ edges[e].second ] += fx[e + 1];
		
			if (sign(fx[e + 1] - 0.5) == 0)
				co++;
		}


		for (int v = 0; v <= m; v++){
			if (sign( deg[v] - floor(deg[v]) ) != 0 and co == 9){
				mxo = max(co, mxo);
				cout << v << " " << deg[v] << endl;
				cout << fx << endl;
					for (int e = 0; e < m; e++)
						cout << '\t' << fx[e + 1] << ' ' << edges[e].first << ' ' << edges[e].second << endl;
				//assert(0);
				break;	
			}
		}
	}

	cout << "--" << mxo << endl;
}


