/*
	Program receivs two lists A, B of extreme points and tests if every 
	element of the list A can be expressed as a convex combination of the 
	elements of the list B.

	Let a \in A, and q be a rational number

		q . x >= \sum_{b in B} l_b b
		\sum_{b in B} l_b = 1
		l_b >= 0, b \in B
*/

#include <iostream>
#include <sstream>
#include <fstream>
#include <vector>
#include <string>
#include <cassert>
#include <iomanip>
#include "gurobi_c++.h"

// Utilities
#define dbg(x)  cout << #x << " = " << x << endl

// Safe handling doubles
const double EPS = 1e-4;
/*
	-1 if x < -EPS
	0 if  -EPS <= x <= EPS
	1 if EPS < x
*/ 
int sign(double x) { return (x > EPS) - (x < -EPS); } 
///





/*
	Combination coefficient
*/
const double __comb_dividend = 4;
const double __comb_divisor = 3;


using namespace std;

class ExtremePoint{
	private:
		vector<double> point;

	public:
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

	if (!getline(is, in))
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

    	if (in.find("/") != in.string::npos){ // fractional point
    		int pos = in.find("/");
    		val = stod(in.substr(0, pos)); // spliting the string in two parts
    		val /= stod(in.substr(pos + 1));
    	}
    	else{
    		val = stod(in);
    	}

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



void BuildModel(GRBModel &model, GRBVar *lambda, GRBVar *frac_point, vector<ExtremePoint> &int_points){
	ExtremePoint p;

	int n = int_points.size(); // number of points
	int d = int_points[0].getDim(); // dimension

	// one variable to each int point
	for (int i = 0; i < n; i++)
		lambda[i] = model.addVar(0.0, 1.0, 0, GRB_CONTINUOUS, "l_" + to_string(i) );

	GRBLinExpr conv;
	for (int i = 0; i < n; i++)
		conv += lambda[i];

	model.addConstr(conv == 1, "conv_comb");

	// one variable to each fractional coordinate
	for (int i = 0; i < d; i++)
		frac_point[i] = model.addVar(0.0, 1.0, 0, GRB_CONTINUOUS, "x_" + to_string(i));


	for (int j = 0; j < d; j++){
		GRBLinExpr comb = 0;
		
		for (int i = 0; i < n; i++)
			comb += int_points[i][j] * lambda[i];

		model.addConstr( comb <= (__comb_dividend/__comb_divisor) *  frac_point[j], "coord_" + to_string(j)); 
	}
}	

int SolveModel(GRBModel &model,  GRBVar *lambda, GRBVar *frac_point, ExtremePoint &fx){
	int d = fx.getDim();

	for (int i = 0; i < d; i++){
		frac_point[i].set(GRB_DoubleAttr_LB, fx[i]);
		frac_point[i].set(GRB_DoubleAttr_UB, fx[i]);
	}

	model.optimize();

	if (model.get(GRB_IntAttr_SolCount) > 0)
		return 1;
	return 0;
}


signed main(int argc, char const *argv[]){

	if (argc < 3){
		cerr << "Usage : frac_points_file int_points_file -verbose" << endl;
		exit(1);
	}

	bool verbose_mode = 0;
	if (argc == 4){
		string vrb = argv[3];
		if (vrb != "-verbose"){
			cerr << "Usage : frac_points_file int_points_file -verbose" << endl;
			exit(1);
		}
		verbose_mode = 1;	
	}

	fstream int_file(argv[2]);
	if (!int_file){
		cerr << "Cant open integral points file" << endl;
		exit(1);
	}


	vector<ExtremePoint> int_points;
	ExtremePoint p;
	
	if (verbose_mode)
		cout << "Integer points" << endl;

	while (int_file >> p){
		if (verbose_mode)
			cout << p << endl;
		
		int_points.push_back(p);
	}

	int n = int_points.size();
	assert(n > 0);

	int d = int_points[0].getDim();
	assert(d > 0);

	GRBEnv env = GRBEnv(true);
	env.set(GRB_IntParam_OutputFlag, 0);
	env.start();

	GRBModel model(env);
	GRBVar lambda[n];
	GRBVar frac_point[d];


	BuildModel(model, lambda, frac_point, int_points);

	fstream frac_file(argv[1]);
	if (!frac_file){
		cerr << "Cant open fractial points file" << endl;
		exit(1);
	}

	
	cout << "Running with coefficient " << __comb_dividend << "/" << __comb_divisor << endl;

	if (verbose_mode)
		cout << "\nFractional points" << endl;



	ExtremePoint fx;
	int cnt = 0;
	while (frac_file >> fx){
		if (SolveModel(model, lambda, frac_point, fx)){
			if (verbose_mode){
				cout << fx << "   ";

				double *sol = model.get(GRB_DoubleAttr_X, lambda, n);

				cout <<  setprecision(3);
				for (int i = 0; i < n; i++)
					cout << sol[i] << ' ';
				cout << endl;

				delete[] sol;
			}
		}
		else{
			cout << "Point " << cnt << " cant be written as a convex combination" << endl;
			cout << fx << endl; 
		}

		cnt++;
	}
}
