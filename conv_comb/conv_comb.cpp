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
#include "gurobi_c++.h"
using namespace std;

// Working with doubles
const double EPS = 1e-4;
int sign(double x) { return (x > EPS) - (x < -EPS); }


// Gurobi enviroment
GRBEnv env = GRBEnv(true);

/*
	Combination coefficients threshold
	If run with -coef dividend divisor these values are overwritten
*/
double __comb_dividend = 4;
double __comb_divisor = 3;

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


/*
	Model to find the smallest coefficient to each point.

	Note: Rebuilding the model everytime is not efficient, but 
	for the current application is ok.
*/
double ConvexComb(double *sol, ExtremePoint &fx, vector<ExtremePoint> &int_points){
	int n = int_points.size(); // number of points
	int d = int_points[0].getDim(); // dimension

	try{
		GRBModel model(env);
		GRBVar lambda[n];

		model.set(GRB_IntParam_Method, 0); // Forcing Primal Simplex Method

		// one variable for each int point (combination coefficient)
		for (int i = 0; i < n; i++)
			lambda[i] = model.addVar(0.0, 1.0, 0, GRB_CONTINUOUS, "l_" + to_string(i) );

		GRBLinExpr conv;
		for (int i = 0; i < n; i++)
			conv += lambda[i];

		model.addConstr(conv == 1, "conv_comb");

		// Model will thy to optmize the coefficient of the combination
		GRBVar coef = model.addVar(0.0, GRB_INFINITY, 1, GRB_CONTINUOUS, "coef");

		// combination constraint
		for (int j = 0; j < d; j++){
			GRBLinExpr comb = 0;
			
			for (int i = 0; i < n; i++)
				comb += int_points[i][j] * lambda[i];

			model.addConstr( comb <= coef * fx[j], "coord_" + to_string(j)); 
		}

		model.optimize();
		assert(model.get(GRB_IntAttr_SolCount) > 0);


		double *opt_sol = model.get(GRB_DoubleAttr_X, lambda, n);

		for (int i = 0; i < n; i++)
			sol[i] = opt_sol[i];

		delete[] opt_sol;

		return model.get(GRB_DoubleAttr_ObjVal);
	
	} catch(GRBException e) {
		cout << "Error code = " << e.getErrorCode() << endl;
		cout << e.getMessage() << endl;
	} catch(...) {
		cout << "Exception during optimization" << endl;
	}
}


signed main(int argc, char const *argv[]){

	if (argc < 3){
		cerr << "Usage : frac_points_file int_points_file -verbose -coef dividend divisor" << endl;
		exit(1);
	}

	bool verbose_mode = 0;
	if (argc > 3){

		for (int i = 3; i < argc; i++){
			string s = argv[i];

			if (s == "-verbose")
				verbose_mode = 1;
			else if (s == "-coef"){ // overwrite default coefficients
				s = argv[i + 1];
				__comb_dividend = stod(s);

				s = argv[i + 2];
				__comb_divisor = stod(s);

				i += 2;
			}
			else{
				cerr << "Usage : frac_points_file int_points_file -verbose -coef dividend divisor" << endl;
				exit(1);
			}
		}
	}

	fstream int_file(argv[2]);
	if (!int_file){
		cerr << "Can't open integral points file" << endl;
		exit(1);
	}

	vector<ExtremePoint> int_points;
	ExtremePoint p;
	
	while (int_file >> p){
		
		int_points.push_back(p);
	}

	int n = int_points.size();
	assert(n > 0);

	int d = int_points[0].getDim();
	assert(d > 0);

	env.set(GRB_IntParam_OutputFlag, 0);
	env.start();

	fstream frac_file(argv[1]);
	if (!frac_file){
		cerr << "Can't open fractial points file" << endl;
		exit(1);
	}
	
	cout << "\nRunning with coefficient " << __comb_dividend << "/" << __comb_divisor << endl;

	if (verbose_mode)
		cout << "\nFractional points" << endl;

	ExtremePoint fx;
	int cnt = 0;
	double sol[n], worst = 0;

	while (frac_file >> fx){
		assert(fx.getDim() == int_points[0].getDim());

		double coef = ConvexComb(sol, fx, int_points);
		worst = max(worst, coef);

		if (sign( coef - (__comb_dividend/__comb_divisor) ) <= 0){ // Convex comb exists
			if (verbose_mode){
				cout << ' ' << coef << " x " << fx << endl;

				cout <<  setprecision(3);
				for (int i = 0; i < n; i++)
					if (sign(sol[i]) > 0)
						cout << '\t' << sol[i] << " x " << int_points[i] << endl;
				cout << endl;	
			}
		}
		else{
			cout << "Point " << cnt << " can't be written as a convex combination with coefficient at most " << __comb_dividend << "/" << __comb_divisor << endl;
			cout << "Best coef: " << coef << endl;
			cout << fx << endl;
		}

		bool ok = 1;
		for (int i = 1; i < fx.getDim(); i += 4){
			double sum = 0;
			for (int j = i; j < i + 4; j++)
				sum += fx[j];
			if (sign(sum - 2.0) != 0)
				ok = 0;

		}


		if (ok){
			double wrst = 0;
			int id = 0;
			for (int i = 1; i < fx.getDim(); i++){
				int x = fx[i];
				fx.point[i] = max(0.0	, fx[i] - 0.5);


				double aux = ConvexComb(sol, fx, int_points);
				if (aux > wrst){
					wrst = aux;
					id = i;
				}

				fx.point[i] = x;
			}

			cout<<wrst<<' '<<id<<endl;

			cnt++;
		}
	}

	cout << "Worst coef found " << worst << endl;
}


// /*
// 	Build the model to check if the convex combination exists
// */
// void BuildModel(GRBModel &model, GRBVar *lambda, GRBVar *frac_point, vector<ExtremePoint> &int_points){
// 	ExtremePoint p;

// 	int n = int_points.size(); // number of points
// 	int d = int_points[0].getDim(); // dimension

// 	try{

// 		// one variable for each int point (combination coefficient)
// 		for (int i = 0; i < n; i++)
// 			lambda[i] = model.addVar(0.0, 1.0, 0, GRB_CONTINUOUS, "l_" + to_string(i) );

// 		GRBLinExpr conv;
// 		for (int i = 0; i < n; i++)
// 			conv += lambda[i];

// 		model.addConstr(conv == 1, "conv_comb");

// 		// one variable to each fractional coordinate
// 		for (int i = 0; i < d; i++)
// 			frac_point[i] = model.addVar(0.0, 1.0, 0, GRB_CONTINUOUS, "x_" + to_string(i));

// 		// Model will thy to optmize the coefficient of the combination
// 		GRBVar coef = model.addVar(0.0, GRB_INFINITY, 1, GRB_CONTINUOUS, "coef");

// 		// combination constraint
// 		for (int j = 0; j < d; j++){
// 			GRBLinExpr comb = 0;
			
// 			for (int i = 0; i < n; i++)
// 				comb += int_points[i][j] * lambda[i];

// 			model.addQConstr( comb <= comb * frac_point[j], "coord_" + to_string(j)); 
// 		}

// 		int d = fx.getDim();

// 		for (int i = 0; i < d; i++){ // set fractional point value
// 			frac_point[i].set(GRB_DoubleAttr_LB, fx[i]);
// 			frac_point[i].set(GRB_DoubleAttr_UB, fx[i]);
// 		}

// 		model.optimize();
// 		// assert(model.get(GRB_IntAttr_SolCount) > 0);

// 		return model.get(GRB_DoubleAttr_ObjVal);
	
// 	} catch(GRBException e) {
// 		cout << "Error code = " << e.getErrorCode() << endl;
// 		cout << e.getMessage() << endl;
// 	} catch(...) {
// 		cout << "Exception during optimization" << endl;
// 	}
// }	

// /*
// 	checks for a given point fx if the convex combination exists
// */
// double SolveModel(GRBModel &model, GRBVar *lambda, GRBVar *frac_point, ExtremePoint &fx){
	
// 	try{
// 		int d = fx.getDim();

// 		for (int i = 0; i < d; i++){ // set fractional point value
// 			frac_point[i].set(GRB_DoubleAttr_LB, fx[i]);
// 			frac_point[i].set(GRB_DoubleAttr_UB, fx[i]);
// 		}

// 		model.optimize();
// 		// assert(model.get(GRB_IntAttr_SolCount) > 0);

// 		return model.get(GRB_DoubleAttr_ObjVal);
// 	} catch(GRBException e) {
// 		cout << "Error code = " << e.getErrorCode() << endl;
// 		cout << e.getMessage() << endl;
// 	} catch(...) {
// 		cout << "Exception during optimization" << endl;
// 	}
// }