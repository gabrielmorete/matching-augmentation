// Jochen 2006/2007

#include <ilcplex/ilocplex.h>
#include <iostream>
#include <vector>
#include <string>

using namespace std;

ILOSTLBEGIN

/********************* ExtPoint ***************************/

class ExtPoint
{
    vector<double> col;

public:
    ExtPoint () {}

    friend istream &operator >> (istream &is, ExtPoint &e); 
    friend ostream &operator << (ostream &os, ExtPoint &e);

    int getDim () {return col.size();}
    
    double &at (int i) {
	return col[i];
    }


};

istream &operator >> (istream &is, ExtPoint &e) 
{
    string str;
    if (!getline(is,str)) {
	cerr << "Filo IO error in ExtPoint::operator >>.\n";
	exit(3);
    }

    // empty col
    e.col.clear();

    stringstream ss(str);
    
    string t_st;
    double t_dbl;
    double num,den;
    char op;

    while (ss >> t_st) {
	istringstream ss2(t_st);
	if (t_st.find("/")==string::npos) {
	    // integer 
	    t_dbl=atof(t_st.c_str());
	} else {
	    // fractional notation a/b
	    ss2 >> num >> op >> den;
	    t_dbl=num/den;
	}
	
	e.col.push_back(t_dbl);
    }
    return is;
}
    
ostream &operator << (ostream &os, ExtPoint &e) 
{
    os << "[ ";
    for (int i=0;i<e.col.size();i++) {
	os << e.col[i] << " ";
    }
    os << "]";
    return os;
}

/******************************** global methods *****************************/

void popModel (IloModel model, 
	       IloNumVarArray x, 
	       IloRangeArray c, 
	       ifstream &ip, 
	       vector<ExtPoint*> &eps) 
{
    // read all extreme points

    ExtPoint *e;

    do {
	e=new ExtPoint();
	ip >> *e;
	// cout << *e << "\n";
	eps.push_back(e);
    } while (!ip.eof());

    int n=eps.size();
    int m=eps[0]->getDim();

    cout << "n: " << n << "\n";
    cout << "m: " << m << "\n";

    // populate model by column

    IloEnv env = model.getEnv();

    // add constraints

    char na[20];
    for (int i=0;i<m;i++) {
	// add an empty constraint 0 <= * <= infty
	sprintf(na,"cover(%d)",i);
	c.add(IloRange(env, 0, IloInfinity,na)); 
    }

    // one more constraint for convex comb
    c.add(IloRange(env,1,1,"convex"));

    // add ip-vert variables

    char temp[10];
    for (int i=0;i<n;i++) {
	// create a new column for this extreme point
	IloNumColumn col(env);

	for (int j=0; j<m; j++) {
	    col+=c[j](eps[i]->at(j));
	}
	// convex comb constraint
	col+=c[m](1);

	// add a non-negative variable for this column
	x.add(IloNumVar(col,0,IloInfinity));
	sprintf(temp,"x%d",i);
	x[i].setName(temp);
    }
    model.add(c);
}

int main (int argc, char **argv) 
{
    char fname_ip[50], fname_frac[50];
    
    if (argc<3) {
	cerr << "Syntax: convcomb ipsols fracsols.\n";
	abort();
    }
    
    strcpy(fname_ip,argv[1]);
    strcpy(fname_frac,argv[2]);
    ifstream ip (fname_ip);
    ifstream frac (fname_frac);

    if (!ip||!frac) {
	cerr << "Something wrong with the file names....\n";
	abort();
    }

    /* initialize cplex environment and create model */

    IloEnv   env;
    IloNumVarArray x(env);
    IloRangeArray c(env);
    IloObjective obj = IloMaximize(env);
    IloModel model(env, "ConvComb");
    vector<ExtPoint*> eps;

    // add objective function
    model.add(obj);

    // add variables for extreme points and constraints
    popModel(model,x,c,ip,eps);

    // keep number integer extreme points in n
    int n=eps.size();
    if (n==0) {
	cerr << "No integer extreme points!\n";
	exit(3);
    }

    // add variable for fractional point
    
    x.add(IloNumVar(obj(1),0,IloInfinity,IloNumVar::Float,"y"));
    IloCplex cplex(model);
    
    // disable cplex output
    cplex.setOut(env.getNullStream());
    
    // now iterate over fractional extreme points 

    int m=eps[0]->getDim();

    int counter=1;
    double val;
    ExtPoint e;
    while(!frac.eof()) {
	frac >> e; 
	
	if (e.getDim()!=m) {
	    cerr << "Dimension of fractional and integer extreme points do not match!\n";
	    exit(3);
	}
    
	for(int j=0;j<m;j++) {
	    c[j].setCoef(x[n],-e.at(j));
	}

	// write LP for testing purposes
	cplex.exportModel("test.lp");

	// solve model
	cplex.solve();
	val=cplex.getObjValue();
	
	if (counter==1) printf("\n\n");
	if ((counter%1)==0) printf("count = %d \t\t alpha/val = %lf \n", counter, val);
	if (val<.66666) {
	    printf("Counter-example: fractional solution number %d has alpha/val %lf.\n", 
		   counter, val);

	}

	counter++;
    }

    env.end();
}
