

#include <cpx_wrap.h>




int main (int argc, char **argv) {

    cpx_wrap lpsolve;

    char fname_ip[50], fname_frac[50], temp[500];
    

    if (argc<3) {
		cerr << "Syntax: convcomb ipsols fracsols.\n";
		abort();
    }
    
    strcpy(fname_ip,argv[1]); // nome do arquivo com pontos inteiros
    strcpy(fname_frac,argv[2]); // nome do arquivo com pontos fracionarios

    ifstream ip (fname_ip); // arquivo de leitura
    ifstream frac (fname_frac);

    if (!ip||!frac) {
		cerr << "Something wromg with the file names....\n";
		abort();
    }

    /* find dimension of extreme points */

    ip.getline(temp, 500);
    char *pos = temp;
    int m = 0;
    int i = 0;
    while(pos[i]) { // Finding the number of elements of the line, dimention of the space
		while((pos[i] == ' ' ) && (pos[i])) {
		    i++;
		}
		if (pos[i]) {
		    m++;
		    while((pos[i]!=' ')&&(pos[i])) {
				i++;
		    }
		}
    }

    // add two extra rows for convcomb
    m += 2;

    /* count number of extreme points */

    int n = 1;

    while(!ip.eof()) {
		ip.getline(temp,500);
		n++;
    }
    n--;

    // have one column more than needed for data
    

    cout << "m: " << m << "\n";
    cout << "n: " << n << "\n";


    /* create ipvert matrix */

    spm_t A(m,n); 

    ip.close();
    ifstream ip2 (fname_ip);
    int j=0;
    int z;
    int num,den;

    do {
        ip2.getline(temp,500);

		if (!ip2.eof()) {
		    istringstream istr (temp);
		    i = 0;

		    cout << j << ": ";
		    
		    while (!istr.eof()) {
				istr >> num; 
			
				if (!istr.eof()) {
				    cout << num << " ";
				    A(i++,j) = -num;
				}
		    }
		    cout << "\n";
		    j++;
		}

    } while (!ip2.eof());
	
    // conv comb rows

    for (int j=0; j < n - 1; j++) {
		A(m-2,j) = 1;
		A(m-1,j) = -1;
    }

    cout << A << "\n";

    /* rhs */

    vec_t b(m);
    for (int i=0; i<m-2; i++) {
	b(i)=0;
    }

    b(m-2)=1;
    b(m-1)=-1;

    /* define c */

    vec_t c(n);
    
    for (int i=0; i<n; i++ ) {
	c(i)=0;
    }
    c(n-1)=1;

    /* now iterate over fractional extreme points */

    double f;
    char st[500];

    while (!frac.eof()) {
	frac.getline(temp,500);    
	
	if (!frac.eof()) {
	    istringstream istr (temp);
	    i=0;

	    while (!istr.eof()) {
		istr >> st;

		if (!istr.eof()) {
		    pos=strchr(st,'/');
		    
		    if (pos>0) {
			sscanf(st,"%d/%d",&num,&den);
			A(i,n-1)=(double)num/den;
		    } else {
			sscanf(st,"%d",&num);
			A(i,n-1)=(double)num;
		    }
		    i++;
		}
	    }

	    /* call CPLEX */
	    
	    IloNumArray x;
	    IloNumArray y;
	    IloNum val;
	    IloAlgorithm::Status status;

	    status=lpsolve.solve(A,b,c,x,y,val);

	    cout << "--------------------------------------------------------------------\n";
	    cout << "CPLEX terminated with status: " << status << "\n";
	    cout << "alpha: " << val << "\n";

	    if (val<.66666) {
		cout << "Counter-example found!\n";
		abort();
	    }

//	    cout << "x: " << x << "\n";
//	    cout << "y: " << y << "\n";
	}
    }
}
