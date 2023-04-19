/*
	10july06
	- using populatebyrow.c for LP-7leaf-t1
	29apr06
	- added LOOP to pick & check file of extreme points, each on one line
	- using populatebyrow.c for LP-9leaf-t2
*/
/*
From: Howard Karloff <howard@research.att.com>
To: jcheriya@math.uwaterloo.ca
Subject: rohitsintegralityratiotester.c
*/

// Input:
//        fractional solution x
//        error parameter epsilon
//        multiplicative factor alpha
// Output:
//        NO if (alpha x) does not dominate any convex combination
//        YES if (alpha(1+epsilon) x) dominates a convex combination
//        FAIL if round upper bound reached

// It uses the following oracle: given a positive cost vector c, find
// a minimum cost integral solution (implemented using CPLEX).

// Assumptions:
//
// entries of a fractional solution are between 0 and 1 and those of
// integral soilutions are either 0 or 1
// dimension dim <= 1000


/* Bring in the CPLEX function declarations and the C library 
   header file stdio.h with the following single include. */

#define NUMROWS    11 /* number of rows, incl. unused row 0, in matrix. */
#define NUMCOLS    22 /* number of columns, incl. unused col 0, in matrix. */
#define NUMNZ      73   /* number of nonzeroes in matrix. */
#define TOL 1.0e-5

#include "cplex.h"  

/* Bring in the declarations for the string functions */

#include <stdlib.h>
#include <string.h>
#include <assert.h>
#include <stdio.h>
#include <math.h>
#include <time.h>

// **************** START OF ROHIT'S DECLARATIONS ******************

//JC added vars
int Nfrac, Nfrac_initial; // Nfrac_initial = no. frac. solns in input file
int znum;              // for reading rational i/p, numerator
int zdenom;            // for reading rational i/p, denominator
char zchar;            // for reading rational i/p, ' ' or '/'

double alpha_initial;	// initial multiplicative factor, as given in input
double alpha_update;	// increment to alpha for "NO" case

double sav_cost[1001];  // for saving cost vector (from last NO to next YES)
double sav_costx,sav_costy;

char sav_filename[100];	// vector (=string) for storing name of .lp file written by cplex
//JC added vars - ends

int dim;             // dimension = number of non-tree edges
double epsilon;      // error parameter
double alpha;        // multiplicative factor

double xx[1001];      // fractional solution
double y[1001];      // integral solution
double ysum[1001];   // cummulative sum of integral solutions
double round = 0.0;    // round number

double c[1001];      // cost vector
double width = 1.0;  // width of the program
double INFTY = 10000000.0;// infinity
int gooddim[1001];    // gooddim[i] = (x[i] != 0)

double RBOUND = 1000000.0; // maximum number of rounds allowed

// **************** END OF ROHIT'S DECLARATIONS ******************

/* Include declaration for functions at end of program */

// int oracle (double *cost, double *y);

#ifndef  CPX_PROTOTYPE_MIN

static int 
   populatebyrow     (CPXENVptr env, CPXLPptr lp),
   populatebycolumn  (CPXENVptr env, CPXLPptr lp),
   populatebynonzero (CPXENVptr env, CPXLPptr lp);

static void
   free_and_null     (char **ptr),
   usage             (char *progname);

#else

static int
   populatebyrow     (),
   populatebycolumn  (),
   populatebynonzero ();

static void
   free_and_null     (),
   usage             ();

#endif

#ifndef  CPX_PROTOTYPE_MIN
int
main (int argc, char **argv)
#else
int
main (argc, argv)
int  argc;
char **argv;
#endif
{
/* Declare pointers for the variables and arrays where we will
   store the optimization results including the status, objective value,
   variable values, dual values, row slacks and variable reduced costs. */

int      solstat;
double   objval;
double   *x = NULL;
double   *ipx = NULL;
double   *pi = NULL;
double   *slack = NULL;
double   *ipslack = NULL;
double   *dj = NULL;
double xval;
char     *ctype = NULL;

double z,costx,costy;
int done;
int seed;

CPXENVptr     env = NULL;
CPXLPptr      lp = NULL;
int           status = 0;
int           LPsolved,IPsolved;
double        LPopt,IPopt;
int           i,j;
int           cur_numrows, cur_numcols;
int           loop;
int           integsoln;

int           trials,numtrials;
double        objcoef[NUMCOLS]; 
double        ratio,maxratio;



setvbuf( stdout, 0, _IONBF, 0); 
setvbuf( stderr, 0, _IONBF, 0);  


   /* Initialize the CPLEX environment */

   env = CPXopenCPLEX (&status);

   /* If an error occurs, the status value indicates the reason for
      failure.  A call to CPXgeterrorstring will produce the text of
      the error message.  Note that CPXopenCPLEXdevelop produces no output,
      so the only way to see the cause of the error is to use
      CPXgeterrorstring.  For other CPLEX routines, the errors will
      be seen if the CPX_PARAM_SCRIND indicator is set to CPX_ON.  */

   if ( env == NULL ) {
   char  errmsg[1024];
      printf ("Could not open CPLEX environment.\n");
      CPXgeterrorstring (env, status, errmsg);
      printf ("%s", errmsg);
      goto TERMINATE;
   }

   /* Turn on output to the screen */

   status = CPXsetintparam (env, CPX_PARAM_SCRIND, CPX_ON);
   if ( status ) {
      printf ( "Failure to turn on screen indicator, error %d.\n", status);
      goto TERMINATE;
   }

   /* Create the problem. */

   lp = CPXcreateprob (env, &status, "rohitsintegralityratiotester");

   /* A returned pointer of NULL may mean that not enough memory
      was available or there was some other problem.  In the case of 
      failure, an error message will have been written to the error 
      channel from inside CPLEX.  In this example, the setting of
      the parameter CPX_PARAM_SCRIND causes the error message to
      appear on stdout.  */

   if ( lp == NULL ) {
      printf( "Failed to create LP.\n");
      goto TERMINATE;
   }

   /* Procedure populatebyrow sets up the original problem.
      Later we will repeatedly change the cost vector (and
      nothing else). */

      status = populatebyrow (env, lp);

   if ( status ) {
      printf ( "Failed to populate problem.\n");
      goto TERMINATE;
   }


   status = CPXsetintparam (env, CPX_PARAM_SCRIND, CPX_OFF);
   if ( status ) {
      printf ("Failure to turn off CPX_PARAM_SCRIND.\n");
      goto TERMINATE;
   }



// *************START OF ROHIT'S CODE ****************************
  
  printf("Enter dimension dim: ");
  scanf("%d",&dim);
  printf("dim = %d\n",dim);
  assert(dim<=1000);
  assert(dim+1==NUMCOLS);

  printf("Enter error parameter epsilon: ");
  scanf("%lf",&epsilon);
  printf("epsilon = %lf\n",epsilon);
  assert(epsilon>0);

  printf("Enter multiplicative factor alpha: ");
  scanf("%lf",&alpha);
  printf("alpha = %lf\n",alpha);
  assert(alpha>=1);
// JC for alpha updating
  alpha_initial = alpha;

  printf("Enter increment to alpha: ");
  scanf("%lf",&alpha_update);
  printf("alpha_update = %lf\n",alpha_update);
  assert(alpha_update>0);

  printf("Enter number of frac'l solns Nfrac: ");
  scanf("%d",&Nfrac);
  printf("Nfrac= %d\n",Nfrac);
  assert(Nfrac>=0);
  Nfrac_initial = Nfrac;

while (Nfrac >= 1) {
  Nfrac -= 1;

  //  get fractional solution and initialize the cost vector

//JC - commentout
//printf("Enter frac'l sol (%d dim) x:\n", dim);

  for (j=1; j<=dim; ++j) {

//JC - changing input from float to rationals znum/zdenom
    scanf("%d",&znum);
    assert((znum>=0.0));
    scanf("%c",&zchar);
    assert((zchar==' ')||(zchar=='/'));
    z = znum ;
    if (zchar == '/') { scanf("%d",&zdenom);
                        assert((zdenom>0.0));
                        z = ((double) znum) / ((double) zdenom);
    }
//  printf("znum=%d, zchar|%c|, zdenom=%d, z=%lf\n", znum, zchar, zdenom, z);

    if (z==0) {               // if fract. sol. has 0 for jth edge
      gooddim[j] = 0;         // set its cost to infinity
      c[j] = INFTY;
    } else {
      gooddim[j] = 1;         // otherwise set its cost to 1
      c[j] = 1.0;
      if (width < (1.0/z)) width = 1.0/z;
    }
    xx[j] = z;
  }
  
//JC - commentout
//for (j=1; j<=dim; ++j) {
//  printf("x[%d]= %lf.\n",j,xx[j]);
//}
//printf("\n");




// ************* main action starts here ****************************

StartLoop:

  round = 0.0;
  for (j=1; j<=dim; ++j) ysum[j] = 0.0;

  while (round<RBOUND) {

    // oracle call (attention: Howard)
    // y = optimum integral solution under cost c
    // Now we need to execute oracle(c,y):
    // oracle(c,y); 
    // I am implementing a procedure which returns the binary
    // vector y minimizing c^T y subject to Ay>=b, where
    // A and b are hardcoded into populatebyrow.

// *************END OF ROHIT'S CODE ****************************
 

   // Now we change the coefficients of the objective function,
   // by executing CPXchgcoef with a reference to row -1.
   CPXchgcoef(env,lp,-1,0,0);   
   for (loop = 1; loop <= NUMCOLS-1; loop++) { // @@@@@@Is "NUMCOLS" right?
      CPXchgcoef(env,lp,-1,loop,c[loop]);   
        // "-1" refers to the objective function 
   }

  

   cur_numrows = CPXgetnumrows (env, lp); 
   cur_numcols = CPXgetnumcols (env, lp);

//JC - commentout
//printf("cur_numrows,cur_numcols= %d %d.\n\n",cur_numrows,cur_numcols);


   x = (double *) malloc (cur_numcols * sizeof(double));
   slack = (double *) malloc (cur_numrows * sizeof(double));
   dj = (double *) malloc (cur_numcols * sizeof(double));
   pi = (double *) malloc (cur_numrows * sizeof(double));

   ctype = (char *) malloc (cur_numcols * sizeof(char));

   if ( x     == NULL ||
        slack == NULL ||
        dj    == NULL ||
        pi    == NULL ||
	ctype == NULL ) {
      status = CPXERR_NO_MEMORY;
      printf ( "Could not allocate memory for solution.\n");
      goto TERMINATE;
   }


   //  Now we convert the problem to an IP. 


   status=CPXchgprobtype(env,lp,CPXPROB_MILP); 
     // ("MIP" has been changed to "MILP")

   for (loop=0; loop<=cur_numcols-1; ++loop) ctype[loop]='B';
     // all variables are binary
   
   status=CPXcopyctype(env,lp,ctype);
   if ( status ) {
      fprintf (stderr, "Failed to copy ctype to make IP.\n");
      goto TERMINATE;
   }

   // Optimize the IP and obtain solution. 

   status = CPXmipopt (env, lp);
   if ( status ) {
      fprintf (stderr, "Failed to optimize IP.\n");
      goto TERMINATE;
   }

   // Write the problem to a file:
   // Howard commented out next 5 lines
/*
   status = CPXwriteprob (env, lp, "rohit.lp", NULL);
   if ( status ) {
      printf ( "Failed to write LP to disk.\n");
      goto TERMINATE;
   }
*/


   ipx = (double *) malloc (cur_numcols * sizeof(double));
   ipslack = (double *) malloc (cur_numrows * sizeof(double));

   if ( ipx     == NULL ||
        ipslack == NULL ) {
      status = CPXERR_NO_MEMORY;
      fprintf (stderr, "Could not allocate memory for solution.\n");
      goto TERMINATE;
   }

   status = CPXgetobjval (env, lp, &objval);
   if ( status ) {

/* jc debug - do not stop here 
      fprintf (stderr,"No MIP objective value available.  Exiting...\n");
      goto TERMINATE;
*/
   }

   IPopt=objval;
// @@@@@@@@@@@@@@@@@@@@
//JC - commentout
// printf ("MIP Solution value  = %lf\n\n", IPopt);

   status=CPXgetmipx (env, lp, ipx, 0, cur_numcols-1);
   if (status) {
     fprintf (stderr, "Failed to get optimal integer x. \n");
     goto TERMINATE;
   }

//JC - commentout
//    printf("ipx[0],ipx[1],ipx[2],ipx[3],ipx[4],ipx[5]= %lf %lf %lf %lf %lf %lf.\n\n",
//    ipx[0],ipx[1],ipx[2],ipx[3],ipx[4],ipx[5]);

   status = CPXgetmipslack (env, lp, ipslack, 0, cur_numrows-1);
   if ( status ) {
      fprintf (stderr, "Failed to get optimal slack values.\n");
      goto TERMINATE; 
   }

//JC - commentout
// printf("cur_numrows,cur_numcols= %d %d.\n\n",cur_numrows,cur_numcols);

   // for (i = 0; i < cur_numrows; i++) {
   //    printf ("Row %d:  Slack = %10lf\n", i, ipslack[i]);
   // }

   for (j = 0; j < cur_numcols; j++) {
      y[j]=ipx[j];  // y, like ipx, is a double
//JC - commentout
//    printf ("Column %d:  Value = %10f\n", j, ipx[j]);
//    printf("y[%d] as returned by CPLEX is %lf.\n",j,y[j]);
// @@@@@
   }

//JC - commentout
//printf("Got to start of MORE OF ROHIT CODE.\n\n");

    // **************** MORE OF ROHIT'S CODE **************8
//JC - commentout
//  printf("round = %lf\n",round);
//  for (j=1; j<=dim; ++j) printf("y[%d]=%lf\n",j,y[j]);
//  printf("\n\n");

    costx = costy = 0;
    for (j=1; j<=dim; ++j) {
      costx += c[j]*xx[j];
      costy += c[j]*y[j];
    }
    if (costy > alpha*costx) {
      printf("NO: ((alpha=%lf) x) does not dominate any convex combination.\n", alpha);
      printf("cost IntSol= %lf, cost FracSol= %lf, ratio= @@@%lf@@@\n\n", costy, costx, costy/costx);

/*	update alpha (add specified increment) and continue with same frac soln xx,
	until we find smallest alpha so xx is accepted, then reset
	alpha and print info for largest alpha s.t. xx was rejected --
	we need to save all relevant info for that event
*/
	sav_costy = costy; sav_costx = costx;
	for (j=1; j<=dim; ++j) sav_cost[j] = c[j] ;

	alpha += alpha_update ;

//JC - commentout
//    return 0;
      goto StartLoop ;
    }
    
    for (j=1; j<=dim; ++j) ysum[j] += y[j];
    round += 1;

    for (j=1; j<=dim; ++j) {
      if (gooddim[j])
	c[j] *= ( 1 + epsilon * ( y[j] / (xx[j]*width) ) );
    }

    done = 1;
    for (j=1; j<=dim; ++j)
      if (ysum[j]/round > alpha*(1+epsilon)*xx[j]) done = 0;
    
    if (done) {
      printf("YES: ((alpha=%lf(1+epsilon=%lf)) x) dominates a convex combination; frac.soln#: %d, round= %d.\n\n",
		alpha, epsilon, Nfrac_initial-Nfrac, (int) round);

/*	if alpha has been updated, then we have been re-trying a frac
	soln xx that has been rejected at least once - so print out all
	the saved info from that rejection, reset alpha and then continue
*/
	if (alpha != alpha_initial) {
	printf("BUT: ((alpha=%lf) x) does not dominate any convex combination; frac.soln#: %d\n",
		alpha-alpha_update, Nfrac_initial-Nfrac);
	printf("cost IntSol= %lf, cost FracSol= %lf, ratio= @@@%lf@@@\n\n",
		sav_costy, sav_costx, sav_costy/sav_costx);
	for (j=1; j<=dim; ++j) printf("x[%d]=%lf \t\t cost[%d]=%lf\n",j,xx[j],j,sav_cost[j]);
	printf("\n");

	// Prepare to write out LP for rejecting cost vector
	// Now we change the coefficients of the objective function,
	// by executing CPXchgcoef with a reference to row -1.
	CPXchgcoef(env,lp,-1,0,0);   
	for (loop = 1; loop <= NUMCOLS-1; loop++) {
	   CPXchgcoef(env,lp,-1,loop,sav_cost[loop]);   
	   // "-1" refers to the objective function 
	}

	sprintf( sav_filename, "treeaug_%d.lp", Nfrac_initial-Nfrac);
	printf("calling cplex to write .lp file: %s\n", sav_filename);
	status = CPXwriteprob (env, lp, sav_filename, NULL);
	if ( status ) {
	   printf ( "Failed to write LP to disk.\n");
	   goto TERMINATE;
	}
	printf("\n");

	// reset alpha
	alpha = alpha_initial ;
	} // if (alpha != alpha_initial) ...

//JC - commentout
//    return 0;
      goto EndLoop ;
    } // if (done) ...
  } // while (round<RBOUND) ...

  printf("FAIL: round upper bound reached -- here is the frac.soln:\n");
  for (j=1; j<=dim; ++j) {
    printf("x[%d]= %lf.\n",j,xx[j]);
  }
  
  printf("\n");

EndLoop: ;
//JC - commentout
//return 1;

}
//JC - end outer while loop
  
   // ******************* END OF SECOND BLOCK OF ROHIT **********

   printf("\nCompleted successfully.\n");

TERMINATE:

   /* Free up the solution */


   free_and_null ((char **) &x);
   free_and_null ((char **) &slack);
   free_and_null ((char **) &dj);
   free_and_null ((char **) &pi);

   /* Free up the problem as allocated by CPXcreateprob, if necessary */

   if ( lp != NULL ) {
      status = CPXfreeprob (env, &lp);
      if ( status ) {
         printf ( "CPXfreeprob failed, error code %d.\n", status);
      }
   }

   /* Free up the CPLEX environment, if necessary */

   if ( env != NULL ) {
      status = CPXcloseCPLEX (&env);

      /* Note that CPXcloseCPLEX produces no output,
         so the only way to see the cause of the error is to use
         CPXgeterrorstring.  For other CPLEX routines, the errors will
         be seen if the CPX_PARAM_SCRIND indicator is set to CPX_ON. */

      if ( status ) {
      char  errmsg[1024];
         printf ( "Could not close CPLEX environment.\n");
         CPXgeterrorstring (env, status, errmsg);
         printf ( "%s", errmsg);
      }
   }
     
//JC - commentout
// printf("\n\nmaxratio=%lf.\n",maxratio);
// printf("\nseed=%d.\n",seed);
// printf("\n\n"); 

   return (status);

}  /* END main */


/* This simple routine frees up the pointer *ptr, and sets *ptr to NULL */

#ifndef  CPX_PROTOTYPE_MIN
static void
free_and_null (char **ptr)
#else
static void
free_and_null (ptr)
char  **ptr;
#endif
{
   if ( *ptr != NULL ) {
      free (*ptr);
      *ptr = NULL;
   }
} /* END free_and_null */  

/* To populate by row, we first create the columns,
and then add the rows. */




#ifndef  CPX_PROTOTYPE_MIN
static int
populatebyrow (CPXENVptr env, CPXLPptr lp)
#else
static int
populatebyrow (env, lp)
CPXENVptr  env;
CPXLPptr   lp;
#endif
{
int      status    = 0;
double   obj[NUMCOLS];
double   lb[NUMCOLS];
double   ub[NUMCOLS];
char     *colname[NUMCOLS];
int      rmatbeg[NUMROWS]; 
int      rmatind[NUMNZ];
double   rmatval[NUMNZ];
double   rhs[NUMROWS];
char     sense[NUMROWS];
char     *rowname[NUMROWS];
int i;
int loop;
CPXchgobjsen (env, lp, CPX_MIN);
/* This loop says, for variables 2..21
that the variable is bounded below by 0 and above by 1.0,
and it gives coef. of obj. function: */
   obj[0]=0; lb[0]=0; ub[0]=0;
   obj[1]=0.000000; lb[1]=0; ub[1]=1.0;
   obj[2]=0.000000; lb[2]=0; ub[2]=1.0;
   obj[3]=0.000000; lb[3]=0; ub[3]=1.0;
   obj[4]=0.000000; lb[4]=0; ub[4]=1.0;
   obj[5]=0.000000; lb[5]=0; ub[5]=1.0;
   obj[6]=0.000000; lb[6]=0; ub[6]=1.0;
   obj[7]=0.000000; lb[7]=0; ub[7]=1.0;
   obj[8]=0.000000; lb[8]=0; ub[8]=1.0;
   obj[9]=0.000000; lb[9]=0; ub[9]=1.0;
   obj[10]=0.000000; lb[10]=0; ub[10]=1.0;
   obj[11]=0.000000; lb[11]=0; ub[11]=1.0;
   obj[12]=0.000000; lb[12]=0; ub[12]=1.0;
   obj[13]=0.000000; lb[13]=0; ub[13]=1.0;
   obj[14]=0.000000; lb[14]=0; ub[14]=1.0;
   obj[15]=0.000000; lb[15]=0; ub[15]=1.0;
   obj[16]=0.000000; lb[16]=0; ub[16]=1.0;
   obj[17]=0.000000; lb[17]=0; ub[17]=1.0;
   obj[18]=0.000000; lb[18]=0; ub[18]=1.0;
   obj[19]=0.000000; lb[19]=0; ub[19]=1.0;
   obj[20]=0.000000; lb[20]=0; ub[20]=1.0;
   obj[21]=0.000000; lb[21]=0; ub[21]=1.0;

colname[0]="dumcol";
colname[1]="x1";
colname[2]="x2";
colname[3]="x3";
colname[4]="x4";
colname[5]="x5";
colname[6]="x6";
colname[7]="x7";
colname[8]="x8";
colname[9]="x9";
colname[10]="x10";
colname[11]="x11";
colname[12]="x12";
colname[13]="x13";
colname[14]="x14";
colname[15]="x15";
colname[16]="x16";
colname[17]="x17";
colname[18]="x18";
colname[19]="x19";
colname[20]="x20";
colname[21]="x21";

rowname[0]="dumrow";
rowname[1]="row1";
rowname[2]="row2";
rowname[3]="row3";
rowname[4]="row4";
rowname[5]="row5";
rowname[6]="row6";
rowname[7]="row7";
rowname[8]="row8";
rowname[9]="row9";
rowname[10]="row10";
status = CPXnewcols (env, lp, NUMCOLS, obj, lb, ub, NULL, colname);
if ( status )  goto TERMINATE;

/* Now add the constraints.  */ 
 
rmatbeg[0]=0;
sense[0]='E';
rhs[0]=0;
rmatbeg[1]=1;
rmatind[1]=1;
rmatval[1]=1.000000;
rmatind[2]=2;
rmatval[2]=1.000000;
rmatind[3]=3;
rmatval[3]=1.000000;
rmatind[4]=4;
rmatval[4]=1.000000;
rmatind[5]=5;
rmatval[5]=1.000000;
rmatind[6]=6;
rmatval[6]=1.000000;
sense[1]='G';
rhs[1]=1.000000;
rmatbeg[2]=7;
rmatind[7]=1;
rmatval[7]=1.000000;
rmatind[8]=7;
rmatval[8]=1.000000;
rmatind[9]=8;
rmatval[9]=1.000000;
rmatind[10]=9;
rmatval[10]=1.000000;
rmatind[11]=10;
rmatval[11]=1.000000;
rmatind[12]=11;
rmatval[12]=1.000000;
sense[2]='G';
rhs[2]=1.000000;
rmatbeg[3]=13;
rmatind[13]=2;
rmatval[13]=1.000000;
rmatind[14]=7;
rmatval[14]=1.000000;
rmatind[15]=12;
rmatval[15]=1.000000;
rmatind[16]=13;
rmatval[16]=1.000000;
rmatind[17]=14;
rmatval[17]=1.000000;
rmatind[18]=15;
rmatval[18]=1.000000;
sense[3]='G';
rhs[3]=1.000000;
rmatbeg[4]=19;
rmatind[19]=3;
rmatval[19]=1.000000;
rmatind[20]=8;
rmatval[20]=1.000000;
rmatind[21]=12;
rmatval[21]=1.000000;
rmatind[22]=16;
rmatval[22]=1.000000;
rmatind[23]=17;
rmatval[23]=1.000000;
rmatind[24]=18;
rmatval[24]=1.000000;
sense[4]='G';
rhs[4]=1.000000;
rmatbeg[5]=25;
rmatind[25]=4;
rmatval[25]=1.000000;
rmatind[26]=9;
rmatval[26]=1.000000;
rmatind[27]=13;
rmatval[27]=1.000000;
rmatind[28]=16;
rmatval[28]=1.000000;
rmatind[29]=19;
rmatval[29]=1.000000;
rmatind[30]=20;
rmatval[30]=1.000000;
sense[5]='G';
rhs[5]=1.000000;
rmatbeg[6]=31;
rmatind[31]=5;
rmatval[31]=1.000000;
rmatind[32]=10;
rmatval[32]=1.000000;
rmatind[33]=14;
rmatval[33]=1.000000;
rmatind[34]=17;
rmatval[34]=1.000000;
rmatind[35]=19;
rmatval[35]=1.000000;
rmatind[36]=21;
rmatval[36]=1.000000;
sense[6]='G';
rhs[6]=1.000000;
rmatbeg[7]=37;
rmatind[37]=6;
rmatval[37]=1.000000;
rmatind[38]=11;
rmatval[38]=1.000000;
rmatind[39]=15;
rmatval[39]=1.000000;
rmatind[40]=18;
rmatval[40]=1.000000;
rmatind[41]=20;
rmatval[41]=1.000000;
rmatind[42]=21;
rmatval[42]=1.000000;
sense[7]='G';
rhs[7]=1.000000;
rmatbeg[8]=43;
rmatind[43]=2;
rmatval[43]=1.000000;
rmatind[44]=3;
rmatval[44]=1.000000;
rmatind[45]=4;
rmatval[45]=1.000000;
rmatind[46]=5;
rmatval[46]=1.000000;
rmatind[47]=6;
rmatval[47]=1.000000;
rmatind[48]=7;
rmatval[48]=1.000000;
rmatind[49]=8;
rmatval[49]=1.000000;
rmatind[50]=9;
rmatval[50]=1.000000;
rmatind[51]=10;
rmatval[51]=1.000000;
rmatind[52]=11;
rmatval[52]=1.000000;
sense[8]='G';
rhs[8]=1.000000;
rmatbeg[9]=53;
rmatind[53]=2;
rmatval[53]=1.000000;
rmatind[54]=3;
rmatval[54]=1.000000;
rmatind[55]=7;
rmatval[55]=1.000000;
rmatind[56]=8;
rmatval[56]=1.000000;
rmatind[57]=13;
rmatval[57]=1.000000;
rmatind[58]=14;
rmatval[58]=1.000000;
rmatind[59]=15;
rmatval[59]=1.000000;
rmatind[60]=16;
rmatval[60]=1.000000;
rmatind[61]=17;
rmatval[61]=1.000000;
rmatind[62]=18;
rmatval[62]=1.000000;
sense[9]='G';
rhs[9]=1.000000;
rmatbeg[10]=63;
rmatind[63]=4;
rmatval[63]=1.000000;
rmatind[64]=6;
rmatval[64]=1.000000;
rmatind[65]=9;
rmatval[65]=1.000000;
rmatind[66]=11;
rmatval[66]=1.000000;
rmatind[67]=13;
rmatval[67]=1.000000;
rmatind[68]=15;
rmatval[68]=1.000000;
rmatind[69]=16;
rmatval[69]=1.000000;
rmatind[70]=18;
rmatval[70]=1.000000;
rmatind[71]=19;
rmatval[71]=1.000000;
rmatind[72]=21;
rmatval[72]=1.000000;
sense[10]='G';
rhs[10]=1.000000;

/* End of constraints. */

status = CPXaddrows (env, lp, 0, NUMROWS, NUMNZ, rhs, sense, rmatbeg,
  rmatind, rmatval, NULL, rowname);
if ( status )  goto TERMINATE;

TERMINATE:

return (status);

}  /* END populatebyrow */

/*
WARNING!  THIS PROGRAM TREATS COLUMN 1 DIFFERENTLY
FROM ALL OTHER COLUMNS!  IS THIS WHAT YOU WANTED?

DON'T FORGET TO CHANGE NUMROWS to 11 (=1+number of constraints),
NUMCOLS TO 22 (=1+number of variables)
and NUMNZ TO 73 (=actual number of nonzeroes)
*/
