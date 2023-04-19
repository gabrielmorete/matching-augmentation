/*
From: Howard Karloff <howard@research.att.com>
To: jcheriya@math.uwaterloo.ca
Subject: populatebyrowgenerator.c
*/

// This program takes an LP represented in matrix format
// and produces the appropriate CPLEX function populatebyrow.

// #include "cplex.h"  
#include <stdlib.h>
#include <string.h>
#include <assert.h>
#include <stdio.h>
#include <math.h>
#include <time.h>

#define NUMCOLS 10001
#define NUMROWS 10001
#define NUMNZ 1000001
// Do I need NUMNZ?
double a[NUMCOLS][NUMROWS];
int rmatbeg[NUMROWS];


main () {
int i,j,m,n,ind;
int sense[NUMROWS]; int objsense;
double rhs[NUMROWS];
double obj[NUMCOLS];

printf("Enter m,n: ");
scanf("%d%d",&m,&n);
printf("m,n= %d %d.\n",m,n);
assert(n<=NUMCOLS-1); // I don't use column 0.
assert(m<=NUMROWS-1); // I don't use row 0.

// First read in the n coefficients of the objective function:
printf("Enter the %d coefficients of the objective function:\n",n);
for (j=1; j<=n; ++j) {
  scanf("%lf",&obj[j]);
}
obj[0]=0; // unused

printf("Enter -1 if goal is MIN and +1 if goal is MAX: ");
scanf("%d",&objsense);
assert((objsense==-1)||(objsense==1));

for (i=1; i<=m; ++i) {
  printf("Enter the %d coefficients of row %d:\n",n,i);
  for (j=1; j<=n; ++j) {
    scanf("%lf",&a[i][j]);
  }
  printf("Enter -1 for <=, 0 for =, 1 for >=: \n");
  scanf("%d",&sense[i]);
  assert((sense[i]==-1)||(sense[i]==0)||(sense[i]==1));
  printf("Enter rhs of constraint %d: ",i);
  scanf("%lf",&rhs[i]);
  printf("Row %d is: ",i);
  for (j=1; j<=n; ++j) {
    printf("%lf ",a[i][j]);
  }
  if (sense[i]==-1) printf("<=");
  if (sense[i]== 0) printf("=");
  if (sense[i]==+1) printf(">=");
  printf("%lf.\n", rhs[i]);
} 


printf("/* To populate by row, we first create the columns,\n");
printf("and then add the rows. */\n");
printf("\n");
printf("#ifndef  CPX_PROTOTYPE_MIN\n");
printf("static int\n");
printf("populatebyrow (CPXENVptr env, CPXLPptr lp)\n");
printf("#else\n");
printf("static int\n");
printf("populatebyrow (env, lp)\n");
printf("CPXENVptr  env;\n");
printf("CPXLPptr   lp;\n");
printf("#endif\n");
printf("{\n");
printf("int      status    = 0;\n");
printf("double   obj[NUMCOLS];\n");
printf("double   lb[NUMCOLS];\n");
printf("double   ub[NUMCOLS];\n");
printf("char     *colname[NUMCOLS];\n");
printf("int      rmatbeg[NUMROWS]; \n"); 
printf("int      rmatind[NUMNZ];\n");
printf("double   rmatval[NUMNZ];\n");
printf("double   rhs[NUMROWS];\n");
printf("char     sense[NUMROWS];\n");
printf("char     *rowname[NUMROWS];\n");
printf("int i;\n");
printf("int loop;\n");
if (objsense==-1)
  printf("CPXchgobjsen (env, lp, CPX_MIN);\n"); 
if (objsense==1)
  printf("CPXchgobjsen (env, lp, CPX_MAX);\n"); 

/*jc 10jul06  - changed 2 to 1*/
printf("/* This loop says, for variables 1..%d\n",n);
printf("that the variable is bounded below by 0 and above by 1.0,\n"); 
printf("and it gives coef. of obj. function: */\n");

printf("   obj[0]=0; lb[0]=0; ub[0]=0;\n"); // dummy column
for (j=1; j<= n; j++) {
printf("   obj[%d]=%lf; lb[%d]=0; ub[%d]=1.0;\n",j,obj[j],j,j);
}

printf("\n");
printf("colname[0]=\"dumcol\";\n");

for (j=1; j<=n; ++j) {
  printf("colname[%d]=\"x%d\";\n",j,j);
}

printf("\n");
printf("rowname[0]=\"dumrow\";\n");

for (i=1; i<=m; ++i) {
  printf("rowname[%d]=\"row%d\";\n",i,i);
}

printf("status = CPXnewcols (env, lp, NUMCOLS, obj, lb, ub, NULL, colname);\n");
printf("if ( status )  goto TERMINATE;\n");
printf("\n");
printf("/* Now add the constraints.  */ \n");
printf(" \n");

ind=0;
rmatbeg[0]=0;
printf("rmatbeg[0]=0;\n");
printf("sense[0]='E';\n");
printf("rhs[0]=0;\n");

for (i=1; i<=m; ++i) {
  rmatbeg[i]=ind+1;
  printf("rmatbeg[%d]=%d;\n",i,rmatbeg[i]);
  for (j=1; j<=n; ++j) {
    if (a[i][j]!=0) {
      ++ind;
      printf("rmatind[%d]=%d;\n",ind,j);
      printf("rmatval[%d]=%lf;\n",ind,a[i][j]);
    }
  }
  // Now add the sense of the constraint and its rhs.
  if (sense[i]==-1) printf("sense[%d]='L';\n",i);
  if (sense[i]== 0) printf("sense[%d]='E';\n",i);
  if (sense[i]==+1) printf("sense[%d]='G';\n",i);
  printf("rhs[%d]=%lf;\n",i,rhs[i]);
}


printf("\n/* End of constraints. */\n");
printf("\n");
printf("status = CPXaddrows (env, lp, 0, NUMROWS, NUMNZ, rhs, sense, rmatbeg,\n");
printf("  rmatind, rmatval, NULL, rowname);\n");
printf("if ( status )  goto TERMINATE;\n");
printf("\n");
printf("TERMINATE:\n");
printf("\n");
printf("return (status);\n");
printf("\n");
printf("}  /* END populatebyrow */\n");

printf("\n\nWARNING!  THIS PROGRAM TREATS COLUMN 1 DIFFERENTLY\n");
printf("FROM ALL OTHER COLUMNS!  IS THIS WHAT YOU WANTED?\n\n");

printf("DON'T FORGET TO CHANGE NUMROWS to %d (=1+number of constraints),\n",1+m);
printf("NUMCOLS TO %d (=1+number of variables)\n",1+n);
printf("and NUMNZ TO %d (=actual number of nonzeroes)\n\n",ind+1);


}
