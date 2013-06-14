/* fitspl.c
   Fit a spline to a set of data points.
   */

/************************************************/
/*                                              */
/*  CMATH.  Copyright (c) 1989 Design Software  */
/*                                              */
/************************************************/

#include "cmath.h"
#if (STDLIBH)
#include <stdlib.h>
#endif
#include <stdio.h>
#include <math.h>

#ifndef NULL
#define  NULL  0
#endif


/*-----------------------------------------------------------------*/

/* Global data */

int    ndg;
double *xg, *yg, *wg;
int    nsg;
double *xsg, *ysg, s1g, s2g;
double *b, *c, *d;


/*-----------------------------------------------------------------*/

#if (PROTOTYPE)

int fitspl (int nd, double x[], double y[], double weight[],
            int ns, double xs[], double ys[],
            double *s1, double *s2,
            double *sums, int *fail)

#else

int fitspl (nd, x, y, weight, ns, xs, ys, s1, s2, sums, fail)
int    nd;
double x[], y[], weight[];
int    ns;
double xs[], ys[], *s1, *s2, *sums;
int    *fail;

#endif

/* Purpose ...
   -------
   Fit a cubic spline to a set of weighted data points.  The sum of
   the residuals squared is used as the measure of fit.  The x-values
   of the spline knots are specified by the user while the y-values
   of the knots and the end slopes are optimized using conjgg().

   Input ...
   -----
   nd      : number of data points. These are numbered 0 .. nd-1.
   x[]     : x-coordinates of the data points
             There is no particular order required.
   y[]     : y-coordinates of the data points
   weight[]: user assigned weights for the data
   ns      : number of knots in the spline 0 ... nspl-1
             This includes the end points.
   xs[]    : x-values of the spline knots.  These must be in
             ascending order.
             xs[j-1] < xs[j], j = 1 ... ns-1
   ys[]    : An initial guess for the y-values of the knots.
   s1      : An initial guess for the slope at xs[0]
   s2      : An initial guess for the slope at xs[n-1]
   sums    : The precision to which the minimum should be found.

   Output ...
   ------
   ys[]    : The fitted y-values of the spline knots.
   s1      : The fitted slope at xs[0].
   s2      : The fitted slope at xs[ns-1].
   sums    : The sum of the residuals squared for the final guess.
   fail    : status indicator
             fail = 1 : illegal values for nd, ns
                        nd < 1, ns < 2
             fail = 2 : xs[] are not in ascending order
             fail = 3 : not all of the weights are positive
             fail = 4 : could not allocate workspace.
             fail = 5 : problems with convergence of the function
                        minimizer

   Workspace ...
   ---------
   Three arrays (b, c, d) are allocated ns double elements each,
   with the array yf[] allocated ns+2 double elements.
   The function minimizer conjgg() then allocates a further
   4n double elements.


   This C code written by ...  Peter & Nigel,
   ----------------------      Design Software,
                               42 Gubberley St,
                               Kenmore, 4069,
                               Australia.

   Version ... 1.0, May 1989
   -------     2.0, Sep 1989 change minimizer to conjgg()
               2.1, 12 dec 89, fixed memory allocation/deallocation

   Notes ...
   -----
   (1) Uses the CMATH routines conjgg(), spline() and seval();

*/

{
int    i, nf;
double *yf;
double reqmin;
int    nfe, nje, numres, flag;

/* equate user arrays with global variables */
nsg = ns; xsg = xs; ysg = ys;           /* spline */
ndg = nd; xg = x; yg = y; wg = weight;  /* data points */

/* clear workspace pointers */
b = (double *) NULL;
c = (double *) NULL;
d = (double *) NULL;
yf = (double *) NULL;

/* Check error conditions */
*fail = 0;

if (ns < 2 || nd < 1)  /* do we have knots and data ? */
   {
   *fail = 1;
   goto LeaveFit;
   }

for (i = 1; i < ns; ++i)  /* are the nots in ascending order */
   {
   if (xs[i-1] >= xs[i])
      {
      *fail = 2;
      goto LeaveFit;
      }
   }

for (i = 0; i < nd; ++i)  /* finite and positive weights */
   {
   if (weight[i] <= 0.0)
      {
      *fail = 3;
      goto LeaveFit;
      }
   }

/* try to allocate workspace */

b = (double *) malloc (ns * sizeof(double));
if (b == NULL)
   {
   *fail = 4;
   goto LeaveFit;
   }
c = (double *) malloc (ns * sizeof(double));
if (c == NULL)
   {
   *fail = 4;
   goto LeaveFit;
   }
d = (double *) malloc (ns * sizeof(double));
if (d == NULL)
   {
   *fail = 4;
   goto LeaveFit;
   }

/* Set up the initial guess and step for the minimizer.
   This example fiddles all of the nodes and the end slopes.
   The user may change this part of the code (and the
   relevant section in L2spl() to fiddle or fix any combination
   of knots and/or slopes. */

nf = ns+2;
yf = (double *) malloc (nf * sizeof(double));
if (yf == NULL)
   {
   *fail = 4;
   goto LeaveFit;
   }
/* pack the fiddle vector */
for (i = 0; i < ns; ++i) yf[i] = ys[i];
yf[ns]   = *s1;
yf[ns+1] = *s2;

reqmin = *sums;
if (reqmin < 1.0e-10) reqmin = 1.0e-10;

/* now fit the spline */
conjgg (L2spl, 1, L2deriv, yf, nf, reqmin, sums, &flag, 5*nf,
        100.0, &numres, &nfe, &nje);

/* unpack fiddle vector */
for (i = 0; i < ns; ++i) ys[i] = yf[i];
*s1 = yf[ns];
*s2 = yf[ns+1];

switch (flag)
   {
   case 0  : break;
   case 1  : *fail = 5; break;
   case 2  : *fail = 5; break;
   case 3  : *fail = 4; break;
   case 4  : *fail = 1; break;
   default : *fail = 1;
   }

LeaveFit:

if (yf != NULL) { free(yf); yf = NULL; }
if (d  != NULL) { free(d); d = NULL; }
if (c  != NULL) { free(c); c = NULL; }
if (b  != NULL) { free(b); b = NULL; }

return 0;
}

/*-----------------------------------------------------------------*/

#if (PROTOTYPE)

double L2spl (int nf, double yf[])

#else

double L2spl (nf, yf)
int    nf;            /* number of knots to fiddle */
double yf[];          /* the current guess         */

#endif

{
double L2norm, t;
int    i, last, flag;

/* fiddle ALL of the knots
   The user may wish to change this for some applications */
nsg = nf - 2;
for (i = 0; i < nsg; ++i) ysg[i] = yf[i];
s1g = yf[nsg];
s2g = yf[nsg+1];

/* evaluate the spline coefficients for the current
   guesses for the knots and end slopes */
spline (nsg, 1, 1, s1g, s2g, xsg, ysg, b, c, d, &flag);
if (flag != 0)  exit(1);

/* now sum the square of the residuals */
L2norm = 0.0;
last   = 1;
for (i = 0; i < ndg; ++i)
   {
   t = seval (nsg, xg[i], xsg, ysg, b, c, d, &last);
   t -= yg[i];
   L2norm += (t * t * wg[i]);
   }

return (L2norm);
}

/*-----------------------------------------------------------------*/

#if (PROTOTYPE)

int L2deriv (int nf, double yf[], double dfdyf[])

#else

int L2deriv (nf, yf, dfdyf)
int nf;
double yf[], dfdyf[];

#endif
{
/* do nothing, but avoid the NOT USED warnings */
dfdyf[nf-1] = yf[nf-1];
return (0);
}

/*-----------------------------------------------------------------*/
