/*  conjgg.c
    Multidimensional minimization using the conjugate gradient
    method.
    (Requires well behaved function function derivatives.)  */

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

/* zeps  = a small number which protects against trying to achieve
           fractional accuracy for a minimum that happens to be
           exactly zero.  */
#define   zeps      EPSILON
#define   tiny      1.0e-20
#define   zero      0.0

/* golden ratio */
#define   gold      1.618034
#define   cgold     0.3819660
#define   glimit    100.0

#define   SIGN(a,b)  ( ((b) >= zero) ? fabs(a) : -fabs(a) )
#define   MAX(a,b)   ( ((a) > (b)) ? (a) : (b) )

/*-----------------------------------------------------------------*/

#if (PROTOTYPE)

int conjgg (double (*f)(int n, double x[]), int mf,
            int (*df)(int n, double x[], double dfdx[]),
            double p[],
            int n, double ftol, double *fret,
            int *flag, int itmax, double xmax,
	    int *iter, int *nfe, int *nje)

#else

int conjgg (f, mf, df, p, n,
            ftol, fret, flag,
	    itmax, xmax,
	    iter, nfe, nje)

double (*f) ();
int    mf, (*df) ();
double p[];
int    n;
double ftol, *fret;
int    *flag, itmax;
double xmax;
int    *iter, *nfe, *nje;

#endif

/* Purpose ...
   -------
   Given a starting point p[], Fletcher-Reeves-Polak-Ribiere
   minimization is performed on a function f(), using its gradient
   as calculated by a function df().

   Input...
   -----
   f     : user defined objective function
           double f (n, x)
           int n;
           double x[];
           {
           ...
           return (double value);
           }
   mf    : method flag for partial derivative calculation
           mf = 0 , user function df() will be used
           mf = 1 , finite differences will be used.  In this case
                    the user supplied function df() need not do
                    anything.
   df    : function for evaluating derivatives
           int df (n, x, dfdx)
           int n;
           double x[], dfdx[];
           {
           ...
           return (0);
           }
   p     : starting point in n-dimensional space
   n     : number of elements in p
   ftol  : convergence tolerance on the function value
   itmax : maximum allowed number of iterations
   xmax  : bounds on parameter values for 1D minimization

   Output...
   ------
   fret  : minimum of f
   iter  : number of iterations performed
   nfe   : number of function evaluations
   nje   : number of derivative evaluations
   flag  : = 0, normal return
           = 1, did not converge within itmax iterations
           = 2, could not bracket a minimum in a line minimization
           = 3, could not allocate workspace
           = 4, invalid user input, n < 1, ftol <= 0.0, p == NULL,
                itmax < 1, xmax <= 0.0.

   Workspace...
   ---------
   uvect    : vector of dimension n
   g, h, xi : vectors of dimension n

   This C code written by ...  Peter & Nigel,
   ----------------------      Design Software,
                               42 Gubberley St,
                               Kenmore, 4069,
                               Australia.

   Version... 1.0, October 1988
   -------    2.0, June    1989   full function prototypes
                                  workspace allocation
                                  n in user function calls
              2.1, August   1989  input checking
                                  derivative estimation

   Notes ...
   -----
   (1) Uses routines linem() and braket() to perform the one-
       dimensional minimizations.
   (2) Adapted from the FORTRAN code FRPRMN in
       W.H. Press et al
       Numerical Recipes. The art of scientific computing.
---------------------------------------------------------------------*/

{  /* begin conjgg() ... */
int    j, Lflag, its;
double fp, gg, dgg, gam, xmin;
double *uvect, *g, *h, *xi;

*flag = 0;
if (n < 1 || ftol <= zero || p == NULL || itmax < 1 ||
    xmax <= 0.0)
    {
    *flag = 4;
    return (0);
    }

uvect = NULL;
g = NULL;
h = NULL;
xi = NULL;

uvect = (double *) malloc(n * sizeof(double));
if (uvect == NULL)
   {
   *flag = 3;
   goto Finish;
   }
g = (double *) malloc(n * sizeof(double));
if (g == NULL)
   {
   *flag = 3;
   goto Finish;
   }
h = (double *) malloc(n * sizeof(double));
if (h == NULL)
   {
   *flag = 3;
   goto Finish;
   }
xi = (double *) malloc(n * sizeof(double));
if (xi == NULL)
   {
   *flag = 3;
   goto Finish;
   }

*nfe = 0;
*nje = 0;
fp = (*f) (n, p);
++(*nfe);
if (mf == 0)
   (*df) (n, p, xi);
else
   partiald (n, f, p, fp, xi, nfe);
++(*nje);

for (j = 0; j < n; ++j)
   {
   g[j] = -xi[j];
   h[j] = g[j];
   xi[j] = h[j];
   }

/* Iterate !  */

for (its = 1; its <= itmax; ++its)
   {
   *iter = its;
   /* Search along direction xi  */
   *fret = linem (f, p, xi, n, ftol, 100, &xmin, xmax, uvect,
                  nfe, &Lflag);
   if (Lflag != 0)
      { /* one dimensional search failed */
      *flag = 2;
      goto Finish;
      }

   /* Check convergence. */
   if ( (fabs((*fret)-fp)) <=
        (ftol * (1.0 + 0.5 * (fabs(*fret) + fabs(fp))) + zeps) )
      {
      /* Normal return */
      *flag = 0;
      goto Finish;
      }

   fp = (*f) (n, p);
   ++(*nfe);
   if (mf == 0)
      (*df) (n, p, xi);
   else
      partiald (n, f, p, fp, xi, nfe);
   ++(*nje);
   gg = zero;
   dgg = zero;

   for (j = 0; j < n; ++j)
      {
      gg += g[j]*g[j];
      /* The following statement for Fletcher-Reeves */
      /* dgg += xi[j]*xi[j]; */
      /* The following statement for Polak-Ribiere */
      dgg += (xi[j] + g[j]) * xi[j];
      }

   if (gg == zero)
      {
      /* Unlikely, but if the gradients are zero then we are done. */
      *flag = 0;
      return (0);
      }

   /* Determine a new direction */
   gam = dgg / gg;
   for (j = 0; j < n; ++j)
      {
      g[j] = -xi[j];
      h[j] = g[j] + gam * h[j];
      xi[j] = h[j];
      }

   }

/* Too many iterations without convergence. */
*flag = 1;

Finish:
/* Clean up work space */
if (xi != NULL) { free(xi); xi = NULL; }
if (h  != NULL) { free(h); h = NULL; }
if (g  != NULL) { free(g); g = NULL; }
if (uvect != NULL) { free(uvect); uvect = NULL; }
return (0);

}  /* end of conjgg() */

/*-----------------------------------------------------------------*/

/* Bracket a minimum along a line.  */

#if (PROTOTYPE)

int braket (double (*f)(int n, double x[]),
            double pvect[], double direct[],
            int n, double *ax, double *bx, double *cx,
            double bound, double *fa, double *fb, double *fc,
            double uvect[], int *nfe, int *flag)

#else

int braket (f, pvect, direct, n,
            ax, bx, cx, bound, fa, fb, fc,
            uvect, nfe, flag)

double  (*f) ();
double  pvect[], direct[];
int     n;
double  *ax, *bx, *cx, bound, *fa, *fb, *fc;
double  uvect[];
int     *nfe, *flag;

#endif

/* Purpose ...
   -------
   Given a function F, a point in N-dimensional space,
   a direction to search and given distinct initial parameter
   values AX and BX, this routine searches in the downhill
   direction (defined by the function evaluated at the initial points)
   and returns new parameter values AX, BX, CX which bracket
   a minimum of the function.

   Input ...
   -----
   f      : user defined objective function that returns
            a double precision value for each n-dimensional point
   pvect  : origin for line along which to search
   direct : direction vector for search
   n      : number of elements in pvect
   ax     : guess for left bracketing parameter
   bx     : guess for right bracketing parameter
   bound  : limit on magnitude of ax, bx, cx (say 1000.0)
   nfe    : number of function evaluations so far

   Output ...
   ------
   ax, bx, cx : values of parameter bracketing a minimum
                such that fc < fb < fa and cx lies between
                ax and bx
   fa, fb, fc : values of the objective function at ax, bx and cx
   nfe        : number of function evaluations
   flag       : = 0, normal return
		= 1, could not bracket within bounds

   Workspace ...
   ---------
   uvect  : n-dimensional points corresponding to parameter u
            where uvect[j] = pvect[j] + u * direct[j]

   This C code written by ...  Peter & Nigel,
   ----------------------      Design Software,
                               42 Gubberley St,
                               Kenmore, 4069,
                               Australia.

   Version ... 1.0, October 1988.
   -------

   Notes ...
   -----
   (1) gold = default ratio by which successive intervals are
              magnified
   (2) glimit = maximum magnification allowed by the parabolic-fit
                step
   (3) Adapted from the FORTRAN code MNBRAK in
       W.H. Press et al
       Numerical Recipes. The art of scientific computing.
----------------------------------------------------------------------*/

{  /* begin braket() */

/* local variables */
int    j;
double axL, bxL, cxL, faL, fbL, fcL;
double u, fu, temp, r, q, ulim;

*flag = 0;
axL = *ax;
bxL = *bx;

for (j = 0; j < n; ++j)  uvect[j] = pvect[j] + axL * direct[j];
faL = (*f) (n, uvect);
for (j = 0; j < n; ++j)  uvect[j] = pvect[j] + bxL * direct[j];
fbL = (*f) (n, uvect);
(*nfe) += 2;

if (fbL > faL)
   {
   /* Switch roles of a and b so that we go downhill in the
      direction from a to b  */
   temp = bxL;  bxL = axL;  axL = temp;
   temp = fbL;  fbL = faL;  faL = temp;
   }

/*  First guess for C  */
cxL = bxL + gold * (bxL - axL);
for (j = 0; j < n; ++j)  uvect[j] = pvect[j] + cxL * direct[j];
fcL = (*f) (n, uvect);
++(*nfe);

/* keep returning here until we bracket */
while (fbL >= fcL)
   {
   /* Compute U by parabolic extrapolation until we bracket.  */
   r = (bxL - axL) * (fbL - fcL);
   q = (bxL - cxL) * (fbL - faL);
   /* Tiny is used to prevent possible division by zero. */
   u = bxL - ((bxL - cxL) * q - (bxL - axL) * r) /
             (2.0 * SIGN(MAX(fabs(q-r),tiny),q-r));
   /*
   temp = fabs(q-r);
   if (tiny > temp) temp = tiny;
   if ((q-r) < zero) temp = -temp;
   u = bxL - ((bxL - cxL) * q - (bxL - axL) * r) / (2.0 * temp);
   */

   for (j = 0; j < n; ++j)  uvect[j] = pvect[j] + u * direct[j];
   /* We won't go farther than ulim.  */
   ulim = bxL + glimit * (cxL - bxL);

   /* Now test various possibilities... */
   if ((bxL - u) * (u-cxL) > zero)
      {
      /* Parabolic U is between B and C, try it */
      fu = (*f) (n, uvect);
      if (fu < fcL)
         { /* Got a minimum between B and C */
         axL = bxL;  faL = fbL;
         bxL = u;    fbL = fu;
         /* exit from this step (iteration) */
         continue;
	 }
      else if (fu > fbL)
	 { /* Got a minimum between A and U  */
         cxL = u;  fcL = fu;
         /* exit from this step */
         continue;
         }
      /* Parabolic fit was no use. Use default magnification. *
      u = cxL + gold * (cxL - bxL);
      for (j = 0; j < n; ++j)  uvect[j] = pvect[j] + u * direct[j];
      fu = (*f) (uvect);
      ++(*nfe);
      }

   elseif ((cxL - u) * (u - ulim) > zero)
      {  Parabolic fit is between C and its allowed limit. */
      fu = (*f) (n, uvect);
      ++(*nfe);
      if (fu < fcL)
	 {
         bxL = cxL;  fbL = fcL;
         cxL = u;    fcL = fu;
         u = cxL + gold * (cxL - bxL);
         for (j = 0; j < n; ++j) uvect[j] = pvect[j] + u * direct[j];
         fu = (*f) (n, uvect);
         ++(*nfe);
         }
      }

   else if ((u - ulim) * (ulim - cxL) >= zero)
      { /* Limit parabolic U to its maximum allowed value. */
      u = ulim;
      for (j = 0; j < n; ++j)  uvect[j] = pvect[j] + u * direct[j];
      fu = (*f) (n, uvect);
      ++(*nfe);
      }

   else
      { /* Reject parabolic U, use default magnification. */
      u = cxL + gold * (cxL - bxL);
      for (j = 0; j < n; ++j)  uvect[j] = pvect[j] + u * direct[j];
      fu = (*f) (n, uvect);
      ++(*nfe);
      }

   /* Eliminate oldest point and continue.  */
   axL = bxL;  faL = fbL;
   bxL = cxL;  fbL = fcL;
   cxL = u;    fcL = fu;

   /* Check limit on parameter values  */
   if (fabs(u) > bound)
      { /* We are out of bounds without bracketing  */
      *flag = 1;
      goto Finish;
      }

   /*  Take another step  */
   }

Finish:
*ax = axL;
*bx = bxL;
*cx = cxL;
*fa = faL;
*fb = fbL;
*fc = fcL;
return(0);

}   /* end of braket() */

/*-----------------------------------------------------------------*/

/* One dimensional function minimizer along a specified line.  */

#if (PROTOTYPE)

double linem (double (*f)(int n, double x[]),
              double pvect[], double direct[],
              int n, double tol, int itmax, double *xmin,
              double bound, double uvect[], int *nfe, int *flag)

#else

double linem (f, pvect, direct, n,
              tol, itmax, xmin, bound,
              uvect, nfe, flag)

double (*f)();
double pvect[], direct[];
int    n;
double tol;
int    itmax;
double *xmin, bound;
double uvect[];
int    *nfe, *flag;

#endif

/* Purpose ...
   -------
   Given a function F, a starting point PVECT and a direction
   to search DIRECT, this routine first brackets and then
   isolates the minimum to a fractional precision of
   about TOL using Brent's method.

   Input ...
   -----
   f        : externally defined objective function f(x)
   pvect    : origin of line along which to search
   direct   : vector direction of search
   n        : number of elements in pvect
   tol      : precision to which the minimum should be found
              For well behaved functions the value of TOL should
              be set greater than the square root of the machine
              precision.
   itmax    : number of iterations allowed
              There is one function evaluation per iteration.
              For a well behaved function, 100 should be plenty.
   bound    : limit on the magnitude of the distance moved along
              the search line (say 1000.0)
   nfe      : number of function evaluations made before entry

   Output ...
   ------
   xmin     : parameter value at minimum
   pvect    : vector "abscissa" at minimum
   direct   : initial direction scaled by xmin
   linem    : minimum value of objective function
   nfe      : total number of function evaluations
   flag     : = 0, normal return
              = 1, exceeded maximum number of iterations
              = 2, out of bounds without bracketing, no valid
                   result is returned

   WorkSpace ...
   ---------
   uvect    : vector of dimension n defining a point
              It has elements
              uvect[j] = pvect[j] + u * direct[j]
              where u is a parameter measuring along
              the line "DIRECT".

   This C code written by ...  Peter & Nigel,
   ----------------------      Design Software,
                               42 Gubberley St,
                               Kenmore, 4069,
                               Australia.

   Version ... 1.0  October, 1988
   -------     2.0  July,    1989  fixed parabolic step selection
                                   to be more reliable
                                   n in user function calls

   Notes ...
   -----
   (1) Adapted from the text
       W.H. Press et al
       Numerical Recipes. The art of scientific computing.
   (2) Does not require derivative information.
-------------------------------------------------------------------*/

{  /* begin linem()  */

int    iter, j, bflag;
double ax, bx, cx, fa, fb, fc;
double a, b, v, w, x, e, fx, fv, fw;
double xm, tol1, tol2,r, q, p, etemp, d;
double u, fu;

*flag = 0;
d = zero;

/*  Bracket the minimum  */

/*  First, guess the bracket  */
ax = 1.0;
bx = 2.0;
/*  then, improve it  */
braket (f, pvect, direct, n, &ax, &bx, &cx, bound,
        &fa, &fb, &fc, uvect, nfe, &bflag);
if (bflag != 0)
   { /* Could not bracket minimum in bounds */
   *flag = 2;
   return (zero);
   }

/* a and b must be in ascending order, though the abscissas
   need not be.  */
if (ax < cx)
   { a = ax;  b = cx; }
else
   { a = cx;  b = ax; }

/* Initializations ... */

x = bx;   /* x == best guess for the minimum  */
w = x;    /* w == next highest function value */
v = x;    /* v == largest function value      */

/* e will be the distance moved on the step before the last */
e = zero;

for (j = 0; j < n; ++j) uvect[j] = pvect[j] + x * direct[j];
fx = (*f) (n, uvect);
++(*nfe);
fv = fx;
fw = fx;

/* Main Loop ... */

for (iter = 1; iter <= itmax; ++iter)
   {
   xm = 0.5 * (a + b);
   tol1 = tol * fabs(x) + zeps;
   tol2 = 2.0 * tol1;

   /* Test done here. */
   if (fabs(x-xm) <= (tol2-0.5*(b-a))) goto Finish;

   if (fabs(e) > tol1)
      {  /* The last step was ok,
            Construct a trial parabolic fit
            for this step.  */
      r = (x - w) * (fx - fv);
      q = (x - v) * (fx - fw);
      p = (x - v) * q - (x - w) * r;
      q = 2.0 * (q - r);
      if (q > zero) p = -p;
      q = fabs(q);
      etemp = e;
      e = d;

      if ((fabs(p) >= fabs(0.5 * q * etemp)) || (p <= q * (a-x))
          || (p >= q * (b-x)) )
         {  /* The papabolic fit is no good;
               Take a golden step into the larger
               of the two segments.  */
         if(x >= xm) e = a - x; else  e = b - x;
         d = cgold * e;
         }
      else
         {  /* The parabolic step is ok, use it.  */
         d = p / q;
         u = x + d;
         if ((u-a < tol2) || (b-u < tol2))  d = SIGN(tol1, (xm-x));
         }
      }
   else
      {  /* The parabolic fit is no good;
            Take a golden step into the larger
            of the two segments.  */
      if(x >= xm) e = a - x; else  e = b - x;
      d = cgold * e;
      }


   /* Arrive here with d computed either from the parabolic fit
      or the golden section  */

   if (fabs(d) >= tol1)
      u = x + d;
   else  /* move at least a little. */
      u = x + SIGN(tol1, d);

   /* This is the one function evaluation per iteration. */
   for (j = 0; j < n; ++j)  uvect[j] = pvect[j] + u * direct[j];
   fu = (*f)(n, uvect);
   ++(*nfe);

   /* Now we have to decide what to do with the function value. */
   if (fu <= fx)
      {  /* Good, we have moved downhill ... */
      if (u >= x) a = x; else b = x;  /* update the bracket */
      v = w;  fv = fw;
      w = x;  fw = fx;
      x = u;  fx = fu;   /* the new minimum */
      }
   else
      {  /* the new point is uphill ... */
      if (u < x) a = u; else b = u;    /* update the bracket */
      if ((fu <= fw) || fabs(w - x) < zeps)
	 {  /* the new point is not as low as x but
               it is better than both v and w  */
         v = w;  fv = fw;
         w = u;  fw = fu;
	 }
      else if ((fu <= fv) || fabs(v - x) < zeps || fabs(v - w) < zeps)
	 {  /* the new point is better than v only */
         v = u;  fv = fu;
	 }
      /* if the new point is no good, forget it */
      }

   }   /* go back for another iteration. */

/* We have exceeded the maximum iteration count;
   Return the best guess for the minimum even
   if we did not achieve the desired tolerance.  */
*flag = 1;

Finish:
/* get out after recording the best estimate for the minimum. */
*xmin = x;
for (j = 0; j < n; ++j)
   {
   direct[j] *= x;
   pvect[j] += direct[j];
   }
return (fx);

}  /* end of linem()  */

/*-----------------------------------------------------------------*/

#if (PROTOTYPE)

int partiald (int n, double (*f)(int n, double x[]),
              double x[], double fp, double dfdx[], int *nfe)

#else

int partiald (n, f, x, fp, dfdx, nfe)
int    n;
double (*f)(), x[], fp, dfdx[];
int    *nfe;

#endif

/* Purpose ...
   -------
   Evaluate the partial derivatives using finite differences.

   Input ...
   -----
   n     : number of elements in the independent variable array
   f     : user supplied function (see conjgg())
   x     : the current position
   fp    : the current function value
   nfe   : current function call count

   Output ...
   ------
   dfdx  : the partial derivatives
           dfdx[j] = d f(x) / dx[j],  j = 0 ... n-1
   nfe   : new function call count

   Version ... 1.0  August 1989
   -------

   This C code written by ...  Peter & Nigel,
   ----------------------      Design Software,
                               42 Gubberley St,
                               Kenmore, 4069,
                               Australia.

*/

{  /* begin partiald() */
double step, xold, fdelx;
int    j;

/* in selecting the step size for the finite differences, we assume
   that the function variables are reasonably well scaled */
step = sqrt(EPSILON) * (1.0 + fabs(fp));

for (j = 0; j < n; ++j)
   {
   xold = x[j];
   x[j] += step;
   fdelx = (*f) (n, x);
   dfdx[j] = (fdelx - fp) / step;
   x[j] = xold;
   }
(*nfe) += n;
return (0);
}

/*-----------------------------------------------------------------*/

