#include <stdio.h>
#include <malloc.h>
#include <assert.h>
#define _USE_MATH_DEFINES 
#include <math.h>
#include "AuxFunctions.h"

/* Definitions */
//#define M_PI_3		1.04719755119659774615		// PI/3
#define M_2PI			6.28318530717958647692		// 2*PI
#define SIDE		0.25		// length of triangle side
#define A1			1.0			// distance from A2 to origin (A1)
#define A2			1.0			// distance from A3 to origin (A1)
#define COS_PI_3	0.5			// cos(pi/3)
#define SIN_PI_3	0.8660254	// sin(pi/3)
//#define TOL_FAST_ATAN		1e-3/SIDE		// Tolerance up to which fast atan can give good results

/* Global variables */
static double gEps = 0.0; // User desired epsilon
static int numSolsFirst = 0; // Number of solutions for first polynom

/* Function declarations */
BOOL NewtonRaphson(poly* p, double boundLo, double boundHi, double* x);

/* Returns the number of sign changes in polynom p, (Descartes rule) */
short FindSignChanges(poly* p)
{
	short i = p->deg-1, Signs = 0;
	double tnew;
	double t = Sign(p->c[p->deg]);
	for (; i >= 0; --i)
	{
		if (IsZero(p->c[i]))
			continue;	// Ignore zero coefficients
		tnew = Sign(p->c[i]);
		if (tnew != t)	// Compare with previous coeff
			Signs++;
		t = tnew;
	}
	return Signs;
}

/* Builds the f(t) polynom, from which the angle will be found, according to the actuators lengths.
   Builds both the theta, and (theta - pi) polynoms (see explanation in report) */
int Build2Polynoms(double L1, double L2, double L3, poly* Poly1, poly* Poly2)
{
	double a, b, c, d, e, f, g;
	// Precalculations
	double L1_2 = L1*L1; // L1^2
	double L2_2 = L2*L2; // L2^2
	double L3_2 = L3*L3; // L3^2
	double L1_4 = L1_2*L1_2; // L1^4
	double L2_4 = L2_2*L2_2; // L2^4
	double L3_4 = L3_2*L3_2; // L3^4
	// Calculate the coefficients (ordered from low to high)
	a = 2000.72+5650.38*L1_4+2578.38*L2_4-3332.68*L3_2+2304.*L3_4+L1_2*(1059.45-5924.76*L2_2-5376.*L3_2)+L2_2*(-3384.12+768.*L3_2);
	b = -2120.03-1499.24*L1_4-2048.*L2_4+1995.32*L3_2+L1_2*(475.9+3547.24*L2_2-548.76*L3_2)+L2_2*(2342.56+548.76*L3_2);
	c = 14714.4+22546.4*L1_4+11282.4*L2_4-18761.4*L3_2+11008.*L3_4+L1_2*(1201.4-22820.8*L2_2-22272.*L3_2)+L2_2*(-18401.3+256.*L3_2);
	d = -10313.-2998.48*L1_4-4096.*L2_4+7537.89*L3_2+L1_2*(2451.04+7094.48*L2_2-1097.52*L3_2)+L2_2*(7281.89+1097.52*L3_2);
	e = 34895.7+28141.6*L1_4+14829.6*L2_4-32571.3*L3_2+15104.*L3_4+L1_2*(-3372.31-27867.2*L2_2-28416.*L3_2)+L2_2*(-32245.5-1792.*L3_2);
	f = -12289.-1499.24*L1_4-2048.*L2_4+5542.56*L3_2+L1_2*(1975.14+3547.24*L2_2-548.76*L3_2)+L2_2*(4939.32+548.76*L3_2);
	g = 26278.1+11245.6*L1_4+6125.62*L2_4-17142.6*L3_2+6400.*L3_4+L1_2*(-3514.26-10971.2*L2_2-11520.*L3_2)+L2_2*(-17228.3-1280.*L3_2);

	// Build 2 polynoms
	Poly1->c[0] = a; Poly1->c[1] = b; Poly1->c[2] = c; Poly1->c[3] = d; Poly1->c[4] = e; Poly1->c[5] = f; Poly1->c[6] = g;
	Poly2->c[0] = g; Poly2->c[1] = -f; Poly2->c[2] = e; Poly2->c[3] = -d; Poly2->c[4] = c; Poly2->c[5] = -b; Poly2->c[6] = a;
	FindPolyDeg(Poly1);
	FindPolyDeg(Poly2);
	if (Poly1->deg <= 0 || Poly2->deg <= 0)
		return -1;	// Error - got a 0 polynom

	// Normalize polynoms
	NormPoly(Poly1);
	NormPoly(Poly2);
	
	return 0;	// All ok
}

/* Finds the Sturm sequence for polynom p
   The delta number of sign changes give the number of roots in the interval.
   Input: boundLo and boundHi - the range to search for roots.
   Output: the number of roots in the given range. */
int Sturm(poly* p, double boundLo, double boundHi)
{
	short i = 2, Signs;
	poly SLow, SHigh;
	poly Seq[MAXDEG+1];
	CopyPolynom(p, &Seq[0]);
	CalcPolyDeriv(p, &Seq[1]);
	for (; i <= MAXDEG; ++i)
	{
		if (PolyDivision(&Seq[i-2], &Seq[i-1], &Seq[i]) != 0)
			return 0;
		PolyScale(&Seq[i], -1);
	}
	for (i = 0; i <= MAXDEG; ++i)
	{
		SLow.c[i] = EvalPoly(&Seq[i], boundLo);
		SHigh.c[i] = EvalPoly(&Seq[i], boundHi);
	}
	FindPolyDeg(&SLow);
	FindPolyDeg(&SHigh);
	Signs = FindSignChanges(&SLow) - FindSignChanges(&SHigh);

	return Signs;
}

/* Find the intervals where the roots are.
   This is a recursive function that divides each interval in two parts,
   and keeps subdividing the sections where roots are found, until reaching
   a size of Eps.
   When an interval where only one root is found is reached, Newton-Raphson method is used to more
   rapidly converge to a solution. If NR goes out of bounds, control is returned to this function to
   refine the search.
   Input: polynom p. boundLo and boundHi - the range to search for roots.
   Output: the roots and num of roots are updated through their pointers. */
void FindRootsInterval(poly* p, double boundLo, double boundHi, double* roots, int* numRoots)
{
	double split, delta = boundHi - boundLo;
	int RootsLeft, RootsRight, div = 2;
	BOOL ContinueFind = TRUE;

	if (delta < gEps)
	{
		roots[(*numRoots)++] = (boundHi + boundLo)/2.;
		return;
	}
	while (IsZero(EvalPoly(p, split = boundLo + delta/div)))	// root at division limit?
	{
		div++; // Try another division point (1/3, 1/4, etc...)
		if (div > 10)
		{	// Avoid too small division. Seems this is a small interval anyways, so return the middle point
			roots[(*numRoots)++] = (boundHi + boundLo)/2.;
			return;
		}
	}
	
	RootsLeft = Sturm(p, boundLo, split);	// Find number of roots on left interval
	RootsRight = Sturm(p, split, boundHi);	// Find number of roots on right interval
	if (RootsLeft > 0)
	{
		if (RootsLeft == 1)
		{	// Try NR
			double x = (boundLo + split)/2.; // Starting point
			//double x = PolyChordIntersection(p, boundLo, split); // Alternative
			if (NewtonRaphson(p, boundLo, split, &x) == TRUE)
			{
				ContinueFind = FALSE;
				roots[(*numRoots)++] = x;
			} // if NR fails will continue with FindRootsInterval
		}
		if (ContinueFind)
			FindRootsInterval(p, boundLo, split, roots, numRoots);
	}
	if (RootsRight > 0)
	{
		ContinueFind = TRUE;
		if (RootsRight == 1)
		{	// Try NR
			double x = (split + boundHi)/2.; // Starting point
			//double x = PolyChordIntersection(p, split, boundHi); // Alternative
			if (NewtonRaphson(p, split, boundHi, &x) == TRUE)
			{
				ContinueFind = FALSE;
				roots[(*numRoots)++] = x;
			}	// if NR fails will continue with FindRootsInterval
		}
		if (ContinueFind)
			FindRootsInterval(p, split, boundHi, roots, numRoots);
	}
}

/*	Finds a root inside an interval using Newton-Raphson's method. If it falls outside the bounds,
	or the max amount of iterations is exceeded it will return FALSE.
	Input: polynom p. boundLo and boundHi - the range to search for roots.
	Output: TRUE if succeded. The root location x. */
BOOL NewtonRaphson(poly* p, double boundLo, double boundHi, double* x)
{
	int i = 0;
	double f, derf, dx, X = *x;
	while (i++ < MAXITER)
	{
		EvalPolyAndDer(p, X, &f, &derf);	// Find the value and derivateve of p at X
		dx = f/derf;
		X -= dx;
		if (X > boundHi || X < boundLo) // out of bounds - go back to FindRootsInterval
			return FALSE;
		if (fabs(dx) < gEps || fabs(f) < gEps) // if change has become small
		{
			//printf("NR iterations %d\n", i);
			*x = X;
			return TRUE;	// Root found
		}		
	}

	return FALSE; // Exceeded number of iterations without finding
}

/* Once the angle has been found, get the XY coords
   For details on the calculations done here see report.
   Input: Actuator lengths (L1, L2, L3)
		  t - a root of the polynom found
		  isPoly2 - it's a solution of the second polynom (must compensate by PI)
		  *angles - angles found so far
		  idx - index of current angle in *angles array
   Output: xi,yi - the coordinates of the triangle for this solution.
		   Returns -1 if no solution can be found, 0 otherwise. */
int GetXYCoords(double L1, double L2, double L3, double t, BOOL isPoly2, double* angles, int idx,
				double *x1, double *y1, double *x2, double *y2, double *x3, double *y3)
{
	int i = 0;
	double a, b, a2, b2, a3, b3, c1, c2, den;
	double ct, st; // cos(theta) and sin(theta). can get away with calculating only one trig function which is expensive
	double theta;
	//if (gEps > TOL_FAST_ATAN) // For large tolerances this approximation seems to be good
	//	theta = 2*FastArcTan(t);
	//else
		theta = 2*atan(t);
	
	if (isPoly2)
		theta += M_PI;
	theta = fmod(theta + M_PI, M_2PI) - M_PI; // restrict x so that -M_PI < x < M_PI

	angles[idx] = theta; // Store angle to check for duplicate later
	if (isPoly2 && numSolsFirst > 0) // Check there's no duplicate solution from poly1 (could happen at -PI/2 or PI/2)
	{
		for (; i < idx; ++i)
			if (Compare(angles[i], theta) == 0) // Duplicate angle (could happen at -PI/2 or PI/2)
				return -1; // Duplicate
	}
	ct = cos(theta);
	//sqrt faster than doing sin(theta)
	if (theta >= 0) // 0 to PI sinus is positive
		st = sqrt(1-ct*ct);
	else
		st = -sqrt(1-ct*ct);

	// B2 related
	a = SIDE*ct;
	a2 = a - A1;
	b2 = SIDE*st;
	// B3 related
	a3 = SIDE*(COS_PI_3*ct - SIN_PI_3*st);	// SIDE*cos(M_PI_3 + theta)
	b = SIDE*(SIN_PI_3*ct + COS_PI_3*st);	// SIDE*sin(M_PI_3 + theta);
	b3 = b - A2;
	
	den	= 2*(a2*b3 - a3*b2);
	if (IsZero(den))
		return -1;	// Invalid solution - division by zero

	c1 = L2*L2-L1*L1-a2*a2-b2*b2;
	c2 = L3*L3-L1*L1-a3*a3-b3*b3;
	// Calculate points
	*x1 = (b3*c1 - b2*c2)/den;
	*y1 = (-a3*c1 + a2*c2)/den;
	*x2 = *x1 + a;
	*y2 = *y1 + b2;
	*x3 = *x1 + a3;
	*y3 = *y1 + b;	    

	return 0;
}

/* Main solver function 
   Input: Actuator lengths (L1, L2, L3)
		  eps - desired tolerance
   Output: xi,yi - the coordinates of the triangle for this solution.
		   numSols - number of solutions found		
		   Returns -1 if no solution can be found, 0 otherwise. */
int GetForwardKinematics(double eps, double L1, double L2, double L3,
                      double *x1, double *y1, double *x2, double *y2, double *x3, double *y3, int *numSols)
{	
	int i = 0, j = 0, n;
	double t[MAXDEG]; // At most can have 6 solutions
	double theta[MAXDEG];
	poly ThePoly, ThePoly2;

	gEps = eps/(2*SIDE);	// Init global epsilon var
	*numSols = 0; // init num of sols to 0

	if (Build2Polynoms(L1, L2, L3, &ThePoly, &ThePoly2) != 0)
		return -1; // Got a 0 polynom
	// Find roots for the first polynom
	FindRootsInterval(&ThePoly, -1, 1, t, numSols);
	numSolsFirst = *numSols;
	// Find roots for the second polynom
	FindRootsInterval(&ThePoly2, -1, 1, t, numSols);
	n = *numSols;
	// Find the final coordiantes for each solution found
	for (i = 0; i < n; ++i, ++j)
	{
		if (GetXYCoords(L1, L2, L3, t[i], (i < numSolsFirst) ? FALSE : TRUE, theta, i,
						&x1[j], &y1[j], &x2[j], &y2[j], &x3[j], &y3[j]) != 0)
		{
			--(*numSols); // Invalid or duplicate solution, remove it
			--j; // Decrease index of solutions, so that in next 'for' it will stay in the same place for xi,yi
		}
	}

	if (*numSols == 0)
		return -1;	// No solution found

	return 0; /* Everything is OK.*/
}

