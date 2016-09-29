#ifndef AUXFUNCTIONS_H
#define AUXFUNCTIONS_H

#include <string.h>
#define M_PI_4     0.785398163397448309616

typedef int BOOL;
#define TRUE 1
#define FALSE 0

#define MAXDEG			6			// Max degree of the polynom
#define MAXITER			20			// Maximum permitted iterations (for NR)

static double gTol = 1e-8;
//static double ZERO[MAXDEG+1] = { 0 }; 

typedef struct
{
	short deg;			 // Actual degree of polynom (could be less than 6)
	double c[MAXDEG+1];  // Coefficients of polynom of degree 6 (ordered from low exp to high)
} poly;

//__inline double Max(double x, double y)
//{
//	return (x > y) ? x : y;
//}

//__inline double Min(double x, double y)
//{
//	return (x < y) ? x : y;
//}

//__inline void Swap(double* x, double* y)
//{
//	double tmp = *x;
//	*x = *y;
//	*y = tmp;
//}

/* Compare to zero with tolerance*/
__inline BOOL IsZero(double var)
{
	return (fabs(var) < gTol);
}

/* Returns 1 if x is positive
		  -1 if x is negative
		   0 if x is 0 */
__inline short Sign(double x)
{
	return (x > 0) - (x < 0);
}

/* Compares two numbers with tolerance
   Returns: 0: x == y
			1: x > y
		   -1: x < y */
__inline short Compare(double x, double y)
{
	double t = x-y;
	if (t > gTol)
		return 1;	// x > y
	else if (t < -gTol)
		return -1;	// x < y
	
	return 0; // equal
}

/* Normalize a polynom p */
__inline void NormPoly(poly* p)
{
	short i = 0;
	double NormVal = p->c[p->deg];
	p->c[p->deg] = 1.0;
	for (; i < p->deg; i++)
		p->c[i] /= NormVal;
}

/* Find the degree of the polynom (in case it's less than ) */
__inline void FindPolyDeg(poly* p)
{
	short i = MAXDEG;
	for (; i >= 0; --i)
	{
		if (!IsZero(p->c[i]))
		{
			p->deg = i;
			return;
		}
	}
	
	p->deg = -1;	// All zeros
}

/* Finds the coefficient with the maximum absolute value */
//__inline double FindMaxAbsCoeff(poly* p, short start, short end)
//{
//	short i = end;
//	double MaxVal = 0.0;
//	for (; i >= start; --i)
//	{
//		double c = fabs(p->c[i]);
//		if (c > MaxVal)
//			MaxVal = c;
//	}
//
//	return MaxVal;
//}

/* Evaluate the value of polynom p, with degree deg at the parameter t */
__inline double EvalPoly(poly* p, double t)
{
  short i;
  double res = p->c[p->deg];
 
  for (i = p->deg-1; i >= 0; i--)
	res = res*t + p->c[i];
  
  return res;
}

/* Calculates the derivative dp of a polynom p */
__inline void CalcPolyDeriv(poly* p, poly* dp)
{
	short i = 0;
	for (; i < p->deg; ++i)
		dp->c[i] = p->c[i+1]*(i+1);
	dp->deg = p->deg-1;
}

/* Evaluates the value of p(t) and dp(t) at the same time
   res value is the value of p(t). The derivative is returned in dres */
_inline void EvalPolyAndDer(poly* p, double t, double* res, double* dres)
{
  short i;
  *res = p->c[p->deg];
  *dres = 0.0;
  for (i = p->deg-1; i >= 0; i--)
  {
	*dres = (*dres)*t + *res;
	*res = (*res)*t + p->c[i];
  }
}

/* Copy a polynom from src to target */
__inline void CopyPolynom(poly* src, poly* target)
{
	memcpy(target->c, src->c, sizeof(double)*(MAXDEG+1));
	target->deg = src->deg;
}

/* Initialize a polynom p to all zeros */
//__inline void InitPolyZero(poly* p)
//{
//	memcpy(p->c, ZERO, sizeof(double)*(MAXDEG+1));
//}

/* Multiply a polynom p by a scalar */
__inline void PolyScale(poly* p, double scalar)
{
	short i;
	for (i = p->deg; i >= 0; i--)
		p->c[i] *= scalar;
}


/* Helper to find the starting guess. The idea is to find the intersection of the
   line between the first and last point of this interval (the chord).
   If this part of the curve is linear enough, this guess will be very close to the real root. */
__inline double PolyChordIntersection(poly* p , double t1, double t2)
{
	double y1, y2;
	if (IsZero(t2-t1))
		return t1;
	y1 = EvalPoly(p, t1);
	y2 = EvalPoly(p, t2);
    if (IsZero(y2-y1))
		return t1;

	return t1 - y1*(t2-t1)/(y2-y1);
}

/* Faster arc tangent function 
   From "Efficient approximations for the arctangent function", Rajan, S. Sichun Wang Inkol, R. Joyal, A., May 2006
   http://nghiaho.com/?p=997 */
//__inline double FastArcTan(double x)
//{
//    return M_PI_4*x - x*(fabs(x) - 1)*(0.2447 + 0.0663*fabs(x));
//}

/* Function declarations */
int PolyDivision(poly* p, poly* d, poly* r);

#endif