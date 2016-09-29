#include <math.h>
#include "AuxFunctions.h"


// Divides p/d. r is the remainder, which is what we are interested in
int PolyDivision(poly* p, poly* d, poly* r) 
{
	// This function can also calculate the quotient q, but it's not needed in this program
	//poly q;
	short i = 0, j;
	double ratio;

	//q.deg = p->deg - d->deg;
	//if (q.deg < 0)
	if (p->deg - d->deg < 0)
		return -1;

	//for (; i <= q.deg; ++i)
	//	q.c[i] = 0.0;
	//InitPolyZero(&q);
	
	CopyPolynom(p, r);

	for (i = p->deg; i >= d->deg; i--)
	{
		 /*q.c[i - d->deg] =*/ ratio = r->c[i] / d->c[d->deg]; 
		 r->c[i] = 0;     
		 for (j = 0; j < d->deg; j++)
			 r->c[i - d->deg + j] -= d->c[j] * ratio;

	}      
	while (! r->c[--(r->deg)]);
	

	return 0;
}