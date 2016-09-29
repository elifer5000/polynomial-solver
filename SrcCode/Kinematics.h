#ifndef KINEMATICS_H
#define KINEMATICS_H

/* Returns 0 if all is ok, something else if error occured.*/
int GetForwardKinematics(double eps, double L1, double L2, double L3,
                      double *x1, double *y1, double *x2, double *y2, double *x3, double *y3, int *numSols);

#endif

