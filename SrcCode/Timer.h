#ifndef TIMER_H
#define TIMER_H

/*
Timer functions (implemented using Windows Timer API functions). 
Usage example:

{
  Timer t;
  double startTime=0.0, endTime=0.0;
  
  t = InitTimer();
  startTime = StartTime(t);
  ...
  endTime = EndTime(t);
  printf("Elapsed milli seconds: %lf\n", GetElapsedTimeInMilliSec(startTime, endTime));
}

*/

#include <windows.h>

typedef LARGE_INTEGER Timer;

Timer InitTimer();
double StartTime(Timer t);                 /* Starting time in micro-second. */
double EndTime(Timer t);                   /* Ending time in micro-second. */
double GetElapsedTimeInSec(double startTime, double endTime);        /* Get elapsed time in second (same as getElapsedTime). */
double GetElapsedTimeInMilliSec(double startTime, double endTime);   /* Get elapsed time in milli-second. */
double GetElapsedTimeInMicroSec(double startTime, double endTime);   /* Get elapsed time in micro-second. */

#endif
