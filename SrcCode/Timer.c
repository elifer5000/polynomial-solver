
#include "Timer.h"
#include <stdlib.h>

Timer InitTimer()
{
  Timer t;
  QueryPerformanceFrequency(&t);
  return t;
}

double StartTime(Timer t)                
{
  LARGE_INTEGER startCount;
  QueryPerformanceCounter(&startCount);

  return startCount.QuadPart * (1000000.0 / t.QuadPart);
}

double EndTime(Timer t)                   
{
  LARGE_INTEGER endCount;
  QueryPerformanceCounter(&endCount);

  return endCount.QuadPart * (1000000.0 / t.QuadPart);
}

double GetElapsedTimeInSec(double startTime, double endTime)               
{
  return (endTime - startTime)*0.000001;
}

double GetElapsedTimeInMilliSec(double startTime, double endTime)         
{
  return (endTime - startTime)*0.001;
}

double GetElapsedTimeInMicroSec(double startTime, double endTime)         
{
  return endTime - startTime;
}
