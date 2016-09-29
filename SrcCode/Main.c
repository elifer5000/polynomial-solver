#include "Timer.h"
#include "PostScriptOutput.h"
#include "Kinematics.h"

#include <assert.h>
#include <math.h>

#pragma warning(disable : 4996) /* This line disables warning 4996, which we are aware of and do intentionally.*/

#define NUM_RUNS 1000   // To get a better estimate of the running time should use the average of many runs

void DrawConfigurationInFile(char* FileName, 
                             double x1, double x2, double x3,
                             double y1, double y2, double y3)
{
	double PPolygonX[3], PPolygonY[3];
	double SegX[2];
	double SegY[2];
	double RGB[3] = {0.0,0.0,1.0};
	FILE* OutputFile = OpenPostScriptOutput(FileName);

	PPolygonX[0] = x1; PPolygonX[1] = x2; PPolygonX[2] = x3;
	PPolygonY[0] = y1; PPolygonY[1] = y2; PPolygonY[2] = y3;

	/*Draw X axis*/
	SegX[0] = 0; SegX[1] = 1;
	SegY[0] = 0; SegY[1] = 0;
	DrawPostScriptPolygon(OutputFile, SegX, SegY, 2,
			0, 0, 0.5,
			RGB);
	/*Draw Y axis*/
	SegX[0] = 0; SegX[1] = 0;
	SegY[0] = 0; SegY[1] = 1;
	DrawPostScriptPolygon(OutputFile, SegX, SegY, 2,
			0, 0, 0.5,
			RGB);

	/*Draw triangle*/
	RGB[0]=1.0;	// turn to purple
	DrawPostScriptPolygon(OutputFile, PPolygonX, PPolygonY, 3,
			1, 1, 0.5,
			RGB);

	RGB[0]=0.0;	// restore to blue
	/*Draw 3 actuator lines*/
	SegX[0] = 0; SegX[1] = x1;
	SegY[0] = 0; SegY[1] = y1;
	DrawPostScriptPolygon(OutputFile, SegX, SegY, 2,
			0, 0, 0.5,
			RGB);
	SegX[0] = 1; SegX[1] = x2;
	SegY[0] = 0; SegY[1] = y2;
	DrawPostScriptPolygon(OutputFile, SegX, SegY, 2,
			0, 0, 0.5,
			RGB);
	SegX[0] = 0; SegX[1] = x3;
	SegY[0] = 1; SegY[1] = y3;
	DrawPostScriptPolygon(OutputFile, SegX, SegY, 2,
			0, 0, 0.5,
			RGB);
	
	ClosePostScriptOutput(OutputFile);
}


int main(int argc, char** argv)
{
	int i, j, ret;
	int testNum = -1; // To save several runs for automated tests

	Timer t;
	double startTime=0.0, endTime=0.0, totalTime=0.0;

	int numSols = 0;
	// It is known that the polynom to be solved is of degree 6, so there are at most 6 solutions.
	// It is better then to use a static array, and return the number of solutions to know how many to read
	double x1Arr[6], y1Arr[6], x2Arr[6], y2Arr[6], x3Arr[6], y3Arr[6]; /* Output array variables for solutions.*/
	double eps=0.0, L1=0.0, L2=0.0, L3=0.0; /* Input variables for function. */
  
	/* Get input arguments. */
	if (argc > 1)
	{
		if (argc < 5)
		{
			printf("Wrong number of parameters.\n Usage: %s [eps L1 L2 L3]\n", argv[0]);
			return -1;
		}
		eps = atof(argv[1]);
		L1 = atof(argv[2]);
		L2 = atof(argv[3]);
		L3 = atof(argv[4]);
		if (argc == 6) // To add test run number to output
			testNum = atoi(argv[5]);
	}
	else
	{
	  /*Get input from console.*/
	  printf("eps = ");
	  scanf("%lf",&eps);
	  printf("L1 = ");
	  scanf("%lf",&L1);
	  printf("L2 = ");
	  scanf("%lf",&L2);
	  printf("L3 = ");
	  scanf("%lf",&L3);
	}

	/* Check input validity.*/
	if (eps < 1.0e-5 || eps > 1.0e-2 ||
		L1 < 0.25 || L1 > 1.0 ||
		L2 < 0.25 || L2 > 1.0 ||
		L3 < 0.25 || L3 > 1.0)
	{
		printf("Argument out of range, exiting\n");
		return -2;
	}


	t = InitTimer();
	startTime = StartTime(t);
	// Runs the same problem a set number of times and finds the average run time
	for (j = 0; j < NUM_RUNS; ++j)
	{
		ret = GetForwardKinematics(eps, L1, L2, L3, x1Arr, y1Arr, x2Arr, y2Arr, x3Arr, y3Arr, &numSols);

		if (ret != 0)
		{
			printf("Something went wrong!\n");
			return -3;
		}
	}
	endTime = EndTime(t);
	totalTime = GetElapsedTimeInMilliSec(startTime, endTime);
	printf("Elapsed milli seconds for %d run(s): %lf\n", NUM_RUNS, totalTime);
	printf("Average per run milli seconds: %lf\n", totalTime/NUM_RUNS);
	
	if (testNum > 0) // For automated tests
		printf("Ran for L1=%.3lf, L2=%.3lf, L3=%.3lf. Eps=%.5lf\n", L1, L2, L3, eps);
	/* Check and print results ||B1-B2||,...,  ||B1-A1||,... */
	for (i = 0; i < numSols; ++i)
	{
		double d2;

		d2 = x1Arr[i]*x1Arr[i] + y1Arr[i]*y1Arr[i];
		printf("||B1-A1||-L1 = %lf\n", sqrt(d2)-L1);
		assert(fabs(sqrt(d2)-L1) < eps);
		if (fabs(sqrt(d2)-L1) >= eps)
			printf("Warning: result doesn't meet tolerance requirement\n");

		d2 = (x2Arr[i]-1.0)*(x2Arr[i]-1.0) + y2Arr[i]*y2Arr[i];
		printf("||B2-A2||-L2 = %lf\n", sqrt(d2)-L2);
		assert(fabs(sqrt(d2)-L2) < eps);
		if (fabs(sqrt(d2)-L2) >= eps)
			printf("Warning: result doesn't meet tolerance requirement\n");

		d2 = x3Arr[i]*x3Arr[i] + (y3Arr[i]-1.0)*(y3Arr[i]-1.0);
		printf("||B3-A3||-L3 = %lf\n", sqrt(d2)-L3);
		assert(fabs(sqrt(d2)-L3) < eps);
		if (fabs(sqrt(d2)-L3) >= eps)
			printf("Warning: result doesn't meet tolerance requirement\n");

		/* Do the same for ||B1-B2||=0.25, ||B2-B3||, ||B3-B1|| :*/
		d2 = (x1Arr[i]-x2Arr[i])*(x1Arr[i]-x2Arr[i]) + (y1Arr[i]-y2Arr[i])*(y1Arr[i]-y2Arr[i]);
		printf("||B1-B2||-0.25 = %lf\n", sqrt(d2)-0.25);
		assert(fabs(sqrt(d2)-0.25) < eps);
		if (fabs(sqrt(d2)-0.25) >= eps)
			printf("Warning: result doesn't meet tolerance requirement\n");

		d2 = (x2Arr[i]-x3Arr[i])*(x2Arr[i]-x3Arr[i]) + (y2Arr[i]-y3Arr[i])*(y2Arr[i]-y3Arr[i]);
		printf("||B2-B3||-0.25 = %lf\n", sqrt(d2)-0.25);
		assert(fabs(sqrt(d2)-0.25) < eps);
		if (fabs(sqrt(d2)-0.25) >= eps)
			printf("Warning: result doesn't meet tolerance requirement\n");

		d2 = (x3Arr[i]-x1Arr[i])*(x3Arr[i]-x1Arr[i]) + (y3Arr[i]-y1Arr[i])*(y3Arr[i]-y1Arr[i]);
		printf("||B3-B1||-0.25 = %lf\n", sqrt(d2)-0.25);
		assert(fabs(sqrt(d2)-0.25) < eps);
		if (fabs(sqrt(d2)-0.25) >= eps)
			printf("Warning: result doesn't meet tolerance requirement\n");
	}  


	/* Drawing to PostScript file.*/
	for (i=0; i<numSols; ++i)
	{
		char fname[14]; /* Filename cannot exceed 14 characters here. */
		if (testNum > 0)
			sprintf(fname, "Solution%d-%d.ps", testNum, i); // To save several runs for automated tests
		else
			sprintf(fname, "Solution%d.ps", i);
		DrawConfigurationInFile(fname,
			x1Arr[i], x2Arr[i], x3Arr[i],
			y1Arr[i], y2Arr[i], y3Arr[i]);
	}

	return 0;
}

