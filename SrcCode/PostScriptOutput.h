#ifndef POSTSCRIPT_OUTPUT_H
#define POSTSCRIPT_OUTPUT_H

#include <stdio.h>

FILE* OpenPostScriptOutput(char *FileName);
void ClosePostScriptOutput(FILE* OutputFile);


/* 
Point should be scaled to be between [-1,1]. 
Gray - value between 0..1
RGBColor - array of three values between 0..1 (can be NULL in which case draws in gray).
*/
void DrawPostScriptPoint(FILE* OutputFile, double x, double y, double Gray, double RGBColor[3]);

/*
Polygon should all be scaled to be between [-1,1].
RGBColor - array of three values between 0..1 (can be NULL).
*/
void DrawPostScriptPolygon(FILE* OutputFile, double *PPolygonX, double *PPolygonY, int PolygonSize,
			int IsPolygon, int Fill, double Gray,
			double RGBColor[3]);


#endif

