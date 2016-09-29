/*This file was adopted from the code of the IRIT library PS filter.   */
/*It contains several PS utilities. In order to print inside the page,*/
/*the data needs to be scaled to [-1,1].                              */

#include <stdio.h>
#include <assert.h>
#include <math.h>
#include <string.h>

#include "PostScriptOutput.h"

#pragma warning(disable : 4996) /* This line disables warning 4996, which we are aware of and do intentionally.*/

/*PS related macros.*/
#define DEFAULT_LINE_WIDTH		0.004
#define DEFAULT_SIZE			7.0
#define BBOX_BNDRY_OFFSET		4   /* Add 4 pixel boundary to bbox. */

typedef enum {
  DRAW_POINT_CROSS,
  DRAW_POINT_FULL_CIRC,
  DRAW_POINT_HOLLOW_CIRC
} DrawPointType;

static DrawPointType 
  GlblDrawPoint = DRAW_POINT_FULL_CIRC;

static int GlblBBoxClip = 0;

static double
  GlblLineWidth = DEFAULT_LINE_WIDTH,
  GlblBBoxClipDomain[4],
  GlblPointScale = 0.01,
  GlblPrintSize = DEFAULT_SIZE;

typedef struct {
  double Min[2];
  double Max[2];
} GMBBBboxStruct;

static GMBBBboxStruct GlblBbox;

static FILE *OutputFile = NULL;

static int ColorTable[][3] = {
  { /* BLACK		*/ 0,   0,   0 },
  { /* BLUE		*/ 0,   0, 255 },
  { /* GREEN		*/ 0, 255,   0 },
  { /* CYAN		*/ 0, 255, 255 },
  { /* RED		*/ 255,   0,   0 },
  { /* MAGENTA 	*/ 255,   0, 255 },
  { /* BROWN		*/ 50,   0,   0 },
  { /* LIGHTBLUE	*/ 0,   0, 255 },
  { /* LIGHTGREEN	*/ 0, 255,   0 },
  { /* LIGHTCYAN	*/ 0, 255, 255 },
  { /* LIGHTRED	*/ 255,   0,   0 },
  { /* LIGHTMAGENTA	*/ 255,   0, 255 },
};


static void RandomRGB(double* RGBColor);

static void DumpOnePoint(FILE *f,
			 double *Pt,
			 int IsVector,
			 double *RGBColor,
			 double Gray);


static void DumpOnePoly(FILE *f,
			double *PPolygonX,
      double *PPolygonY,
      int PolygonSize,
			int IsPolygon,
			int Fill,
			double Gray,
			double *RGBColor);

static void DumpDataTrailerForPostScript(FILE *OutputFile);
static void DumpDataHeaderForPostScript(FILE *OutputFile);


/*********************************************************/
static void DumpOnePoint(FILE *f,
			 double *Pt,
			 int IsVector,
			 double *RGBColor,
			 double Gray)
{
  double
  Size = GlblPointScale;

  if (RGBColor)
  	fprintf(f, "%f %f %f setrgbcolor\n",
		RGBColor[0],
		RGBColor[1],
		RGBColor[2]);
  else
	  fprintf(f, "%f setgray\n", (1.0 - Gray));

  if (IsVector) {
	  fprintf(f, "0 0 %6.4f %6.4f 0 lw\n", Pt[0], Pt[1]);
  }

  switch (GlblDrawPoint) {
	  case DRAW_POINT_FULL_CIRC:
	    /* Dump a point as a full circle. */
	    fprintf(f, "newpath %6.4f %6.4f %6.4f 0 360 arc fill\n",
		    Pt[0], Pt[1], Size); 
	    break;
	  case DRAW_POINT_HOLLOW_CIRC:
	    /* Dump a point as a hollowed circle. */
	    fprintf(f, "newpath %6.4f %6.4f %6.4f 0 360 arc fill\n",
		    Pt[0], Pt[1], Size); 
	    if (RGBColor)
		    fprintf(f, "1.0 1.0 1.0 setrgbcolor\n");
	    else
		    fprintf(f, "1.0 setgray\n");
	    fprintf(f, "newpath %6.4f %6.4f %6.4f 0 360 arc fill\n",
		      Pt[0], Pt[1], Size * 0.7); 
	    break;
	case DRAW_POINT_CROSS:
	    /* Dump a point as a cross. */
	    fprintf(f, "%6.4f %6.4f %6.4f %6.4f  %6.4f lw\n",
		    Pt[0] - Size, Pt[1],
		    Pt[0] + Size, Pt[1],
		    Size * 0.2);
	    fprintf(f, "%6.4f %6.4f %6.4f %6.4f  %6.4f lw\n",
		    Pt[0], Pt[1] - Size,
		    Pt[0], Pt[1] + Size,
		    Size * 0.2);
	    break;
    }
}

static void DumpOnePoly(FILE *f,
			double *PPolygonX,
      double *PPolygonY,
      int PolygonSize,
			int IsPolygon,
			int Fill,
			double Gray,
			double *RGBColor)
{
  int i;

  if (RGBColor) {
    fprintf(f, "%f %f %f setrgbcolor\n",
      RGBColor[0] ,
      RGBColor[1] ,
      RGBColor[2] );
  }
  else
    fprintf(f, "%f setgray\n", (1.0 - Gray));

  for (i=0; i<PolygonSize; ++i) {
    if (i == 0) {
      fprintf(f, "mark\n");
    }
    fprintf(f, "%6.4f %6.4f\n",
      PPolygonX[i], PPolygonY[i]);
  }

  if (IsPolygon && i > 0) {
    fprintf(f, "%6.4f %6.4f\n",
      PPolygonX[0], PPolygonY[0]);
  }

  fprintf(f, Fill ? "pf\n" : "pl\n");
}

void RandomRGB(double* RGBColor)
{
  static int Color = 0;
  Color = (++Color)%12;
  RGBColor[0] = ColorTable[Color][0];
  RGBColor[1] = ColorTable[Color][1];
  RGBColor[2] = ColorTable[Color][2];
}



static void DumpDataHeaderForPostScript(FILE *OutputFile)
{
  static char *PSHeader[] = {
	"%!PS-Adobe-2.0 EPSF-1.2",
	"BBOX",
	"TITLE",
	"%%Creator: irit2ps",
	"DATE",
	"%%EndComments",
	"",
	"",
	"gsave",
	"",
	"72 72 scale",	/* First scale to inches, se we can speak in inches. */
	"TRANSFORM",
	"",
	"/pl {",
	"    newpath",
	"    moveto",
        "    {counttomark 0 ne {lineto} {exit} ifelse} loop",
	"    pop % the mark",
	"    stroke",
	"} def",
	"",
	"/pf {",
	"    newpath",
	"    moveto",
        "    {counttomark 0 ne {lineto} {exit} ifelse} loop",
	"    pop % the mark",
	"    fill",
	"} def",
	"",
	"/cl {",
	"    newpath",
	"    moveto",
        "    curveto",
	"    stroke",
	"} def",
	"",
	"/cf {",
	"    newpath",
	"    moveto",
        "    curveto",
	"    fill",
	"} def",
	"",
	"/ln { % Line",
	"    newpath",
	"    moveto",
	"    lineto",
	"    stroke",
	"} def",
	"",
	"/lw { % Line with Width control",
	"    newpath",
	"    setlinewidth",
	"    moveto",
	"    lineto",
	"    stroke",
	"} def",
	"",
	"/ldg { % Line with Depth control - Gray",
	"    newpath",
	"    setgray",
	"    moveto",
	"    lineto",
	"    stroke",
	"} def",
	"",
	"/ldc { % Line with Depth control - Color",
	"    newpath",
	"    setrgbcolor",
	"    moveto",
	"    lineto",
	"    stroke",
	"} def",
	"",
	"/cw { % Curve with Width control",
	"    newpath",
	"    setlinewidth",
	"    moveto",
	"    curveto",
	"    stroke",
	"} def",
	"",
	"/cdg { % Curve with Width control - Gray",
	"    newpath",
	"    setgray",
	"    moveto",
	"    curveto",
	"    stroke",
	"} def",
	"",
	"/cdc { % Curve with Width control - Color",
	"    newpath",
	"    setrgbcolor",
	"    moveto",
	"    curveto",
	"    stroke",
	"} def",
	"",
	"1 setlinecap",
	"1 setlinejoin",
	"WIDTH",
	"",
	NULL
  };
  double
	  PrintWidth = GlblPrintSize,
	  PrintHeight = GlblPrintSize;
    int i, IPrintLeft, IPrintRight, IPrintBot, IPrintTop,
	    IPrintCntrX = (int) (72.0 * 4.0),
	    IPrintCntrY = (int) (72.0 * (1.0 + PrintHeight / 2.0)),
	    IPrintWidth2 = (int) (72.0 * PrintWidth / 2),
	    IPrintHeight2 = (int) (72.0 * PrintHeight / 2);


  if (GlblBBoxClip) {
    IPrintLeft = IPrintCntrX + ((int) (IPrintWidth2 * GlblBBoxClipDomain[0]));
	  IPrintRight = IPrintCntrX + ((int) (IPrintWidth2 * GlblBBoxClipDomain[2]));
	  IPrintBot = IPrintCntrY + ((int) (IPrintHeight2 * GlblBBoxClipDomain[1]));
	  IPrintTop = IPrintCntrY + ((int) (IPrintHeight2 * GlblBBoxClipDomain[3]));
  }
  else {
    IPrintLeft = IPrintCntrX + ((int) (IPrintWidth2 * GlblBbox.Min[0]));
    IPrintRight = IPrintCntrX + ((int) (IPrintWidth2 * GlblBbox.Max[0]));
    IPrintBot = IPrintCntrY + ((int) (IPrintHeight2 * GlblBbox.Min[1]));
    IPrintTop = IPrintCntrY + ((int) (IPrintHeight2 * GlblBbox.Max[1]));
  }

  for (i = 0; PSHeader[i] != NULL; i++) {
    if (strcmp(PSHeader[i], "BBOX") == 0) {
      fprintf(OutputFile, "%%%%BoundingBox: %d %d %d %d\n",
        IPrintLeft - BBOX_BNDRY_OFFSET,
        IPrintBot - BBOX_BNDRY_OFFSET,
        IPrintRight + BBOX_BNDRY_OFFSET,
        IPrintTop + BBOX_BNDRY_OFFSET);
    }
  	else if (strcmp(PSHeader[i], "DATE") == 0) {
      fprintf(OutputFile, "%%%%CreationDate: %s\n", "dd:mm:yy");
  	}
	  else if (strcmp(PSHeader[i], "TITLE") == 0)
	    fprintf(OutputFile, "%%%%Title: %s\n", NULL);
	  else if (strcmp(PSHeader[i], "TRANSFORM") == 0) {
	    fprintf(OutputFile, "4.0 %f translate\t\t%% center image\n",
		    1.0 + PrintHeight / 2.0);
	    fprintf(OutputFile, "%f %f scale\t\t%% maps -1..1 to this domain\n",
		    PrintWidth / 2.0, PrintHeight / 2.0);
	  }
	  else if (strcmp(PSHeader[i], "WIDTH") == 0)
	    fprintf(OutputFile, "%f setlinewidth\n",
		    GlblLineWidth > 0.0 ? GlblLineWidth
					: -GlblLineWidth * DEFAULT_LINE_WIDTH);
	  else
	    fprintf(OutputFile, "%s\n", PSHeader[i]);
  }
}

static void DumpDataTrailerForPostScript(FILE *OutputFile)
{
  int i;

  static char *PSTrailer[] = {
      "",
      "showpage",
      "grestore",
      NULL
  };

  for (i = 0; PSTrailer[i] != NULL; i++)
    fprintf(OutputFile, "%s\n", PSTrailer[i]);
}


/*   Public Functions */
static double Tx = 0.0;
static double Ty = 0.0;
static double Scale = 1.0;

FILE* OpenPostScriptOutput(char *FileName)
{
  FILE* OutputFile = NULL;
  if (FileName != NULL) {
	  if ((OutputFile = fopen(FileName, "w")) == NULL) {
      fprintf(stderr, "Failed to open \"%s\".\n", FileName);
      assert(0);
	  }
  }
  else
  	OutputFile = stdout;

  DumpDataHeaderForPostScript(OutputFile);

  return OutputFile;
}

void ClosePostScriptOutput(FILE* OutputFile)
{
  DumpDataTrailerForPostScript(OutputFile);
  fclose(OutputFile);
}



void DrawPostScriptPoint(FILE* OutputFile, double x, double y, double Gray, double RGBColor[3])
{
    double Pt[2];
    Pt[0] = x;
    Pt[1] = y;

    DumpOnePoint(OutputFile, Pt, 0, RGBColor, Gray);
}


void DrawPostScriptPolygon(FILE* OutputFile, double *PPolygonX, double *PPolygonY, int PolygonSize,
			int IsPolygon, int Fill, double Gray,
			double RGBColor[3])
{
  DumpOnePoly(OutputFile,
			PPolygonX,
      PPolygonY,
      PolygonSize,
			IsPolygon,
			Fill,
			Gray,
			RGBColor);
}
