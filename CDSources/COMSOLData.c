/*
 *  COMSOLData.h
 *
 *  Class representing the data from the COMSOL file.
 *  We store the data in the dStore arrays, one column
 *  per entry in dStore.
 *  Beyond that the data may have a hidden internal
 *  organisation as a rectangular array of points.
 *  In that case they have been organised with the
 *  most rapidly varying dimension in dStore[0] and the
 *  least rapidly in dStore[nDim-1].
 *  Information about this organization is NOT stored
 *  in the file. It can be added by calling the appropriate
 *  methods or you can try to deduce the structure by
 *  processing the data looking for patterns.
 *
 *  To make this easier for Fred to use directly I am
 *  tranlating it to pure C.
 *  To support 2D and 3D specialised sub-classes in C
 *  I have to do my own function overloading, which I do
 *  by ensuring that all methods in the base class have
 *  names that begin CD.
 *
 *  BCollett 2/14.
 *  BCollett 7/24/15 Upgrade ParseHeader to support reading files with
 *  expressions other than E.E<dir>. At least I need to support names
 *  without periods in them.
 */
#include <string.h>
#include <stdio.h>
#include <stdlib.h>
#include <math.h>

//#include "Debug.h"
#include "COMSOLData.h"

//
//  These two kludgy globals are used to pass filenames to the binary
//  writer if you wish.
//
const char* gFieldFileName = NULL;
const char* gModelFileName = NULL;


/*
 *  Class vars. deal with the issues of errors.
 *  These co-operate with the CDError enum.
 *  Converted to file statics.
 */
static const char* sgCDErrorMsg[] = { "",
                    "Unable to open input file ",
                    "Incomplete header. Error at line ",
                    "Storage allocation failed at expression ",
                    "Failed to allocate storage for expression names.",
                    "Failed to open output file "
                    };
static int sgCDErrorExtraInfo[] = { 0, 2, 1, 1, 0, 2};
static int sgCDErrorVal = 0;
static char sgCDErrorStr[256];

/****************************************************************/
//
//  ctors
//  Actually plain functions that play the role of constructors.
//  Don't need two-step in C since allocation and initialisation
//  are clearly separate.
//
/****************************************************************/
/*
 *  Init  is passed a filename. It opens the file, parses the header
 *  and then sucks the data into storage.
 */
CDError CDInit(CDData* dp, const char* fname)
{
  CDError theErr;
  int line, expr, nExpr, i;
  size_t nChar;
  //
  //  Let's try to open the file for reading.
  //
  FILE* ifp = fopen(fname, "rt");
  if (ifp == NULL) {
    fprintf(stderr, "CDInit: Failed to open file %s.", fname);
    return kCDCantOpenIn;
  }
  gFieldFileName = fname;
  //
  //  Store copy of file name.
  //
  nChar = strlen(fname);
  dp->mFileName = (char *) malloc((nChar+2) * sizeof(char));
  if (dp->mFileName == NULL) {
    fprintf(stderr, "CDInit: No space for file name %s.", fname);
    return kCDAllocFailed;
  }
  strncpy(dp->mFileName, fname, nChar);
  dp->mFileName[nChar] = 0;
  //
  //  Read in the header.
  //
  theErr = CDParseHeader(ifp, dp);
  if (theErr != kCDNoErr) {
      sgCDErrorVal = dp->mNHeadline;
      return kCDIncompleteHeader;
  }

  //
  //  Get space for data.
  //
  nExpr = dp->mNExpression + dp->mNDimension;
  dp->mDStore = (double**) malloc(nExpr * sizeof(double *));
  if (dp->mDStore == NULL) {
    fprintf(stderr, "CDInit: Failed to get space for %d expressions of data.",
            nExpr);
    return kCDAllocFailed;
  }
  for (i = 0; i < nExpr; i++) {
    dp->mDStore[i] = (double *) malloc(dp->mNLine * sizeof(double));
    if (NULL == dp->mDStore[i]) return kCDAllocFailed;
  }
  dp->mRange = (CDRange *) malloc(nExpr * sizeof(CDRange));
  //
  //  Read in the data, collecting range info as we go.
  //  Note that I split the arrays up as I pull them in.
  //
  for (expr = 0; expr < nExpr; expr++) {
    dp->mRange[expr].mMax = -DBL_MAX;
    dp->mRange[expr].mMin = DBL_MAX;
  }
  for (expr = 0; expr < nExpr; expr++) {
    printf("At %d have %d from %g to %g.\n",
           expr,
           dp->mRange[expr].mNVal,
           dp->mRange[expr].mMin,
           dp->mRange[expr].mMax);
  }
  for (line = 0; line < dp->mNLine; line++) {
    for (expr = 0; expr < nExpr; expr++) {
      double v;
      fscanf(ifp, "%lg", &v);
      dp->mDStore[expr][line] = v;
//      printf("%f,",v);
      if (expr < dp->mNDimension) {
        if (v < dp->mRange[expr].mMin) {
          dp->mRange[expr].mMin = v;
        } else if (v > dp->mRange[expr].mMax) {
          dp->mRange[expr].mMax = v;
        }
      }
    }
//  printf("\n");
  }
  for (expr = 0; expr < nExpr; expr++) {
    printf("At %d have %d from %lg to %lg by %lg\n",
           expr,
           dp->mRange[expr].mNVal,
           dp->mRange[expr].mMin,
           dp->mRange[expr].mMax,
           dp->mRange[expr].mDelta);
  }
  CDAnalyse(dp);
  return kCDNoErr;
}
//
//  CDFinish must be called after any call to CDInit once you are done
//  the data to return the storage.
//
void CDFinish(CDData* dp)
{
  int line;
  if (dp->mNExpression > 0) {
    free(dp->mExprNames);
  }
  if (dp->mNLine > 0) {
    for (line = 0; line < dp->mNDimension + dp->mNExpression; line++) {
      if (NULL != dp->mDStore[line]) {
        free(dp->mDStore[line]);
      }
    }
    free(dp->mDStore);
    free(dp->mRange);
  }
}

/****************************************************************/
//
//        Accessors
//
/****************************************************************/
/*
 *  GetValueAt takes a trio of coordinates in either index
 *  space or real space and returns the appropriate value.
 *  These do careful error checking and return a NaN if the
 *  range check fails.
 */
double CDGetValueAtIndex(CDData* dp, unsigned int dim, unsigned int index[3])
{
  int idx;
  if (dim >= dp->mNDimension + dp->mNExpression) {
    return nan("");
  }
  if ((index[0] > dp->mRange[0].mNVal) || (index[1] > dp->mRange[1].mNVal) ||
      (index[2] > dp->mRange[2].mNVal)) {
    return nan("");
  }
  idx = (index[2]*dp->mRange[1].mNVal+index[1])*dp->mRange[0].mNVal+index[0];
  return dp->mDStore[dim][idx];
}
double CDGetValueAtPoint(CDData* dp, unsigned int dim, double coord[3])
{
  int index[3], uindex[3], i, j, idxl, idxu;
  double min, max;
  if (dim >= dp->mNDimension + dp->mNExpression) {
    return nan("");
  }
  for (i = 0; i < 3; i++) {
    if ((coord[i] < dp->mRange[i].mMin) || (coord[i] > dp->mRange[i].mMax)) {
      printf("Out of range\n");
      return nan("");
    }
    if (dp->mRange[i].mNVal > 1) {
      index[i] = (coord[i] - dp->mRange[i].mMin) / dp->mRange[i].mDelta;
    } else {
      index[i] = 0;
    }
  }
  //
  //  At this point I have found the indices for the three dimensions. These
  //  are normally the indices of the coord BELOW the given coord. Use this
  //  to check that the desired point (ccord) lies inside the box formed
  //  by the indexed coords. NOTE can't do this before we have all three
  //  indices.
  //
  idxl = (index[2]*dp->mRange[1].mNVal+index[1])*dp->mRange[0].mNVal+index[0];
  for (i = 0; i < 3; i++) {
    if (dp->mRange[i].mNVal > 1) {
      for (j = 0; j < 3; j++) {
        uindex[j] = index[j];
      }
      uindex[i] = index[i]+1;
      idxu = (uindex[2]*dp->mRange[1].mNVal+uindex[1])*dp->mRange[0].mNVal+uindex[0];
      min = dp->mDStore[i][idxl];
      max = dp->mDStore[i][idxu];
      if ((coord[i] > max) || (coord[i] < min)) {
        fprintf(stderr,
                "Coord %d out of range at [%g,%g,%g], range %g->%g.\n",
                i, coord[0], coord[1], coord[2], min, max);
        fprintf(stderr, "[%d,%d,%d]=%d->[%d,%d,%d]=%d\n\n",
                index[0], index[1], index[2], idxl,
                uindex[0], uindex[1], uindex[2], idxu);
      }
    }
  }
  return dp->mDStore[dim][idxl];
}

/****************************************************************/
//
//        Write functions
//
/****************************************************************/
//
//  This writes the contents of the data store to a set of files
//  with names in the format basename_<exprname>.txt
//
CDError CDWriteBinaryTo(CDData* dp, const char* basename)
{
  char fname[256];
  FILE* ofp;
  int e;
  for (e = 0; e < dp->mNExpression + dp->mNDimension; e++) {
    sprintf(fname, "%s_%s.bin", basename, dp->mExprNames[e]);
    ofp = fopen(fname, "wb");
    if (NULL == ofp) {
      strncpy(sgCDErrorStr, fname, 255);
      return kCDCantOpenOut;
    }
    printf("Write %d doubles to %s\n", dp->mNLine, fname);
    fwrite((const char*) dp->mDStore[e], sizeof(double), dp->mNLine, ofp);
    fclose(ofp);
    printf("Closed %s\n", fname);
  }
  return kCDNoErr;
}
//
//  WriteErrorOn writes a human readable error message to
//  a FILE.
//
void CDWriteErrorOn(FILE* ofp, CDError theErr)
{
  switch (sgCDErrorExtraInfo[theErr]) {
    case 0:
      fprintf(ofp,"%s\n", sgCDErrorMsg[theErr]);
      break;

    case 1:
      fprintf(ofp,"%s%d.\n", sgCDErrorMsg[theErr], sgCDErrorVal);
      break;

    case 2:
      fprintf(ofp,"%s%s.\n", sgCDErrorMsg[theErr], sgCDErrorStr);
      break;

    default:
      break;
  }
}
/****************************************************************/
//
//        Internal Helpers
//
/****************************************************************/
/*
 *  Get info from the file header.
 *  Start by pulling out the info about the numbers of
 *  dimensions, lines, and expressions in the file.
 *  Then add the ability to extract the expression
 *  names for use in building filenames.
 */
CDError CDParseHeader(FILE* ifp, CDData* dp)
{
  char headline[255];
  char option[32];
  static char modelName[256];
  char ch;
  int expr;
  //
  //  File is open. Read the header. Should be a 9 line
  //  header of comments beginning %.
  //  Among them should be lines
  //    Dimension: <number of dimensions, 2 or 3>
  //    Nodes: <number of lines of data in file>
  //    Expressions: <number of expressions. Have 3+nE # per line>
  //  All of those need to be parsed.
  //  Every line but the last has the format
  //  % <option>: <argument>
  //  where, depending on the option, <argument> may be a number of a
  //  string.
  //
  dp->mNHeadline = 0;
  while ((ch = fgetc(ifp)) == '%') {
    dp->mNHeadline++;
    fscanf(ifp, "%s", option);
    printf("Found header option %s on line %d.\n", option, dp->mNHeadline);
    if (strcmp(option, "Dimension:") == 0) {
      fscanf(ifp, "%d", &(dp->mNDimension));
    } else if (strcmp(option, "Nodes:") == 0) {
      fscanf(ifp, "%d", &(dp->mNLine));
    } else if (strcmp(option, "Expressions:") == 0) {
      fscanf(ifp, "%d", &(dp->mNExpression));
    } else if (strcmp(option, "Model:") == 0) {
      fscanf(ifp, "%s", modelName);
      gModelFileName = modelName;
    } else if (dp->mNHeadline == 9) {
      //
      //  Parse the line for variable names.
      //  Expect nDim 1 char names followed by
      //  nExpr names of form <var>.<comp> <units>
      //
      int nName = dp->mNExpression + dp->mNDimension;
      char* nameBuff = (char *) malloc(nName * 16 * sizeof(char));
      if (NULL == nameBuff) return kCDNameAllocFailed;
      dp->mExprNames = (char **) malloc(nName * sizeof(char *));
      if (NULL == dp->mExprNames) return kCDNameAllocFailed;
      for (expr = 0; expr < nName; expr++) {
        dp->mExprNames[expr] = &(nameBuff[16 * expr]);
      }
      strcpy(dp->mExprNames[0], option); // option already holds first name
      for (expr = 1; expr < dp->mNDimension; expr++) {
        fscanf(ifp, "%s", dp->mExprNames[expr]);
      }
      for (; expr < nName; expr++) {
        unsigned int nfirst;
        //
        //  These are more difficult because there are two strings
        //  for each (expression and unit) and the name is only
        //  the second part of the first string.
        //  Altered BC 7/24/15 to support a wider range of expression so
        //  we can read files containing the potential.
        //
//        fscanf(ifp, "%s", option);
//        printf("%s\n", option);
//        nfirst = (unsigned int) strcspn(option, ".");
//        printf("%s\n", &option[nfirst+1]);
//        strcpy(dp->mExprNames[expr], &option[nfirst+1]);
        fscanf(ifp, "%s", dp->mExprNames[expr]);
        fscanf(ifp, "%s", option);  // Skip over units
      }
//      if (gDebug) {
        printf("Found expressions\n");
        for (expr = 0; expr < nName; expr++) {
          printf("%s\n", dp->mExprNames[expr]);
        }
      }
//    }
    fgets(headline, 250, ifp);
    printf("Skipped over %s\n", headline);
  }
  ungetc(ch, ifp);
  return kCDNoErr;
}
//
//  Analyse figures out the structure of the data and saves
//  the info in the mRange structures.
//
void CDAnalyse(CDData* dp)
{
  int nRep[3], d;
  double v;
  int nPoint;
  /*
   *  Let's see if we can figure out the grid structure of the file.
   *  Any dimension can either be fixed or can vary.
   *  For dimensions that do vary, the lower the dimension number
   *  the faster the variation should be.
   *  We found the ranges of the dimensions as we read them in,
   *  so figure out which ones are active.
   */
  for (d = 0; d < dp->mNDimension; d++) {
    if (dp->mRange[d].mMax - dp->mRange[d].mMin > 0) {
      dp->mRange[d].mActive = true;
      v = dp->mDStore[d][0];
      nRep[d] = 0;
      while (dp->mDStore[d][nRep[d]] == v) {
        nRep[d]++;
      }
    } else {
      dp->mRange[d].mActive = false;
      nRep[d] = dp->mNLine;
    }
  }
  //
  //  Now the slow actives are set with the fastest changing one set to 1.
  //  Work backwards over the dimensions figuring our how many different
  //  values each takes.
  //
  nPoint = dp->mNLine;
  for (d = dp->mNDimension-1; d >= 0; d--) {
    if (dp->mRange[d].mActive) {
      dp->mRange[d].mNVal = nPoint / nRep[d];
      nPoint = nRep[d];
    } else {
      dp->mRange[d].mNVal = 1;
    }
    if (dp->mRange[d].mNVal > 1) {
      dp->mRange[d].mDelta = (dp->mRange[d].mMax - dp->mRange[d].mMin)/
                                  (dp->mRange[d].mNVal - 1);
    } else {
      dp->mRange[d].mDelta = 0.0;
    }
  }

}
