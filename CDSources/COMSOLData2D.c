/*
 *  COMSOLData2D.cpp
 *  COMSOLReader2
 *
 *  This is a specialisation of the COMSOLReader class that only
 *  supports 2D data files containing 2D fields. It checks that
 *  the file has two active spatial dimensions though it expects
 *  three real spatial dimensions.
 *  It expects 2 field dimensions and expects them to be the same
 *  ones as are active in the spatial field.
 *
 *  Created by Brian Collett on 3/11/14.
 *  NOTE that this is NOT as general as it could be. It is pretty much
 *  intended to read 2D slices of 3D axi-symmetric models.
 *  Copyright (c) 2014 Brian Collett. All rights reserved.
 *  BCollett 3/14/14 OBSOLETE; Merged with ComsolData3D.
 */

#include <string.h>
#include <stdlib.h>
#include <math.h>
#include "COMSOLData2D.h"

//
//  Define this if you want to bounds check every value.
//
#define CD2BoundsCheck 1

//
//  Init fills in the data structure using the information in the file.
//  BEWARE: CDInit allocates a lot of storage. We MUST ensure that that
//  storage gets disposed of before we leave. This means that once CDInit
//  has run ALL error must go through an exit code path, hence the rare
//  use of a goto and label.
//
CDError CD2Init(CD2Data* dp, const char* fname)
{
  CDError theErr;
  int dim, nDim;
  int inactiveDim = -1;
  int nInactive = 0;
  //
  //  Start by constructing a CDData from the file.
  //
  CDData cData;
  theErr = CDInit(&cData, fname);
  if (theErr != kCDNoErr) {
    goto ErrorExit;
  }
  //
  //  Have a real file nicely parsed out. See if it is a 2D file.
  //  We need three dimensions and two expressions.
  //  I did them on separate lines to make debugging a little easier.
  //
  if (cData.mNDimension != 3) {
    fprintf(stderr, "Expected three dimensions, found %d.\n",
            cData.mNDimension);
    theErr = kCDBadStructure;
    goto ErrorExit;
  }
  if (cData.mNExpression != 2) {
    fprintf(stderr, "Expected two expressions, found %d.\n",
            cData.mNExpression);
    theErr = kCDBadStructure;
    goto ErrorExit;
  }
  //
  //  Count how many dimensions are active and, if only one, figure
  //  out which is not. Meanwhile fill in info for active dims and
  //  mark cData so not thrown away.
  //
  for (nDim = 0, dim = 0; dim < 3; dim++) {
    if (!cData.mRange[dim].mActive) {
      nInactive++;
      if (nInactive > 1) {
        fprintf(stderr, "Should be only one inactive dimension, found %d.\n",
                nInactive);
        theErr = kCDBadStructure;
        goto ErrorExit;
      }
      inactiveDim = dim;
    } else {
      dp->mNVal[nDim] = cData.mRange[dim].mNVal;
      dp->mMin[nDim] = cData.mRange[dim].mMin;
      dp->mMax[nDim] = cData.mRange[dim].mMax;
      dp->mDelta[nDim] = cData.mRange[dim].mDelta;
      nDim++;
    }
  }
  dp->mFieldVals[0] = cData.mDStore[3];
  cData.mDStore[3] = NULL;
  dp->mFieldVals[1] = cData.mDStore[4];
  cData.mDStore[4] = NULL;
  if (nInactive != 1) {
    fprintf(stderr, "Should be exactly one inactive dimension, found %d.\n",
            nInactive);
    theErr = kCDBadStructure;
    goto ErrorExit;
  }
  //
  //  Now can check the field components.
  //
  switch (inactiveDim) {
    case 0:
      //
      //  x inActive. Check for Ey and Ez.
      //
      if (strcmp(cData.mExprNames[3], "es.Ey") != 0) {
        fprintf(stderr, "First expression '%s' should be 'es.Ey'.\n",
                cData.mExprNames[3]);
        theErr = kCDBadStructure;
        goto ErrorExit;
      }
      if (strcmp(cData.mExprNames[4], "es.Ez") != 0) {
        fprintf(stderr, "Second expression '%s' should be 'es.Ey'.\n",
                cData.mExprNames[4]);
        theErr = kCDBadStructure;
        goto ErrorExit;
      }
      break;
      
    case 1:
      //
      //  y inActive. Check for Ex and Ez.
      //
      if (strcmp(cData.mExprNames[3], "es.Ex") != 0) {
        fprintf(stderr, "First expression '%s' should be 'es.Ex'.\n",
                cData.mExprNames[3]);
        theErr = kCDBadStructure;
        goto ErrorExit;
      }
      if (strcmp(cData.mExprNames[4], "es.Ez") != 0) {
        fprintf(stderr, "Second expression '%s' should be 'es.Ez'.\n",
                cData.mExprNames[4]);
        theErr = kCDBadStructure;
        goto ErrorExit;
      }
      break;
      
    case 2:
      //
      //  z inActive. Check for Ex and Ey.
      //
      if (strcmp(cData.mExprNames[3], "es.Ex") != 0) {
        fprintf(stderr, "First expression '%s' should be 'es.Ex'.\n",
                cData.mExprNames[3]);
        theErr = kCDBadStructure;
        goto ErrorExit;
      }
      if (strcmp(cData.mExprNames[4], "es.Ey") != 0) {
        fprintf(stderr, "Second expression '%s' should be 'es.Ey'.\n",
                cData.mExprNames[4]);
        theErr = kCDBadStructure;
        goto ErrorExit;
      }
      break;
      
    default:
      break;
  }
ErrorExit:
  //
  //  This makes sure that we dispose of the remaining storage in CData.
  //  ALL code paths exit through here (hence the label).
  //
  CDFinish(&cData);
  return theErr;
}
//
//  Finish tidies up after us, releasing our storage. Every call of
//  COMSOLDataInit should be balanced by a call to COMSOLDataFinish.
//
void CD2Finish(CD2Data* dp)
{
  int dim;
  for (dim = 0; dim < 3; dim++) {
    if (NULL != dp->mFieldVals[dim]) {
      free(dp->mFieldVals[dim]);
    }
  }
}
//
//  Accessors.
//
double CD2GetValueAtIndex(CD2Data* dp, unsigned int dim, unsigned int index[2])
{
  if (dim > 1) {
    return nan("");
  }
  if ((index[0] > dp->mNVal[0]) || (index[1] > dp->mNVal[1])) {
    return nan("");
  }
  int idx = ((index[1])*dp->mNVal[1])+index[0];
  return dp->mFieldVals[dim][idx];
  
}
//
//  Only have to do 2D interpolation with 2D data.
//
double CD2GetValueAtPoint(CD2Data* dp, unsigned int dim, double coord[2])
{
  int index[2], i;
  double rc[2], irc[2];                 // Reduced coords and inverses
  double minc[2];                       // Minima of surrounding box
  int idx00, idx01, idx10, idx11;
  double c0, c1, c; // Interpolation steps
  if (dim > 1) {
    return nan("");
  }
  for (i = 0; i < 2; i++) {
    index[i] = (int) ((coord[i] - dp->mMin[i]) / dp->mDelta[i]);
#ifdef CD2BoundsCheck
    if ((coord[i] < dp->mMin[i]) || (coord[i] > dp->mMax[i])) {
      fprintf(stderr, "CD2GetValueAtPoint:Coordinate %d out of range\n", i);
      return nan("");
    }
    if (index[i] >= dp->mNVal[i]) {
      fprintf(stderr, "CD2GetValueAtPoint:Index %d out of range\n", i);
      return nan("");
    }
#endif
  }
  //
  //  At this point I have found the indices for the three dimensions. These
  //  are normallly the indices of the coord BELOW the given coord. Use this
  //  to compute the array indices (idx's) of the eight points that surround
  //  the cell containing the coordinate.
  //  They have form idx<z><y><x>.
  //  NOTE can't do this before we have all three
  //  indices.
  //
  idx00 = (index[1])*dp->mNVal[0] + index[0];
  idx01 = (index[1])*dp->mNVal[0] + index[0] + 1;
  idx10 = (index[1]+1)*dp->mNVal[0] + index[0];
  idx11 = (index[1]+1)*dp->mNVal[0] + index[0] + 1;
  //
  //  Next have to find where in each dimension of the box the coordinate is.
  //  This produces a set of reduced coords expressed as a fraction of the
  //  way from the lower to upper point.
  //
  for (i = 0; i < 2; i++) {
    minc[i] = dp->mMin[i] + index[i] * dp->mDelta[i];
    rc[i] = (coord[i] - minc[i])/dp->mDelta[i];
    irc[i] = 1.0 - rc[i];
#ifdef CDBoundsCheck
    if ((rc < 0.0) || (rc > 1.0)) {
      fprintf(stderr, "CD2GetValueAtPoint:Reduced coord %d out of range.\n", i);
      return nan("");
    }
#endif
  }
  //
  //  Now we can do the interpolation.
  //
  c0 = irc[0]*dp->mFieldVals[dim][idx00] + rc[0]*dp->mFieldVals[dim][idx01];
  c1 = irc[0]*dp->mFieldVals[dim][idx10] + rc[0]*dp->mFieldVals[dim][idx11];
  c = irc[1] * c0 + rc[1] * c1;
  return c;
}

double CD2GetErAtPoint(CD2Data* dp, double coord[2])
{
  return CD2GetValueAtPoint(dp, 0, coord);
}

double CD2GetEzAtPoint(CD2Data* dp, double coord[2])
{
  return CD2GetValueAtPoint(dp, 1, coord);
}

bool CD2GetEAtPoint(CD2Data* dp, double coord[3], double* EField)
{
  int index[2], i;
  double rc[2], irc[2];                 // Reduced coords and inverses
  double minc[2];                       // Minima of surrounding box
  int idx00, idx01, idx10, idx11;
  double c0, c1; // Interpolation steps
  for (i = 0; i < 2; i++) {
    index[i] = (int) ((coord[i] - dp->mMin[i]) / dp->mDelta[i]);
#ifdef CD2BoundsCheck
    if ((coord[i] < dp->mMin[i]) || (coord[i] > dp->mMax[i])) {
      fprintf(stderr, "CD2GetValueAtPoint:Coordinate %d out of range\n", i);
      return false;
    }
    if (index[i] >= dp->mNVal[i]) {
      fprintf(stderr, "CD2GetValueAtPoint:Index %d out of range\n", i);
      return false;
    }
#endif
  }
  //
  //  At this point I have found the indices for the three dimensions. These
  //  are normallly the indices of the coord BELOW the given coord. Use this
  //  to compute the array indices (idx's) of the eight points that surround
  //  the cell containing the coordinate.
  //  They have form idx<z><y><x>.
  //  NOTE can't do this before we have all three
  //  indices.
  //
  idx00 = (index[1])*dp->mNVal[0] + index[0];
  idx01 = (index[1])*dp->mNVal[0] + index[0] + 1;
  idx10 = (index[1]+1)*dp->mNVal[0] + index[0];
  idx11 = (index[1]+1)*dp->mNVal[0] + index[0] + 1;
  //
  //  Next have to find where in each dimension of the box the coordinate is.
  //  This produces a set of reduced coords expressed as a fraction of the
  //  way from the lower to upper point.
  //
  for (i = 0; i < 2; i++) {
    minc[i] = dp->mMin[i] + index[i] * dp->mDelta[i];
    rc[i] = (coord[i] - minc[i])/dp->mDelta[i];
    irc[i] = 1.0 - rc[i];
#ifdef CDBoundsCheck
    if ((rc < 0.0) || (rc > 1.0)) {
      fprintf(stderr, "CD2GetValueAtPoint:Reduced coord %d out of range.\n", i);
      return nan("");
    }
#endif
  }
  //
  //  Now we can do the interpolations.
  //
  c0 = irc[0]*dp->mFieldVals[0][idx00] + rc[0]*dp->mFieldVals[0][idx01];
  c1 = irc[0]*dp->mFieldVals[0][idx10] + rc[0]*dp->mFieldVals[0][idx11];
  EField[0] = irc[1] * c0 + rc[1] * c1;
  c0 = irc[0]*dp->mFieldVals[1][idx00] + rc[0]*dp->mFieldVals[1][idx01];
  c1 = irc[0]*dp->mFieldVals[1][idx10] + rc[0]*dp->mFieldVals[1][idx11];
  EField[1] = irc[1] * c0 + rc[1] * c1;
  return true;
  
}
//
//  Last one treats the field as a defining slice for an axi-symmetric
//  field and returns fully 3D values from 3D points.
//
bool CD2AxGetEAtPoint(CD2Data* dp, double coord[3], double* EField)
{
  double coord2D[2];
  double field2D[2];
  double sinval = 0.0, cosval = 0.0;
  //
  //  Start by mapping from #d point to 2D point.
  //
  double r = sqrt(coord[0]*coord[0] + coord[1]*coord[1]);
  if (r > 0.0) {
    sinval = coord[1]/r;
    cosval = coord[0]/r;
  }
  coord2D[0] = r;
  coord2D[1] = coord[2];
  //
  //  Get the 2D field.
  //
  if (!CD2GetEAtPoint(dp, coord2D, field2D)) {
    return false;
  }
  //
  //  Map back to 3D.
  //
  EField[0] = field2D[0] * cosval;
  EField[1] = field2D[0] * sinval;
  EField[2] = field2D[1];
  return true;
}
