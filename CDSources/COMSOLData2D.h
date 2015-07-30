/*
 *  COMSOLData2D.h
 *  COMSOLReader2
 *
 *  This is a specialisation of the COMSOLReader class that only
 *  supports 2D data files containing 2D fields. It checks that
 *  the file has two active spatial dimensions though it may well
 *  have three real spatial dimensions.
 *  It expects 2 field dimensions and expects them to be the same
 *  ones as are active in the spatial field.
 *
 *  Created by Brian Collett on 3/11/14.
 *  Copyright (c) 2014 Brian Collett. All rights reserved.
 */

#ifndef __COMSOLReader2__COMSOLData2D__
#define __COMSOLReader2__COMSOLData2D__

#include "COMSOLData.h"

//
//  Because we understand the structure of this kind of file much
//  better we can have a streamlined class. For example, we no longer
//  need to store the coordinate data.
//
typedef struct CD2DataTag {
  //
  //  Three sets of data arrays, one for each dimension of either the
  //  field or the coordinate space.
  //
  double* mFieldVals[2];    // Field data
  unsigned int mNVal[3];    // Number of different indices in this dim
  double mMax[3];           // Max values of each coord.
  double mMin[3];           // and corresponding minima
  double mDelta[3];         // Deltas needed for coord conversion
} CD2Data;
//
//  Init fills in the data structure using the informaiton in the file.
//
CDError CD2Init(CD2Data* dp, const char* fname);
//
//  Finish tidies up after us, releasing our storage. Every call of
//  CD3Init should be balanced by a call to CD3Finish.
//
void CD2Finish(CD2Data* dp);
//
//  Accessors.
//  First ones treat the field as a simple 2D slice.
//
double CD2GetValueAtIndex(CD2Data* dp, unsigned int dim, unsigned int index[2]);
double CD2GetValueAtPoint(CD2Data* dp, unsigned int dim, double coord[2]);
double CD2GetErAtPoint(CD2Data* dp, double coord[2]);
double CD2GetEzAtPoint(CD2Data* dp, double coord[2]);
bool CD2GetEAtPoint(CD2Data* dp, double coord[2], double* EField);
//
//  Last one treats the field as a defining slice for an axi-symmetric
//  field and returns fully 3D values from 3D points.
//
double CD2AxGetValueAtPoint(CD2Data* dp, unsigned int dim, double coord[3]);
bool CD2AxGetEAtPoint(CD2Data* dp, double coord[3], double* EField);
//
#endif /* defined(__COMSOLReader2__COMSOLData2D__) */
