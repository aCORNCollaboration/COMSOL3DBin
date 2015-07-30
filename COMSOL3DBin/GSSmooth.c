//
//  GSSmooth.c
//  COMSOL3DBin
//
//  Routine to Gauss-Seidel smooth a 3D electric field array.
//  Builds a type array from a geometry file.
//
//  Created by Brian Collett on 7/29/15.
//  Copyright (c) 2015 Brian Collett. All rights reserved.
//

#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#include <stdint.h>
#include <inttypes.h>
#include <math.h>
#include "assert.h"
#include "GSSmooth.h"
#include "CD3List.h"
#include "Geometries.h"

uint8_t* NewTypeArray(uint32_t nVal[3]);
void SmoothPrintOn(CD3Data* d, uint8_t* type, FILE* ofp);

int GSSmooth(const char* fname, CD3Data* dp, int nPass)
{
  int errCode = 0, pass, comp;
  uint8_t* pointType = NULL;
  CD3List gList;
  double err = 0.0;
  uint64_t idx, idy, idz, index, rIndex;
  double newVal;
  double terr, wa, wx, wy, wz;
  double* a = dp->mField;
  /*
   *  Start by computing offsets for the neighbours. The x neigbours are
   *  always 1 higher and 1 lower in index but the y neighbours are a
   *  distance mNVal[0] (the number of possible x values) away while the
   *  z neighbours are even further.
   *  These are all multiplied by three because there are 3 components
   *  per vector.
   */
  const uint64_t dx = 3;
  const uint64_t dy = 3 * dp->mNVal[0];
  const uint64_t dz = 3 * dp->mNVal[0] * dp->mNVal[1];
  assert(NULL != dp);
  assert(NULL != a);
  //
  //  First make sure that we have a 3D leaf array.
  //
  if  (dp->mNSubField > 0) {
    fprintf(stderr, "Field is not a leaf.\n");
    return kCDNotLeaf;
  }
  //
  // Have to make sure that the data are suitable for 4-way
  // averaging.
  // First, the data must be 3D!
  //
  if (dp->mStride != 0) {
    fprintf(stderr, "Attempt to average 2D field %s.\n", dp->mFieldName);
    return kCDNot4Fold;
  }
  //
  //  Now build the point type array.
  //
  pointType = NewTypeArray(dp->mNVal);
  if (NULL == pointType) {
    fprintf(stderr, "Attempt to get storage for type array failed.\n");
    return kCDAllocFailed;
  }
//  SmoothPrintOn(dp, pointType, stdout);
  //
  //  Read in the geometry
  //
  CD3ListInit(&gList);
  if (CD3ListReadGeom(&gList, fname)) {
    CD3ListAddGeomTo(&gList, pointType, dp);
  } else {
    fprintf(stderr, "Cannot read geometry from file %s.\n", fname);
    return kCDBadGeom;
  }
  //
  //  And run the smooth.
  //
  /*
   *  Weightings to compensate for non-isotropic grid.
   */
  wa = 1.0 / (1.0/(dp->mDelta[0]*dp->mDelta[0]) +
                    1.0/(dp->mDelta[1]*dp->mDelta[1]) +
                    1.0/(dp->mDelta[2]*dp->mDelta[2]));
  wx = wa / (2.0 *dp->mDelta[0]*dp->mDelta[0]);
  wy = wa / (2.0 *dp->mDelta[1]*dp->mDelta[1]);
  wz = wa / (2.0 *dp->mDelta[2]*dp->mDelta[2]);
  /*
   *  Now we do the fancy red-black scanning.
   *  This is made a little more complex because we have to
   *  do all three components at each point.
   */
  for (pass = 0; pass < nPass; pass++) {
    err = 0.0;
    for (idz = 0; idz < dp->mNVal[2]; idz ++) {  // Red
      for (idy = 0; idy < dp->mNVal[1]; idy ++) {  // Red
        for (idx = ((idy + idz) & 1); idx < dp->mNVal[0]; idx += 2) {  // Red
          index = idz * dz + idy * dy + idx*dx;
          rIndex = index / 3;
          switch (pointType[rIndex]) {
            case 0:
              break;
              
            case 1:
              for (comp = 0; comp < 3; comp++) {
                newVal = wx*(a[index+dx+comp]+a[index-dx+comp])+
                wy*(a[index+dy+comp]+a[index-dy+comp])+
                wz*(a[index+dz+comp]+a[index-dz+comp]);
                terr = newVal - a[index+comp];
                err += terr * terr;
                //              fprintf(gDebugFile, "V[%lld]=%6.3lf -> %6.3lf\n", index, mData[index], newVal);
                a[index+comp] = newVal;
              }
              break;
              
            default:
              fprintf(stderr, "Unknown element type at index %" PRIu64 ".\n", idx);
          }
        }
      }
    }
//    SmoothPrintOn(dp, pointType, stdout);
    for (idz = 0; idz < dp->mNVal[2]; idz ++) {  // Black
      for (idy = 0; idy < dp->mNVal[1]; idy ++) {  // Black
        for (idx = ((idy + idz + 1) & 1); idx < dp->mNVal[0]; idx += 2) {  // Black
          index = idz * dz + idy * dy + idx*dx;
          rIndex = index / 3;
          switch (pointType[rIndex]) {
            case 0:
              break;
              
            case 1:
              for (comp = 0; comp < 3; comp++) {
                newVal = wx*(a[index+dx+comp]+a[index-dx+comp])+
                wy*(a[index+dy+comp]+a[index-dy+comp])+
                wz*(a[index+dz+comp]+a[index-dz+comp]);
                terr = newVal - a[index+comp];
                err += terr * terr;
                //              fprintf(gDebugFile, "V[%lld]=%6.3lf -> %6.3lf\n", index, mData[index], newVal);
                a[index+comp] = newVal;
              }
              break;
              
            default:
              fprintf(stderr, "Unknown element type at index %" PRIu64 ".\n", idx);
          }
        }
      }
    }
    fprintf(stderr, "Pass %d error = %lf.\n", pass, err);
  }
//  SmoothPrintOn(dp, pointType, stdout);
  return errCode;
}
/*
 *  Create a new type array with its guts set to active and its border
 *  to inactive.
 */
uint8_t* NewTypeArray(uint32_t nVal[3]) {
  uint8_t* pointType = NULL;
  uint64_t arraySize;
  assert(NULL != nVal);
  uint32_t ix, iy, iz;
  uint32_t ixMax, iyMax, izMax;
  //
  //  Now get space for the point type array.
  //  Build size a bit at a time to avoid overflow.
  //
  arraySize = nVal[0];
  arraySize *= nVal[1];
  arraySize *= nVal[2];
  pointType = (uint8_t*) malloc(arraySize * sizeof(uint8_t));
  if (NULL == pointType) {
    return NULL;
    
  }
  ixMax = nVal[0]-1;
  iyMax = nVal[1]-1;
  izMax = nVal[2]-1;
  /*
   *  Set the guts to active and the border to inactive.
   */
  memset(pointType, 1, arraySize);
  for (iy = 0; iy <= iyMax; iy++) {     // Mark z ends inactive
    for (ix = 0; ix <= ixMax; ix++) {
      pointType[(izMax * nVal[1] + iy)*nVal[0]+ix] = 0;
      pointType[(0 * nVal[1] + iy)*nVal[0]+ix] = 0;
    }
  }
  for (iz = 0; iz <= izMax; iz++) {     // Mark y ends inactive
    for (ix = 0; ix <= ixMax; ix++) {
      pointType[(iz * nVal[1] + 0)*nVal[0]+ix] = 0;
      pointType[(iz * nVal[1] + iyMax)*nVal[0]+ix] = 0;
    }
  }
  for (iz = 0; iz <= izMax; iz++) {     // Mark x ends inactive
    for (iy = 0; iy <= iyMax; iy++) {
      pointType[(iz * nVal[1] + iy)*nVal[0]+0] = 0;
      pointType[(iz * nVal[1] + iy)*nVal[0]+ixMax] = 0;
    }
  }
  return pointType;
}

void SmoothPrintOn(CD3Data* d, uint8_t* type, FILE* ofp)
{
  int i,j,k, c;
  double* a = d->mField;
  for (k = 0; k < d->mNVal[2]; k++) {
    fprintf(ofp, "k=%d\n", k);
    for (j = 0; j < d->mNVal[1]; j++) {
      for (i = 0; i < d->mNVal[0]; i++) {
        fprintf(ofp, "%d", type[((k*d->mNVal[1] + j)*d->mNVal[0] +
                                  i)]);
      }
      fprintf(ofp, "   ");
      for (c = 0; c < 3; c++) {
        for (i = 0; i < d->mNVal[0]; i++) {
          fprintf(ofp, "%6.3f ", a[((k*d->mNVal[1] + j)*d->mNVal[0] +
                                    i)*3 + c]);
        }
        fprintf(ofp, "   ");
      }
      fprintf(ofp, "\n");
    }
  }
}
