//
//  CD3Smoothable.h
//  COMSOLSmooth
//
//  This is a class that represents a COMSOL field dataset with an additional
//  control array (the pointType array) and that supports Gauss-Seidel
//  smoothing of the field data.
//  By default the pointType array allows all interior points to be smoothed
//  but keeps the boundaries fixed. Optionally you can add geometry information
//  to freeze interior points as well.
//  In part this acts as a wrapper for the C only CD3Data pseudo-class.
//
//  Created by Brian Collett on 8/5/15.
//  Copyright (c) 2015 Brian Collett. All rights reserved.
//
#include <stdio.h>
#include <stdint.h>
#include <string.h>
#include <math.h>
#include <inttypes.h>
#include "assert.h"
#include "COMSOLData3D.h"
#include "CD3Smoothable.h"
#include "CD3List.h"

//
//  Helper class.
//
PointArray::PointArray() : mNPoint(0), mArray(nullptr)
{
  
}
PointArray::PointArray(uint64_t nx, uint64_t ny, uint64_t nz) :
mNPoint(0), mArray(nullptr)
{
  SetSize(nx, ny, nz);
}

PointArray::~PointArray()
{
  if (nullptr != mArray) {
    delete mArray;
  }
}
//
//  SetSize either from explicit size or from a CD3Data.
//
CDError PointArray::SetSize(const CD3Data& d)
{
  return SetSize(d.mNVal[0], d.mNVal[1], d.mNVal[2]);
}

CDError PointArray::SetSize(uint64_t nx, uint64_t ny, uint64_t nz)
{
  uint64_t ix, iy, iz;
  mNVal[0] = nx;
  mNVal[1] = ny;
  mNVal[2] = nz;
  mNPoint = nx * ny * nz;
  mArray = new uint8_t[mNPoint];
  if (nullptr == mArray) {
    fprintf(stderr, "PointArray::SetSize failed to allocate space.\n");
    return kCDAllocFailed;
  }
  uint64_t ixMax = nx-1;
  uint64_t iyMax = ny-1;
  uint64_t izMax = nz-1;
  /*
   *  Set the guts to active and the border to inactive.
   */
  memset(mArray, 1, mNPoint);
  for (iy = 0; iy <= iyMax; iy++) {     // Mark z ends inactive
    for (ix = 0; ix <= ixMax; ix++) {
      mArray[(izMax * mNVal[1] + iy)*mNVal[0]+ix] = 0;
      mArray[(0 * mNVal[1] + iy)*mNVal[0]+ix] = 0;
    }
  }
  for (iz = 0; iz <= izMax; iz++) {     // Mark y ends inactive
    for (ix = 0; ix <= ixMax; ix++) {
      mArray[(iz * mNVal[1] + 0)*mNVal[0]+ix] = 0;
      mArray[(iz * mNVal[1] + iyMax)*mNVal[0]+ix] = 0;
    }
  }
  for (iz = 0; iz <= izMax; iz++) {     // Mark x ends inactive
    for (iy = 0; iy <= iyMax; iy++) {
      mArray[(iz * mNVal[1] + iy)*mNVal[0]+0] = 0;
      mArray[(iz * mNVal[1] + iy)*mNVal[0]+ixMax] = 0;
    }
  }
  return kCDNoErr;
}

//
//  Accessors.
//
uint8_t& PointArray::operator()(uint64_t idx0)
{
  assert(idx0 < mNPoint);
  return mArray[idx0];
}

uint8_t& PointArray::operator()(uint64_t idx1, uint64_t idx2, uint64_t idx3)
{
  assert(idx1 < mNVal[0]);
  assert(idx2 < mNVal[1]);
  assert(idx3 < mNVal[2]);
  return mArray[(idx3 * mNVal[1] + idx2) * mNVal[0] + idx1];
}

const char* CD3Smoothable::gDefName ="";

CD3Smoothable::CD3Smoothable() : mPointType()
{
  int i;
  //
  //  Filling in default elems in the CD3Data.
  //
  mData.mType = kCD3Error;
  for (i = 0; i < 3; i++) {
    mData.mNVal[i] = 0;
    mData.mMin[i] = DBL_MAX;
    mData.mMax[i] = -DBL_MAX;
    mData.mDelta[i] = nan("");
  }
  mData.mStride = 0;
  mData.mNSubField = 0;
  for (i = 0; i < kNSub; i++) {
    mData.mSubField[0] = nullptr;
  }
  mData.mField = NULL;
  mData.mFieldName = gDefName;
}
//
//  2-step construction for error handling.
//
CDError CD3Smoothable::ReadBinary(const char* filename) {
  //
  //  Start by reading in the data if possible.
  //
  FILE* ifp = fopen(filename, "rb");
  if (nullptr == ifp) {
    fprintf(stderr,
            "CD3Smoothable::ReadBinary failed to open file %s for reading.\n",
            filename);
    return kCDCantOpenIn;
  }
  if (!CD3ReadBinary(&mData, ifp)) {
    fprintf(stderr,
            "CD3Smoothable::ReadBinary failed to read file %s.\n",
            filename);
    fclose(ifp);
    return kCDCantOpenIn;
  }
  fclose(ifp);
  mData.mFieldName = filename;
  //
  //  If we got here then we have a valid mData and can use it to build
  //  ourselves a point array.
  //
  return mPointType.SetSize(mData);
}

CD3Smoothable::~CD3Smoothable()
{
  CD3Finish(&mData);
}
//
//  Can add geometry info, possibly several.
//
CDError CD3Smoothable::AddGeometry(const char* fname)
{
  CD3List gList;
  //
  CD3ListInit(&gList);
  if (CD3ListReadGeom(&gList, fname)) {
    //
    //  AddGeometryTo runs over every point in the array asking the
    //  geometry whether the point is inside and thus should be
    //  made inactive.
    //  NOTE that it runs over it in real space and index space
    //  at the same time.
    //
    Point3D p;
    uint32_t ix, iy, iz;
    int check = 0;
    for (iz = 0, p.m[2] = mData.mMin[2]; iz < mData.mNVal[2]; iz++, p.m[2] += mData.mDelta[2]) {
      for (iy = 0, p.m[1] = mData.mMin[1]; iy < mData.mNVal[1]; iy++, p.m[1] += mData.mDelta[1]) {
        for (ix = 0, p.m[0] = mData.mMin[0]; ix < mData.mNVal[0]; ix++, p.m[0] += mData.mDelta[0]) {
          if (fabs(p.m[0]) < 0.1) {
            check++;
          }
          if (CD3ListPointIn(&gList, &p, mData.mDelta[0])) {
            mPointType((iz * mData.mNVal[1] + iy) * mData.mNVal[0] + ix) = 0;
          }
        }
      }
    }
    return kCDNoErr;
  }
  fprintf(stderr, "Cannot read geometry from file %s.\n", fname);
  return kCDBadGeom;
}

//
CDError CD3Smoothable::Smooth(int nPass)
{
  int pass, comp;
  double err = 0.0;
  uint64_t idx, idy, idz, index, rIndex;
  double newVal;
  double terr, wa, wx, wy, wz;
  double* a = mData.mField;
  /*
   *  Start by computing offsets for the neighbours. The x neigbours are
   *  always 1 higher and 1 lower in index but the y neighbours are a
   *  distance mNVal[0] (the number of possible x values) away while the
   *  z neighbours are even further.
   *  These are all multiplied by three because there are 3 components
   *  per vector.
   */
  const uint64_t dx = 3;
  const uint64_t dy = 3 * mData.mNVal[0];
  const uint64_t dz = 3 * mData.mNVal[0] * mData.mNVal[1];
  assert(NULL != a);
  //
  //  First make sure that we have a 3D leaf array.
  //
  if  (mData.mNSubField > 0) {
    fprintf(stderr, "Field is not a leaf.\n");
    return kCDNotLeaf;
  }
  //
  // Have to make sure that the data are suitable for 4-way
  // averaging.
  // First, the data must be 3D!
  //
  if (mData.mStride != 0) {
    fprintf(stderr, "Attempt to average 2D field %s.\n", mData.mFieldName);
    return kCDNot4Fold;
  }
  //
  //  Now build the point type array.
  //
  //
  //  Read in the geometry
  //  SmoothPrintOn(dp, pointType, stdout);
  //
  //  And run the smooth.
  //
  /*
   *  Weightings to compensate for non-isotropic grid.
   */
  wa = 1.0 / (1.0/(mData.mDelta[0]*mData.mDelta[0]) +
              1.0/(mData.mDelta[1]*mData.mDelta[1]) +
              1.0/(mData.mDelta[2]*mData.mDelta[2]));
  wx = wa / (2.0 *mData.mDelta[0]*mData.mDelta[0]);
  wy = wa / (2.0 *mData.mDelta[1]*mData.mDelta[1]);
  wz = wa / (2.0 *mData.mDelta[2]*mData.mDelta[2]);
  /*
   *  Now we do the fancy red-black scanning.
   *  This is made a little more complex because we have to
   *  do all three components at each point.
   */
  for (pass = 0; pass < nPass; pass++) {
    err = 0.0;
    for (idz = 0; idz < mData.mNVal[2]; idz ++) {  // Red
      for (idy = 0; idy < mData.mNVal[1]; idy ++) {  // Red
        for (idx = ((idy + idz) & 1); idx < mData.mNVal[0]; idx += 2) {  // Red
          index = idz * dz + idy * dy + idx*dx;
          rIndex = index / 3;
          switch (mPointType(rIndex)) {
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
    for (idz = 0; idz < mData.mNVal[2]; idz ++) {  // Black
      for (idy = 0; idy < mData.mNVal[1]; idy ++) {  // Black
        for (idx = ((idy + idz + 1) & 1); idx < mData.mNVal[0]; idx += 2) {  // Black
          index = idz * dz + idy * dy + idx*dx;
          rIndex = index / 3;
          switch (mPointType(rIndex)) {
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
  return kCDNoErr;
}
//
//  Write our data back to a file.
//
CDError CD3Smoothable::WriteBinary(const char* filename)
{
  FILE* ofp = fopen(filename, "wb");
  if (ofp == NULL) {
    fprintf(stderr, "Failed to open %s for writing.", filename);
    return kCDCantOpenOut;
  }
  CDError err = WriteBinary(ofp);
  fclose(ofp);
  return err;
}

CDError CD3Smoothable::WriteBinary(FILE* ofp)
{
  //
  //  Write field to file.
  //
  if (!CD3WriteBinary(&mData, ofp)) {
    fprintf(stderr, "CD3Smoothable::WriteBinary failed.\n");
    return kCDBadWrite;
  }
  return kCDNoErr;
}
