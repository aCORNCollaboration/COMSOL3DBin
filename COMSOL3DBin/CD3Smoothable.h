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

#ifndef __COMSOLSmooth__CD3Smoothable__
#define __COMSOLSmooth__CD3Smoothable__

#include <stdio.h>
#include <stdint.h>
#include "COMSOLData3D.h"

//
//  Helper class.
//
class PointArray {
protected:
  //
  //  An array of little flags arranged as a 3D array.
  //
  uint64_t mNPoint;
  uint64_t mNVal[3];    // Index limits in the three directions
  uint8_t* mArray;
public:
  PointArray();         // SUpports 2-step construction
  PointArray(uint64_t nx, uint64_t ny, uint64_t nz);
  ~PointArray();
  //
  //  SetSize either from explicit size or from a CD3Data.
  //
  CDError SetSize(uint64_t nx, uint64_t ny, uint64_t nz);
  CDError SetSize(const CD3Data& d);
  //
  //  Accessors.
  //
  uint8_t& operator()(uint64_t idx0);
  uint8_t& operator()(uint64_t idx1, uint64_t idx2, uint64_t idx3);
};

class CD3Smoothable {
protected:
  //
  //  Class var is just a default blank name.
  //
  static const char* gDefName;
  //
  //  Instance vars.
  //  We start with a CD3Data and then add our pointArray
  //
  CD3Data mData;
  PointArray mPointType;
public:
  CD3Smoothable();
  virtual ~CD3Smoothable();
  //
  //  2-step construction for error handling.
  //
  CDError ReadBinary(const char* filename);
  //
  //  Can add geometry info, possibly several.
  //
  CDError AddGeometry(const char* filename);
  //
  //  Smooth does the real work
  //
  CDError Smooth(int nPass);
  //
  //  Write our data back to a file.
  //
  CDError WriteBinary(const char* filename);
  CDError WriteBinary(FILE* ofp);
};

#endif /* defined(__COMSOLSmooth__CD3Smoothable__) */
