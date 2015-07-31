/*
 *  COMSOLData3D.cpp
 *  COMSOLReader2
 *
 *  This is a specialisation of the COMSOLReader class that only
 *  supports 3D data files containing 3D fields. It checks that
 *  the file not merely has three dimensions but that all three
 *  spatial dimensions are active, that is have more than one value
 *  of that coordinate represented. This means that GetValueAtPoint
 *  can be somewhat simpler.
 *
 *  Created by Brian Collett on 2/22/14.
 *  Copyright (c) 2014 Brian Collett. All rights reserved.
 *  BCollett 3/12/14 Add support for binary storage of 3D files.
 *  They have a header containing information on the structure and
 *  then a single array of doubles as the data organized
 *  data[((((index[2]*dp->mNVal[1])+index[1])*dp->mNVal[1])+index[0])+comp]
 *  where comp is 0 for x, 1 for y, and 2 for z. Restructured mField to match.
 *  BCollett 3/14/14 Merge with 2D version. We keep a flag that tells us
 *  how our internal data are stored and and that allows us to forward
 *  field requests to the right code path.
 *  NOTE that the two types of field now supported are fully 3D fields in
 *  boxes aligned with the Cartesian axes and 2D fields axisymmetric about
 *  the z axis only.
 *  BCollett7/29/15 Add some helper methods to clip points to the bounds
 *  of the field and to map indices to coords and vice-versa.
 */

#ifndef __COMSOLData3D__
#define __COMSOLData3D__

#include <stdint.h>
#include "COMSOLData.h"

//
//  This distinguishes between our internal formats.
//
typedef enum CD3TypeTag {
  kCD3Data2 = 0,        // 2D axisymmetric data
  kCD3Data3,            // Full 3D data
  kCD3Unused,           // Actual data part unused, only bounds valid
  kCD3Error             // Oh Dear!
} CD3TypeTag;

//
//  Because we understand the structure of this kind of file much
//  better we can have a streamlined class. For example, we no longer
//  need to store the coordinate data.
//
//  How many sub-fields can we have?
//
#define kNSub 20
//
typedef struct CD3DataTag {
  //
  //  Three sets of data arrays, one for each dimension of either the
  //  field or the coordinate space.
  //
  CD3TypeTag mType;                     // Tells us how our data are stored
  unsigned int mNVal[3];                // Number of different indices this dim
  double mMin[3];                       // and corresponding minima
  double mMax[3];                       // Max values of each coord.
  double mDelta[3];                     // Deltas needed for coord conversion
  int mStride;                          // Used only for 2D data.
  int mNSubField;                       // Number of subfields
  const struct CD3DataTag* mSubField[kNSub];  // Stored here
  double* mField;                       // Field data
  const char* mFieldName;
} CD3Data;
//
//  Have a second structure that we use as the header for a binary file.
//  It incorporates extras such as a data offset, space for a filename,
//  and, of course, a magic number.
//
extern uint32_t gCD3Magic;
//
typedef struct CD3HeadTag {
  uint32_t magic;
  uint32_t dataOffset;
  char modelName[64];
  char fileName[64];
  CD3Data dp;
  char filler[0];           // On disk will be stored as 256 bytes.
} CD3Header;

#if defined(__cplusplus)
extern "C" {
#endif
//
//  Init fills in the data structure using the information in a
//  text file produced by COMSOL.
//
CDError CD3Init(CD3Data* dp, const char* fname);
//
//  InitFEMM does the same for a file extracted from FEMM
//  and written by Matlab or Octave.
//
CDError CD3InitFEMM(CD3Data* dp, const char* fname);
//
//  Finish tidies up after us, releasing our storage. Every call of
//  CD3Init should be balanced by a call to CD3Finish.
//
void CD3Finish(CD3Data* dp);
//
//  File operations.
//  CD3ReadBinary fills in the data structure with info from binary file.
//
bool CD3ReadBinary(CD3Data* dp, FILE* ifp);
//
//  CD3WriteBinary saves a complete data structure in a binary file
//  ready to be read in by CD3ReadBinary.
//
bool CD3WriteBinary(CD3Data* dp, FILE* ofp);
//
//  Accessors.
//  First checks whether a point is inside this field.
//
bool PtInBounds(const CD3Data* dp, const double coord[3]);
//
//  This understands that there are two internal formats and
//  redirects the requests to the correct internal helpers.
//
bool CD3GetEAtPoint(const CD3Data* dp, const double coord[3], double* EField);
//
//  This can tell you what file the data for a particular point
//  came from.
//
const char* CD3GetNameAtPoint(const CD3Data* dp, const double coord[3]);
/*
 *  I am not currently using these.
 *
 double CD3GetValueAtIndex(const CD3Data* dp, unsigned int dim,
 const unsigned int index[3]);
 double CD3GetValueAtPoint(const CD3Data* dp, unsigned int dim,
 const double coord[3]);
 double CD3GetExAtPoint(const CD3Data* dp, const double coord[3]);
 double CD3GetEyAtPoint(const CD3Data* dp, const double coord[3]);
 double CD3GetEzAtPoint(const CD3Data* dp, const double coord[3]);
 */
//
//  Force a point to the nearest point inside the bounds.
//
void CD3ClipPt(const CD3Data* dp, const double coord[3], double newCoord[3]);

//
//  Map checks whether a point is in range and if it is fills in the
//  indices of the corresponding entry and returns true, else false.
//
bool CD3Map(const CD3Data* dp, const double coord[3], uint32_t newIndices[3]);

//
//  This translates an index trio into a single index.
//
uint64_t CD3IndexAt(const CD3Data* dp, uint32_t ix, uint32_t iy, uint32_t iz);

#if defined(__cplusplus)
}
#endif


#endif /* defined(__COMSOLReader2__COMSOLData3D__) */
