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
 *  To help understand this grid structure I do some
 *  analysis as I read the values in. I track the max
 *  and min values on each dimension.
 *
 *  To make this easier for Fred to use directly I am
 *  translating it to pure C.
 *  To support 2D and 3D specialised sub-classes in C
 *  I have to do my own function overloading, which I do
 *  by ensuring that all methods in the base class have
 *  names that begin CD.
 *
 */
#ifndef COMSOLDATA_H
#define COMSOLDATA_H

//#include <sys/cdefs.h>
#include <stdio.h>
#include <float.h>
#include <stdbool.h>

#if defined(__cplusplus)
extern "C" {
#endif

//
//  These two kludgy globals are used to pass filenames to the binary
//  writer if you wish.
//
extern const char* gFieldFileName;
extern const char* gModelFileName;

/*
 *  CDError is an enumerated type for the various
 *  errors that can be returned during the parsing
 *  of a data file.
 */
typedef enum CDErrorTag {
  kCDNoErr = 0,
  kCDCantOpenIn,
  kCDIncompleteHeader,
  kCDAllocFailed,
  kCDNameAllocFailed,
  kCDCantOpenOut,
  kCDBadStructure,    // Structure of file not match class
  kCDNotLeaf,
  kCDNot4Fold,
  kCDBadGeom,
  kCDNoArgs,
  kCDBadWrite,
  kCDXYCompatFail,
  kCDError
} CDError;

//
//  Little helper structure that collects range information
//  about a dimension. Each dimension has a max value, a min
//  value, a number of different values in the range, and the
//  resulting increment from one value to the next (delta).
//
typedef struct CDRangeTag {
  double mMin;
  double mMax;
  double mDelta;
  unsigned int mNVal;
  bool mActive;
} CDRange;

//
//  Main structure to represent the data.
//
typedef struct CDDataTag {
  /*
   *  Instance vars.
   */
  int mNDimension;          // Number of dimensions in data in file
  int mNLine;               // Number of data lines in file
  int mNExpression;         // Number of expressions in file
  int mNHeadline;           // Number of lines parsed in header
  char** mExprNames;        // Array of expression names
  double** mDStore;         // Array of arrays of data. First coords.
  CDRange* mRange;          // Array of range info for each dimension
  char* mFileName;
} CDData;

//
//  Instead of a constructor, CDRange has an initialisation function.
//
void CDRangeInit(CDRange* r);
//
//  Init fills in the data structure using the information in the file.
//
CDError CDInit(CDData* dp, const char* fname);
//
//  Finish tidies up after us, releasing our storage. Every call of
//  CDInit should be balanced by a call to CDFinish.
//
void CDFinish(CDData* dp);
//
//  Accessors.
//
double CDGetValueAtIndex(CDData* dp, unsigned int dim, unsigned int index[3]);
double CDGetValueAtPoint(CDData* dp, unsigned int dim, double coord[3]);
//
//  This writes the contents of the data store to a set of files
//  with names in the format basename_<exprname>.txt
//
CDError CDWriteBinaryTo(CDData* dp, const char* basename);
//
//  Class method WriteErrorOn writes a human readable error message to
//  a FILE.
//
void CDWriteErrorOn(FILE* ofp, CDError theErr);
//
//    Internal helpers.
//
CDError CDParseHeader(FILE* ifp, CDData* dp);
//
//  Analyse figures out the structure of the data and saves
//  the info in the mRange structures.
//
void CDAnalyse(CDData* dp);
//
#if defined(__cplusplus)
}
#endif

#endif // COMSOLDATA_H
