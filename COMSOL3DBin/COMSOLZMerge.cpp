//
//  COOMSOLZMerge.cpp
//  COMSOLZMerge
//
//  Driver for the COMSOLZMerge tool that combines two .bin files
//  into one so long as they are z-adjacent.
//  It requires that the two files have the same x and y limits and
//  that all three deltas match. It can deal with overlap
//  between the files.
//
//  COMSOLZMerge <binfile1.bin> <binfile2.bin> ...
//
//  will produce binfile_ZMin-ZMax.bin.
//
//  Created by Brian Collett on 8/5/15.
//  Copyright (c) 2015 Brian Collett. All rights reserved.
//
#include <inttypes.h>
#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#include <math.h>
#include "COMSOLData3D.h"
#include "CD3List.h"

CDError ProcessArguments(int argc, const char** argv);
CDError MergeData(const CD3Data* ind1, const CD3Data* ind2, CD3Data* outd);
const CD3Data* XYCompatible(const CD3Data* ind1, const CD3Data* ind2);
bool NearlyEqual(double v1, double v2);
#define MIN(x,y) (((x) < (y)) ? (x) : (y))

static const int kMaxNFiles = 2;  // Not sure which version of C this needs
//#define kMaxNFiles 2

int gNFile = 0;
const char* gFilenames[kMaxNFiles];
char gOutName[256];

int main(int argc, const char * argv[])
{
  CDError theErr;
  CD3Data cData1;
  CD3Data cData2;
  CD3Data cData3;
  printf("COMSOLZMerge\n");
  //
  //  Process arguments.
  //
  theErr = ProcessArguments(argc, argv);
  if (theErr != kCDNoErr) {
    fprintf(stderr, "Processing terminated with error %d.\n", theErr);
    return theErr;
  }
  //
  //  Read in the input files.
  //
  FILE* ifp = fopen(gFilenames[0], "rb");
  if (nullptr == ifp) {
    fprintf(stderr, "Failed to open file %s for reading.\n", gFilenames[0]);
    return kCDCantOpenIn;
  }
  if (!CD3ReadBinary(&cData1, ifp)) {
    fprintf(stderr, "Error %d: Failed to read file %s.\n",
            theErr, gFilenames[0]);
    fclose(ifp);
    return kCDCantOpenIn;
  }
  fclose(ifp);
  cData1.mFieldName = gFilenames[0];
  printf("Loaded file %s.\n", gFilenames[0]);
  //
  ifp = fopen(gFilenames[1], "rb");
  if (nullptr == ifp) {
    fprintf(stderr, "Failed to open file %s for reading.\n", gFilenames[1]);
    return kCDCantOpenIn;
  }
  if (!CD3ReadBinary(&cData2, ifp)) {
    fprintf(stderr, "Error %d: Failed to read file %s.\n",
            theErr, gFilenames[1]);
    fclose(ifp);
    return kCDCantOpenIn;
  }
  fclose(ifp);
  cData2.mFieldName = gFilenames[1];
  printf("Loaded file %s.\n", gFilenames[1]);
  //
  //  Merge data into empty CData.
  //
  if ((cData1.mType != kCD3Data3) || (cData2.mType != kCD3Data3)) {
    fprintf(stderr, "Data to merge must be fully 3D.\n");
    return kCDError;
  }
  theErr = MergeData(&cData1, &cData2, &cData3);
  if (theErr != kCDNoErr) {
    fprintf(stderr, "Error %d: Unable to merge files.\n",
            theErr);
    return theErr;
  }
  printf("Files merged successfully.\n");
  printf("Output name assigned %s.\n", cData3.mFieldName);
  //
  //  Write results.
  //
  FILE* ofp = fopen(cData3.mFieldName, "wb");
  if (nullptr == ofp) {
    fprintf(stderr, "Error %d: Unable to open file %s for writing.\n",
            theErr, cData3.mFieldName);
    return kCDBadWrite;
  }
  if (CD3WriteBinary(&cData3, ofp)) {
    printf("Merged data written to %s.\n", cData3.mFieldName);
  }
  return 0;
}
//
//  Try to merge two CD3Data into a third.
//
CDError MergeData(const CD3Data* ind1, const CD3Data* ind2, CD3Data* outd)
{
  //
  //  Make sure inputs are compatible.
  //
  const CD3Data* lowz = XYCompatible(ind1, ind2);
  if (nullptr == lowz) {
    fprintf(stderr, "Merge Error: Input files incompatible.\n");
    return kCDXYCompatFail;
  }
  //
  //  If we get here then the two files are compatible.
  //  Now we need to figure out how many planes we take from
  //  the low file and how many from the high.
  //  Let's make up a rule that we always take all of the
  //  lowz file and then some of the hiz.
  //
  const CD3Data* hiz = (lowz == ind1) ? ind2 : ind1;
  unsigned int nz = lowz->mNVal[2];
  unsigned int kmin = (lowz->mMax[2] - hiz->mMin[2]) / lowz->mDelta[2];
  nz += hiz->mNVal[2] - kmin;
  //
  //  Start filling in out. We know the ranges, deltas, and numbers
  //  of values.Copy from lowz and then adjust the z vals.
  //
  outd->mType = kCD3Data3;
  for (int i = 0; i < 3; i++) {
    outd->mMin[i] = lowz->mMin[i];
    outd->mDelta[i] = lowz->mDelta[i];
    outd->mMax[i] = lowz->mMax[i];
    outd->mNVal[i] = lowz->mNVal[i];
  }
  outd->mMax[2] = hiz->mMax[2];
  outd->mNVal[2] = nz;
  outd->mNSubField = 0;
  //
  //  Get storage for data and copy.
  //
  uint64_t size = outd->mNVal[0];
  size *= outd->mNVal[1];
  size *= outd->mNVal[2];     // Avoids overflow
  size *= 3;    // 3 doubles at each point!
  outd->mField = (double *) malloc(size * sizeof(double));
  if (nullptr == outd->mField) {
    fprintf(stderr, "Cannot allocate %" PRIu64 " bytes of space for output.\n", size);
    return kCDAllocFailed;
  }
  uint64_t lowsize = lowz->mNVal[0];
  lowsize *= lowz->mNVal[1];
  lowsize *= lowz->mNVal[2];     // Avoids overflow
  lowsize *= 3;    // 3 doubles at each point!
  memcpy(outd->mField, lowz->mField, lowsize);
  uint64_t hisize = hiz->mNVal[0];
  hisize *= hiz->mNVal[1];
  hisize *= (hiz->mNVal[2] - kmin);
  hisize *= 3;    // 3 doubles at each point!
  uint64_t histart = hiz->mNVal[0];
  histart *= hiz->mNVal[1];
  histart *= kmin;
  histart *= 3;    // 3 doubles at each point!
  memcpy(outd->mField+lowsize, hiz->mField + histart, hisize * sizeof(double));
  //
  //  Make up a name.
  //
  strcpy(gOutName, lowz->mFieldName);
  char* baseEnd = strrchr(gOutName, '.');
  sprintf(baseEnd, "%5.2f-%5.2f.bin",lowz->mMin[2], hiz->mMax[2]);
  outd->mFieldName = gOutName;
  return kCDNoErr;
}
//
//  Make sure data are compatible.
//  This means that x and y ranges and deltas match.
//  Returns pointer to lower z data or null if fail.
//
const CD3Data* XYCompatible(const CD3Data* ind1, const CD3Data* ind2)
{
  //
  //  First make sure that x and y ranges and all deltas match.
  //
  for (int i = 0; i < 2; i++) {
    if (!NearlyEqual(ind1->mMin[i], ind2->mMin[i])) {
      fprintf(stderr, "%c mins (%f and %f) don't match.\n",
              (i == 0) ? 'x' : 'y', ind1->mMin[i], ind2->mMin[i]);
      return nullptr;
    }
    if (!NearlyEqual(ind1->mMax[i], ind2->mMax[i])) {
      fprintf(stderr, "%c maxs (%f and %f) don't match.\n",
              (i == 0) ? 'x' : 'y', ind1->mMax[i], ind2->mMax[i]);
      return nullptr;
    }
    if (!NearlyEqual(ind1->mDelta[i], ind2->mDelta[i])) {
      fprintf(stderr, "%c deltas (%f and %f) don't match.\n",
              (i == 0) ? 'x' : 'y', ind1->mDelta[i], ind2->mDelta[i]);
      return nullptr;
    }
  }
  if (!NearlyEqual(ind1->mDelta[2], ind2->mDelta[2])) {
    fprintf(stderr, "z deltas (%f and %f) don't match.\n",
            ind1->mDelta[2], ind2->mDelta[2]);
    return nullptr;
  }
  //
  //  Then make sure that z coords form a touching or overlapping
  //  sequence.
  //  First figure out which has lower z.
  //
  const CD3Data* lowz = ind1;
  const CD3Data* hiz = ind2;
  if (ind2->mMin[2] <= ind1->mMin[2]) {
    lowz = ind2;
    hiz = ind1;
  }
  if (lowz->mMax[2] < hiz->mMin[2]) {
    fprintf(stderr, "z coords (%f,%f) and (%f,%f) neither touch nor overlap.\n",
            lowz->mMin[2], lowz->mMax[2], hiz->mMin[2], hiz->mMax[2]);
    return nullptr;
  }
  return lowz;
}
//
//  This is a weak comparison between two floating point numbers. It
//  treats two numbers as equal so long as they don't differ by more
//  than some fraction.
//
static const double fract = 1e-6;
bool NearlyEqual(double v1, double v2)
{
  double tol = fract * MIN(fabs(v1), fabs(v2));
  bool res = fabs(v1 - v2) < tol;
  return fabs(v1 - v2) < tol;
}

//
//  ProcessArguments
//  Interprets anything beginning with '-' as an option and
//  tries to use it. Anything else is collected in the file list.
//
CDError ProcessArguments(int argc, const char** argv)
{
  int argn;
  if (argc < 2) {
    fprintf(stderr, "No arguments given.\n");
    return kCDNoArgs;
  }
  for (argn = 1; argn < argc; argn++) {
    if (argv[argn][0] == '-') {
      switch (argv[argn][1]) {
          
          
        default:
          fprintf(stderr, "Ignored unknown option %s.\n", argv[argn]);
          break;
      }
    } else {
      if (gNFile < kMaxNFiles) {
        gFilenames[gNFile++] = argv[argn];
      } else {
        fprintf(stderr, "Too many input files. Ignoring %s.\n", argv[argn]);
      }
    }
  }
  return kCDNoErr;
}
