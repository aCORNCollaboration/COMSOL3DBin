//
//  main.c
//  COMSOL3DBin
//
//  Converts a text COMSOL data file containing a 3D grid
//  section of field into a binary format file.
//
//  COMSOL3D2Bin [-a] [-c] [-f] [-s:<geomfile.txt>] [-n:<nPass>] <textfile.txt>
//
//  will produce textfile.bin.
//  -c  Move to checking phase after build phase.
//  -a  Four-fold average (only for 3D input files)
//  -f  Process a FEMM input file rather than a
//      COMSOL file--input order is altered.
//  -n  Set number of smoothing passes (only meaningful if -s present)
//  -s  Use the geometry info to GS smooth the data.
//
//  Created by Brian Collett on 3/13/14.
//  Copyright (c) 2014 Brian Collett. All rights reserved.
//  BCollett 4/21/14 add the ability to average files with
//  four-fold symmetry to improve the noise in the data.
//  BCollett 7/24/14 Add support for FEMM input files output
//  by MATLAB or Octave.
//  BCollett 7/29/15 Add support for Gauss-Seidel smoothing controled
//  by a geometry file.
//  BCollett 7/30/15 Add an argument to set the number of smoothing
//  passes.
//

#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#include <math.h>
#include "COMSOLData3D.h"
#include "CD3List.h"

int ProcessArguments(int argc, const char** argv);
int QuadAverage(CD3Data* cd);
void DoCheck(const char* name);
int DoFile(const char* filename);

//static const int kMaxNFiles = 20;   Not sure which version of C this needs
#define kMaxNFiles 20

bool gDoAverage = false;
bool gCheckFile = false;
bool gFEMMFile = false;
int gNFile = 0;
int gNPass = 1;
const char* gGeomFilename = NULL;
const char* gFilenames[kMaxNFiles];

int main(int argc, const char * argv[])
{
  int fileNum = 0;
  const char* filename = NULL;
  CDError theErr;
  printf("COMSOL3D2Bin\n");
  //
  //  Process arguments.
  //
  theErr = ProcessArguments(argc, argv);
  if (theErr != kCDNoErr) {
    fprintf(stderr, "Processing terminated with error %d.\n", theErr);
    return theErr;
  }
  //
  //  Work through the input files.
  //
  while (fileNum < gNFile) {
    filename = gFilenames[fileNum++];
    theErr = DoFile(filename);
    fprintf(stderr, "Processing file %s terminated with error %d.\n",
            filename, theErr);
  }
  return 0;
}
//
//  Handle a single file.
//
int DoFile(const char* filename)
{
  char outName[256];
  char* ext;
  FILE* ofp;
  CD3Data cData;
  CDError theErr;
  //
  //  Read the file in.
  //
  if (gFEMMFile) {
    theErr = CD3InitFEMM(&cData, filename);
  } else {
    theErr = CD3Init(&cData, filename);
  }
  if (theErr != kCDNoErr) {
    fprintf(stderr, "Error %d: Failed to read file %s.\n",
            theErr, filename);
    return 2;
  }
  //
  //  If desired do the average.
  //  Note that this requires us to throw away the raw data and
  //  to create a new copy of the output array that will become
  //  the final copy.
  //
  if (gDoAverage) {
    if ((theErr = QuadAverage(&cData)) != kCDNoErr) {
      fprintf(stderr, "Error %d: Failed to average file %s.\n",
              theErr, filename);
      return 4;
    }
  }
  if (NULL != gGeomFilename) {
    GSSmooth(gGeomFilename, &cData, gNPass);
  }
  //
  //  Construct output file name.
  //
  strncpy(outName, filename, 254);
  ext = strrchr(outName, '.');
  if (gDoAverage) {
    strcpy(ext, "_av.bin");
  } else {
    strcpy(ext, ".bin");
  }
  ofp = fopen(outName, "wb");
  if (ofp == NULL) {
    fprintf(stderr, "Failed to open %s for writing.", outName);
    return 3;
  }
  //
  //  Write field to file.
  //
  if (!CD3WriteBinary(&cData, ofp)) {
    fprintf(stderr, "Binary write failed.\n");
    return 4;
  }
  CD3Finish(&cData);
  fclose(ofp);
  //
  //  Now if checking read it back in.
  //
  if (gCheckFile) {
    DoCheck(outName);
  }
  return 0;
}
//
//  This allows you to probe the resulting file.
//
void DoCheck(const char* name)
{
  CD3Data cData;
  double x,y,z;
  double coord[3];
  double field[3];
  FILE* ifp = fopen(name, "rb");
  if (NULL == ifp) {
    fprintf(stderr, "Failed to open %s for reading.\n", name);
    return;
  }
  if (CD3ReadBinary(&cData, ifp)) {
    for (; ; ) {
      printf("Enter 3 coord values (x<-100 to stop).\n");
      scanf("%lf %lf %lf", &x, &y, &z);
      if (x < -100.0) {
        break;
      }
      coord[0] = x;
      coord[1] = y;
      coord[2] = z;
      if (CD3GetEAtPoint(&cData, coord, field)) {
        printf("<%lf,%lf,%lf> -> [%lf,%lf,%lf]\n",
               x,y,z,field[0],field[1],field[2]);
      } else {
        printf("Point <%lf,%lf,%lf> out of bounds.\n",
               x,y,z);
      }
    }
  }
}
 //
 //  This does the averaging.
 //  It first has to build a new data store.
 //
 int QuadAverage(CD3Data* dp)
{
  uint32_t i, j, k;   // x, y, z indices
  uint32_t jmid, imid;
  double eps;
  //
  //  First just a quick check that we are sane. This must be a leaf
  //  field.
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
  // The x and y ranges must match
  // and the ranges must be centered on zero.
  // These can be inexact comparisons.
  //
  eps = 1.0e-6;
  if (fabs(dp->mMin[0] - dp->mMin[1]) > eps) {
    fprintf(stderr, "dp->mMin[0] %12.8g != dp->mMin[1] %12.8g\n",
            dp->mMin[0], dp->mMin[1]);
    return kCDNot4Fold;
  }
  if (fabs(dp->mMin[0] + dp->mMax[0]) > eps) {  // Yes PLUS, opposite signs
    fprintf(stderr, "dp->mMin[0] %g != dp->mMax[0] %g\n",
            dp->mMin[0], dp->mMax[0]);
    return kCDNot4Fold;
  }
  if (fabs(dp->mMax[0] - dp->mMax[1]) > eps) {
    fprintf(stderr, "dp->mMax[0] %g != dp->mMax[1] %g\n",
            dp->mMax[0], dp->mMax[1]);
    return kCDNot4Fold;
  }
  if (fabs(dp->mMin[1] + dp->mMax[1]) > eps) { // Ditto
    fprintf(stderr, "dp->mMin[1] %g != dp->mMax[1] %g\n",
            dp->mMin[1], dp->mMax[1]);
    return kCDNot4Fold;
  }
  //
  //  x and y index ranges have separate mid-points.
  //  Location depends on odd or even.
  //
  imid = dp->mNVal[0];
  imid = ((imid % 1) == 0) ? imid/2 : (imid-1)/2;
  jmid = dp->mNVal[1];
  jmid = ((jmid % 1) == 0) ? jmid/2 : (jmid-1)/2;
  //
  //  Then run over the arrays copying the averages into place.
  //
  for (k = 0; k < dp->mNVal[2]; k++) {
    uint32_t idxk = k*dp->mNVal[1];
    for (j = jmid;  j < dp->mNVal[1]; j++) {
      uint32_t jn = (dp->mNVal[1] - 1) - j;
      uint32_t idxkjp = (idxk + j)*dp->mNVal[0];
      uint32_t idxkjn = (idxk + jn)*dp->mNVal[0];
      for (i = imid;  i < dp->mNVal[0]; i++) {
        //
        //  Four indices for the positive and negative versions of j amd i
        //
        uint32_t in = (dp->mNVal[0] - 1) - i;
        uint32_t idxpp = (idxkjp + i)*3;
        uint32_t idxpn = (idxkjp + in) * 3;
        uint32_t idxnp = (idxkjn + i)*3;
        uint32_t idxnn = (idxkjn + in) * 3;
        //
        //  x components.
        //
        double av = 0.25 * (dp->mField[idxpp+0] + dp->mField[idxnp+0] -
                            dp->mField[idxpn+0] - dp->mField[idxnn+0]);
        dp->mField[idxpp+0] = dp->mField[idxnp+0] = av;
        dp->mField[idxpn+0] = dp->mField[idxnn+0] = -av;
        //
        //  y components.
        //
        av = 0.25 * (dp->mField[idxpp+1] + dp->mField[idxpn+1] -
                     dp->mField[idxnp+1] - dp->mField[idxnn+1]);
        dp->mField[idxpp+1] = dp->mField[idxpn+1] = av;
        dp->mField[idxnp+1] = dp->mField[idxnn+1] = -av;
        //
        //  z components.
        //
        av = 0.25 * (dp->mField[idxpp+2] + dp->mField[idxpn+2] +
                     dp->mField[idxnp+2] + dp->mField[idxnn+2]);
        dp->mField[idxpp+2] = dp->mField[idxnp+2] = dp->mField[idxpn+2] =
        dp->mField[idxnn+2] = av;
      }
    }
  }
  return kCDNoErr;
}
//
//  ProcessArguments
//  Interprets anything beginning with '-' as an option and
//  tries to use it. Anything else is collected in the file list.
//
int ProcessArguments(int argc, const char** argv)
{
  int iVal, argn;
  if (argc < 2) {
    fprintf(stderr, "No arguments given.\n");
    return 1;
  }
  for (argn = 1; argn < argc; argn++) {
    if (argv[argn][0] == '-') {
      switch (argv[argn][1]) {
        case 'c':
          gCheckFile = true;
          break;

        case 'a':
          gDoAverage = true;
          break;

        case 'f':
          gFEMMFile = true;
          break;

        case 'n':
          if (argv[argn][2] == ':') {
            if (sscanf(&argv[argn][3], "%d", &iVal) == 1) {
              gNPass = iVal;
            } else {
              fprintf(stderr, "Failed to find valid number of passes in argument %s\n", argv[argn]);
            }
          }
          break;

        case 's':
          if (argv[argn][2] == ':') {
            gGeomFilename = &argv[argn][3];
          }
          break;

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
  return 0;
}
/*
 //
 //  This does the averaging.
 //  It computes the four-way average by doing two successive
 //  two-way averages, one in x and the other in y. This way
 //  the file can have different grid spacings in the two directions.
 //
 int QuadAverage(CD3Data* dp)
 {
 uint32_t i, j, k;   // x, y, z indices
 uint32_t mid;
 //
 //  First just a quick check that we are sane. This must be a leaf
 //  field.
 //
 if  (dp->mNSubField > 0) {
 fprintf(stderr, "Field is not a leaf.\n");
 return kCDNotLeaf;
 }
 //
 // Have to make sure that the data are suitable for 4-way
 // averaging. The x and y ranges must be centered on zero.
 //
 if (dp->mMin[0] != -dp->mMax[0]) {
 fprintf(stderr, "dp->mMin[0] %g != dp->mMax[0] %g\n",
 dp->mMin[0], dp->mMax[0]);
 return kCDNot4Fold;
 }
 if (dp->mMin[1] != -dp->mMax[1]) {
 fprintf(stderr, "dp->mMin[1] %g != dp->mMax[1] %g\n",
 dp->mMin[1], dp->mMax[1]);
 return kCDNot4Fold;
 }
 //
 //  First average over y. Note that I only output HALF of
 //  the data here. The other half are redundant.
 //
 //  Location depends on odd or even.
 //
 mid = dp->mNVal[0];
 mid = ((mid % 1) == 0) ? mid/2 : (mid-1)/2;
 //
 //  Then run over the arrays copying the averages into place.
 //
 for (k = 0; k < dp->mNVal[2]; k++) {
 uint32_t idxk = k*dp->mNVal[1];
 for (j = mid;  j < dp->mNVal[1]; j++) {
 uint32_t jn = (dp->mNVal[1] - 1) - j;
 uint32_t idxkjp = (idxk + j)*dp->mNVal[0];
 uint32_t idxkjn = (idxk + jn)*dp->mNVal[0];
 for (i = mid;  i < dp->mNVal[0]; i++) {
 //
 //  Four indices for the positive and negative versions of j amd i
 //
 uint32_t in = (dp->mNVal[0] - 1) - i;
 uint32_t idxpp = (idxkjp + i)*3;
 uint32_t idxpn = (idxkjp + in) * 3;
 uint32_t idxnp = (idxkjn + i)*3;
 uint32_t idxnn = (idxkjn + in) * 3;
 //
 //  x components.
 //
 double av = 0.25 * (dp->mField[idxpp+0] + dp->mField[idxnp+0] -
 dp->mField[idxpn+0] - dp->mField[idxnn+0]);
 dp->mField[idxpp+0] = dp->mField[idxnp+0] = av;
 dp->mField[idxpn+0] = dp->mField[idxnn+0] = -av;
 //
 //  y components.
 //
 av = 0.25 * (dp->mField[idxpp+1] + dp->mField[idxpn+1] -
 dp->mField[idxnp+1] - dp->mField[idxnn+1]);
 dp->mField[idxpp+1] = dp->mField[idxpn+1] = av;
 dp->mField[idxnp+1] = dp->mField[idxnn+1] = -av;
 //
 //  z components.
 //
 av = 0.25 * (dp->mField[idxpp+2] + dp->mField[idxpn+2] +
 dp->mField[idxnp+2] + dp->mField[idxnn+2]);
 dp->mField[idxpp+2] = dp->mField[idxnp+2] = dp->mField[idxpn+2] =
 dp->mField[idxnn+2] = av;
 }
 }
 }
 return kCDNoErr;
 }

 Here is the original QuadAverage code. It needs only a single pass
 over the data but requires that the grids be identical in the x and y
 directions.
 */
