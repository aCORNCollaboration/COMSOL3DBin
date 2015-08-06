//
//  COMSOLSmooth.cpp
//  COMSOLSmooth
//
//  Gauss-Seidel smooth a COMSOL binary file optionally with
//  geometry control.
//  Without a geometry file this sets all the points on the boundary
//  to Dirichlet and smooths all the interior points as many times
//  as desired.
//  With a Geometry file additional interior points can be set to
//  Dirichlet bounds and excluded from the smoothing.
//
//  COMSOLSmooth [-g:<geomfile.txt>] [-n:<nPass>] <comsolfile.bin> ...
//
//  will produce comsolfile_sm.bin.
//  -g  Use the geometry info to GS smooth the data.
//  -n  Set number of smoothing passes (only meaningful if -g present)
//  You can pass multiple files and all will be processed the same way.
//
//  Created by Brian Collett on 8/5/15.
//  Copyright (c) 2015 Brian Collett. All rights reserved.
//
#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#include <math.h>
#include "COMSOLData3D.h"
#include "CD3Smoothable.h"
#include "CD3List.h"

CDError ProcessArguments(int argc, const char** argv);
CDError DoFile(const char* filename);

static const int kMaxNFiles = 20;   // Fine in C++

int gNFile = 0;
int gNPass = 1;
const char* gGeomFilename = NULL;
const char* gFilenames[kMaxNFiles];

int main(int argc, const char * argv[])
{
  int fileNum = 0;
  const char* filename = NULL;
  CDError theErr;
  printf("COMSOLSmooth\n");
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
CDError DoFile(const char* filename)
{
  char outName[256];
  char* ext;
  CD3Smoothable cData;
  CDError theErr;
  //
  //  Read the file in.
  //
  theErr = cData.ReadBinary(filename);
  if (kCDNoErr != theErr) {
    fprintf(stderr, "Error %d: Failed to read file %s.\n",
            theErr, filename);
    return kCDCantOpenIn;
  }
  printf("Loaded file %s.\n", filename);
  //
  //  Smooth the file.
  //  This is done in-place.
  //
  if (nullptr != gGeomFilename) {
    theErr = cData.AddGeometry(gGeomFilename);
    if (kCDNoErr != theErr) {
      fprintf(stderr, "Error %d: Failed to read geometry file %s.\n",
              theErr, gGeomFilename);
      return kCDCantOpenIn;
    }
    printf("Loaded geometry file %s.\n", gGeomFilename);
  }
  cData.Smooth(gNPass);
  //
  //  Construct output file name.
  //
  strncpy(outName, filename, 254);
  ext = strrchr(outName, '.');
  strcpy(ext, "_sm.bin");
  theErr = cData.WriteBinary(outName);
  return theErr;
}
//
//  ProcessArguments
//  Interprets anything beginning with '-' as an option and
//  tries to use it. Anything else is collected in the file list.
//
CDError ProcessArguments(int argc, const char** argv)
{
  int iVal, argn;
  if (argc < 2) {
    fprintf(stderr, "No arguments given.\n");
    return kCDNoArgs;
  }
  for (argn = 1; argn < argc; argn++) {
    if (argv[argn][0] == '-') {
      switch (argv[argn][1]) {
        case 'g':
          if (argv[argn][2] == ':') {
            gGeomFilename = &argv[argn][3];
          }
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
