//
//  COOMSOLTxt2Bin.cpp
//  COMSOLTxt2Bin
//
//  Driver for the simplified COMSOLTxt2Bin tool that ONLY converts
//  COMSOL .txt or FEMM .txt files into my .bin files.
//
//  COMSOL3D2Bin [-f] <textfile1.txt> <textfile2.tx> ...
//
//  will produce textfile.bin.
//  -f  Process a FEMM input file rather than a
//      COMSOL file--input order is altered.
//
//  Created by Brian Collett on 8/4/15.
//  Copyright (c) 2015 Brian Collett. All rights reserved.
//

#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#include <math.h>
#include "COMSOLData3D.h"
#include "CD3List.h"

CDError ProcessArguments(int argc, const char** argv);
CDError DoFile(const char* filename);

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
  printf("COMSOLTxt2Bin\n");
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
    return kCDCantOpenIn;
  }
  //
  //  Construct output file name.
  //
  strncpy(outName, filename, 254);
  ext = strrchr(outName, '.');
  strcpy(ext, ".bin");
  ofp = fopen(outName, "wb");
  if (ofp == NULL) {
    fprintf(stderr, "Failed to open %s for writing.", outName);
    return kCDCantOpenOut;
  }
  //
  //  Write field to file.
  //
  if (!CD3WriteBinary(&cData, ofp)) {
    fprintf(stderr, "Binary write failed.\n");
    return kCDBadWrite;
  }
  CD3Finish(&cData);
  fclose(ofp);
  return kCDNoErr;
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
        case 'f':
          gFEMMFile = true;
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
