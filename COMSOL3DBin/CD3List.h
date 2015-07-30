//
//  CD3List.h
//  Smooth3D
//
//  This is the pure C version of the CD3List.
//
//  A CD3List is a way of putting geometry information into a Scalar3D
//  file. It reads from a text file a simple description of the geometry
//  of those regions in the Scalar3D array that constitute boundaries,
//  initially only Dirichlet boundaries.
//  The systax supported is a simple one with lines that begin with a
//  command name followed by a set of n numbers that the command can
//  interpret as it pleases.
//  The set of supported commands is
//
//  icyl <direction> <min> <max> <radius> <potential>
//
//  The first line in the file should begin BCGeom.
//  Lines that begin with # are comments and completely ignored.
//
//  The file is compiled into an internal representation (a simple list)
//  I have forced 2-step construction because reading from a file can fail.
//
//  Created by Brian Collett on 7/24/15.
//  Copyright (c) 2015 Brian Collett. All rights reserved.
//

#ifndef __CD3List__
#define __CD3List__

#include "GSSmooth.h"
#include "COMSOLData3D.h"
#include "Geometries.h"


typedef struct CD3ListTag {
  //
  //  The list of geometries.
  //
  Geom* mGList;
  //
  //  Temps used for reading.
  //
  char* mLineBuff;
  double mArg[20];
  int mNArg;
} CD3List;

//
//
//	Constructor and destructor.
//
void CD3ListInit(CD3List* l);
void CD3ListFinish(CD3List* l);
//
//  Load from a file.
//
bool CD3ListReadGeom(CD3List* l, const char* inFilename);
//
//  Make printable
//
void CD3ListPrintOn(CD3List* l, FILE* ofp);
//
//  Write into a Smoothable.
//
bool CD3ListAddGeomTo(CD3List* l, uint8_t* type, CD3Data* d);
//
//


#endif /* defined(__CD3List__) */
