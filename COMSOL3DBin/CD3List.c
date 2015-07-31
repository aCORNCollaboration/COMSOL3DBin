//
//  CD3List.cpp
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
#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#include "assert.h"
#include "GSSmooth.h"
#include "COMSOLData3D.h"
#include "CD3List.h"
//
//	Helper functions used internally.
//
//	GetArgList scans the rest of the current line for numbers which
//	it puts into the argument array.
//
static int GetArgList(CD3List* l);
//
//  BuildICylinder wants six arguments to construct an interior cylinder
//
static bool BuildICylinder(CD3List* l);

//
//  CD3List
//	Constructor and destructor.
//
void CD3ListInit(CD3List* l)
{
  l->mGList = NULL;
  l->mLineBuff = (char*) malloc(512*sizeof(char));
  l->mNArg = 0;
}

void CD3ListFinish(CD3List* l)
{
  if (NULL != l->mGList) {
    free(l->mGList);
  }
  free(l->mLineBuff);
}

//
//  Load geometry from a file.
//
bool CD3ListReadGeom(CD3List* l, const char* inFilename)
{
  char* cmd;
  /*
   *  See if we can open the file.
   */
  FILE* ifp = fopen(inFilename, "rt");
  if  (NULL == ifp) {
    fprintf(stderr, "Cannot open geometry file %s.\n",
            inFilename);
    return false;
  }
  /*
   *  Is it a geometry file?
   */
  fgets(l->mLineBuff, 511, ifp);
//  fprintf(stderr, "read head into mLineBuff alloc at %p\n", mLineBuff);
  if (strncmp(l->mLineBuff, "BCGeom", 6) != 0) {
    fprintf(stderr, "File %s is not a geometry file.\n",
            inFilename);
    fclose(ifp);
    return false;
  }
  /*
   *  Read lines and do the trivial parsing.
   */
  while (fgets(l->mLineBuff, 511, ifp) != NULL) {
//    fprintf(stderr, "read line into mLineBuff alloc at %p\n", mLineBuff);
    if (l->mLineBuff[0] != '#') {
      //
      //  Line should begin with a command.
      //
      cmd = strtok(l->mLineBuff, " ,\t\r\n");
      GetArgList(l);
      if (strncmp(cmd, "icyl", 4) == 0) {
        if (l->mNArg != 8) {
          fprintf(stderr, "Expecting 8 arguments for an icyl, found %d.",
                  l->mNArg);
          fclose(ifp);
          return false;
        }
        BuildICylinder(l);
      } else {
        fprintf(stderr, "Ignored nknown command %s.\n", cmd);
      }
    }
  }
  fclose(ifp);
  return true;
}

//
//	Helper functions used internally.
//
//	GetArgList scans the rest of the current line for numbers which
//	it puts into the argument array.
//
int GetArgList(CD3List* l) {
  char* arg;
  l->mNArg = 0;
  while ((arg = strtok(NULL, " ,\t\r\n")) != NULL) {
    if (sscanf(arg, "%lf", &((l->mArg)[l->mNArg])) == 1) {
      ++(l->mNArg);
    } else {
      fprintf(stderr, "Unknown argument %s at argument #%d.\n",
              arg, l->mNArg);
      break;
    }
    if (l->mNArg > 18) {
      fprintf(stderr, "Too many arguments in %s.\n",
              l->mLineBuff);
      break;
    }
  }
  return l->mNArg;
}

//
//  BuildICylinder wants eight arguments to construct an interior cylinder
//  icyl <xmin> <xmax> <ymin> <ymax> <zmin> <zmax> <radius>
//
bool BuildICylinder(CD3List* l)
{
  int i, act = -1, nActive = 0;
  double diff[3];
  bool active[3];
  static char sNames[4] = "xyz";
  Geom* newCyl = NULL;
  /*
   *  Verify the arguments.
   *  Make sure that min and max are ordered or same.
   */
  for (i = 0; i < 3; i++) {
    diff[i] = l->mArg[i+3] - l->mArg[i];
    if (diff[i] < 0.0) {
      fprintf(stderr, "Error building ICylinder %cMin=%lf > %cMax=%lf\n",
              sNames[i],l->mArg[i],sNames[i],l->mArg[i+3]);
      return false;
    } else if (diff[i] == 0) {
      active[i] = false;
    } else {
      active[i] = true;
      ++nActive;
      act = i;
    }
  }
  /*
   *  Make sure only ONE pair is different
   */
  if (nActive != 1) {
    fprintf(stderr, "Error building ICylinder max and min are all same.\n");
    return false;
  }
  //
  //  If get here then points are different but only in 1 dimension and
  //  we stored the direction in act. Build the ICylinder and add it
  //  to the list.
  //
  newCyl = (Geom*) malloc(sizeof(Geom));
  assert(NULL != newCyl);
  ICylinderInit(newCyl, act, l->mArg);
  newCyl->mNext = l->mGList;
  l->mGList = newCyl;
  return true;
}

//
//  Make printable
//
void CD3ListPrintOn(CD3List* l, FILE* ofp)
{
  Geom* g = NULL;
  fprintf(ofp, "CD3List begin\n");
  for (g = l->mGList; NULL != g; g = g->mNext) {
    GeomPrintOn(g, ofp);
  }
  fprintf(ofp, "CD3List end\n");
}
//
//  Tell caller whether a point is inside the geometry (and
//  thus inactive) or not.
//
int CD3ListPointIn(CD3List* l, Point3D* p)
{
  Geom* g = NULL;
  int test;
  for (g = l->mGList; NULL != g; g = g->mNext) {
    test = false;
    switch (g->mId) {
      case kSD3ICyl:
        test = ICylinderPointIn(g, p);
        break;

      default:
        fprintf(stderr, "Unsupported geometry type %d.\n", g->mId);
        break;
    }
    if (test) {   // As soon as we see inside we are done
      return true;
    }
  }
  return false;
}
/*
 *  Write into a Smoothable.
 *
bool CD3ListAddGeomTo(CD3List* l, uint8_t* type, CD3Data* d)
{
  Geom* g = NULL;
  for (g = l->mGList; NULL != g; g = g->mNext) {
    switch (g->mId) {
      case kSD3ICyl:
        ICylinderAddTo(g, type, d);
        break;

      default:
        fprintf(stderr, "Unsupported geometry type %d.\n", g->mId);
        break;
    }
  }
  return true;
}
*/
