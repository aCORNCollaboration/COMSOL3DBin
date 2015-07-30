//
//  Geometries.h
//  Smooth3D
//
//  Pure C version for COMSOL3DBin.
//  This is somewhat different as it has to use switches and the internal
//  ID instead of inheritance to distinguish different geometries.
//
//  These objects represent geometries that can be added to the pointType
//  array of a Smoothable3D. They are designed to be stored in lists
//  (in a CD3List) and know how to alter the values in the pointType
//  array of a storable to mark some regions as boundaries.
//
//  Created by Brian Collett on 7/27/15.
//  Copyright (c) 2015 Brian Collett. All rights reserved.
//

#ifndef __Geometries__
#define __Geometries__

#include <stdio.h>
#include "COMSOLData3D.h"

//
//	Enum for the different kinds of OpenGL command that can be
//	found in a .gla file.
//
typedef enum {
  kSD3Empty = 0,
  kSD3ICyl,
  kSD3End,     // end of geometry. Allows sharing geom file with other.
  kSD3Error
} SD3Command;

//
//  Structures holding representations of individual geometries.
//  Since the base class is abstract no-one can ever make one.
typedef struct GeomTag {
  int mId;
  struct GeomTag* mNext;
  int mAxis;        // 0-2
  double mParams[20];
} Geom;
//
int GeomInit(Geom* g, int id);
int GeomFinish(Geom* g);
//
//  Make printable
//
void GeomPrintOn(Geom* g, FILE* ofp);
//
//  Constructor for an ICylinder.
//
void ICylinderInit(Geom* g, int axis, double* args);
//
//  Add to a type array.
//
int ICylinderAddTo(Geom* g, uint8_t* type, CD3Data* d);

//
//  Access names for the parameters for a cylinder
//
enum {
  kICylXMin = 0,
  kICylYMin,
  kICylZMin,
  kICylXMax,
  kICylYMax,
  kICylZMax,
  kICylRadius,
  kICylPotential,
  kICylError
};

#endif /* defined(__Geometries__) */
