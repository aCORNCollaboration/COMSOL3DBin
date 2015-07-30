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
//  BCollett 7/30/15 Alter the way that Geometries work. Instead of being
//  able to alter a pointType array all they know how to do is tell the
//  world whether a point is inside or outside the geometry. The array
//  takes care of passing in relevant points and updating the arrays.
//  To make this more efficient each Geometry adds the idea of a BoundingBox
//  that can be queried by the outside world.
//
//  Created by Brian Collett on 7/27/15.
//  Copyright (c) 2015 Brian Collett. All rights reserved.
//

#ifndef __Geometries__
#define __Geometries__

#include <stdio.h>
#include "COMSOLData3D.h"

//
//  Some simple structures to make it easier to work with 3D objects.
//  These parallel the C++ versions in my more usual

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
