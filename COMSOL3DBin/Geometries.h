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
//  that can be queried by the outside world. Restructured the Geom
//  parameters into named groups at the same time.
//  BCollett 8/1/15 Change the way axis info is stored in Cyls and Torii.
//  Instead of 1 mAxis variable have three mIdx0, mIdx1, mIdx2 that are
//  the indices into the coordinate arrays organized so that idx2 points
//  along the axis of the cylinder and idx0 and 1 form a right-handed
//  coordinate system.
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
typedef struct Point3DTag {
  double m[3];
} Point3D;

typedef struct Vector3DTag {
  double m[3];
} Vector3D;
//
//  Access macros
//
#define mX(p) (p)->m[0]
#define mY(p) (p)->m[1]
#define mZ(p) (p)->m[2]

typedef struct Box3DTag {
  Point3D mMin;
  Point3D mMax;
} Box3D;
//
//	Enum for the different kinds of OpenGL command that can be
//	found in a .gla file.
//
typedef enum {
  kSD3Empty = 0,
  kSD3ICyl,     // Volume inside a cylinder aligned with axes.
  kSD3Torus,   // Volume inside a square-sectioned torus
  kSD3End,      // end of geometry. Allows sharing geom file with other.
  kSD3Error
} SD3Command;

//
//  Structures holding representations of individual geometries.
//
typedef struct GeomTag {
  int mId;
  struct GeomTag* mNext;
  Box3D mBounds;
  Point3D mMin;
  Point3D mMax;
  int mIdx0, mIdx1, mIdx2;  // Idx2 is axis of symmetry if present
  double mR1Squared;        // Store radii as squares
  double mR2Squared;        // Torii have two.
} Geom;
//
int GeomInit(Geom* g, int id);
int GeomFinish(Geom* g);
//
//  Make printable
//
void GeomPrintOn(Geom* g, FILE* ofp);
//
//  Individual actual geometries support a few more functions.
//
//
//  Constructor for an ICylinder.
//
void ICylinderInit(Geom* g, int axis, double* args);
//
//  ICylinderPointIn returns true of the point falls inside
//  the cylinder itself.
//
int ICylinderPointIn(Geom* g, Point3D* p, double tol);
//
//  Constructor for an Torus.
//
void TorusInit(Geom* g, int axis, double* args);
//
//  ICylinderPointIn returns true of the point falls inside
//  the cylinder itself.
//
int TorusPointIn(Geom* g, Point3D* p, double tol);
//
//  Old form functionality removed 7/30/15
//  Add to a type array.
//
//int ICylinderAddTo(Geom* g, uint8_t* type, CD3Data* d);

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
