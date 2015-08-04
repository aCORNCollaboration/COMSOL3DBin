//
//  Geometries.cpp
//  Smooth3D
//
//  Pure C version for working with COMSOL3dBin.
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
//#include "Smooth3D.h"
#include <math.h>
#include <stdlib.h>
#include "assert.h"
#include "Geometries.h"
//#include "Smoothable3D.h"
//
//  Geometries.
//
int GeomInit(Geom* g, int id)
{
  g->mId = id;
  g->mNext = NULL;
  return 0;
}

int GeomFinish(Geom* g)
{
  if (g->mNext) {
    free(g->mNext);
  }
  return 0;
}
//
//  Make printable
//
void GeomPrintOn(Geom* g, FILE* ofp)
{
  assert(NULL != ofp);
  assert(NULL != g);
  switch(g->mId) {
    case kSD3ICyl:
      fprintf(ofp, "ICylinder from (%f, %f, %f) to (%f, %f, %f)\n",
              g->mMin.m[0], g->mMin.m[1], g->mMin.m[2],
              g->mMax.m[0], g->mMax.m[1], g->mMax.m[1]);
      fprintf(ofp, "Axis indices (%d, %d, %d) radius %f.\n",
              g->mIdx0, g->mIdx1, g->mIdx2, sqrt(g->mR1Squared));
      break;

    case kSD3Torus:
      fprintf(ofp, "Torus from (%f, %f, %f) to (%f, %f, %f)\n",
              g->mMin.m[0], g->mMin.m[1], g->mMin.m[2],
              g->mMax.m[0], g->mMax.m[1], g->mMax.m[1]);
      fprintf(ofp, "Axis indices (%d, %d, %d) radius1 %f radius2 %f.\n",
              g->mIdx0, g->mIdx1, g->mIdx2,
              sqrt(g->mR1Squared), sqrt(g->mR2Squared));
      break;
      
    default:
      fprintf(ofp, "Raw Geometry ID = %d\n", g->mId);
      break;
  }
}
/*
ICylinder::ICylinder(int idx, double *args) : Geom(kSD3ICyl)
{
  mAxis = idx;
  for (int i = 0; i < 3; i++) {
    mMin[i] = args[i];
    mMax[i] = args[i+3];
  }
  mRadius = args[6];
  mPotential = args[7];
}mMax[
//
//  Make printable
void ICylinder::PrintOn(FILE* ofp)
{
}
 */
//
//  Constructor for an ICylinder.
//  Copy the parameters into their proper places and compute
//  the bounding box.
//  Params xMin, yMin, zMin, xMax, yMax, zMax, radius
//
void ICylinderInit(Geom* g, int axis, double* args)
{
  int i;
  GeomInit(g, kSD3ICyl);
  for (i = 0; i < 3; i++) {
    g->mMin.m[i] = args[i];
    g->mMax.m[i] = args[i+3];
  }
  switch (axis) {
    case 0:
      //
      //  Axis is x, lowest is y, mid is z
      //
      g->mIdx0 = 1;
      g->mIdx1 = 2;
      g->mIdx2 = 0;
      break;
      
    case 1:
      //
      //  Axis is y, lowest is z, mid is x
      //
      g->mIdx0 = 2;
      g->mIdx1 = 0;
      g->mIdx2 = 1;
      break;
      
    case 2:
      //
      //  Axis is z, lowest is x, mid is y
      //
      g->mIdx0 = 0;
      g->mIdx1 = 1;
      g->mIdx2 = 2;
      break;
      
    default:
      fprintf(stderr, "INVALID AXIS %d\n", axis);
      g->mIdx0 = 0;
      g->mIdx1 = 1;
      g->mIdx2 = 2;
  };
  g->mR1Squared = args[6] * args[6];
}
//
//  ICylinderPointIn returns true of the point falls inside
//  the cylinder itself.
//
int ICylinderPointIn(Geom* g, Point3D* p, double tol)
{
  double dx0, dx1, rsq;
  double  tolsq = tol * tol;
  //
  //  Start by testing axial position.
  //
  if ((p->m[g->mIdx2] < g->mMin.m[g->mIdx2]) ||
      (p->m[g->mIdx2] > g->mMax.m[g->mIdx2])) {
    return false;
  }
  //
  //  We're OK in the x direction. Try the radius.
  //
  dx0 = p->m[g->mIdx0] - g->mMin.m[g->mIdx0];
  dx1 = p->m[g->mIdx1] - g->mMin.m[g->mIdx1];
  rsq = dx0*dx0 + dx1*dx1;
  if  (rsq < g->mR1Squared + tolsq) {
    return true;
  }
  return false;
}

//
//  Constructor for a Torus.
//  Copy the parameters into their proper places and compute
//  the bounding box.
//  Params xMin, yMin, zMin, xMax, yMax, zMax, radius1, radius2
//
void TorusInit(Geom* g, int axis, double* args)
{
  int i;
  GeomInit(g, kSD3Torus);
  for (i = 0; i < 3; i++) {
    g->mMin.m[i] = args[i];
    g->mMax.m[i] = args[i+3];
  }
  switch (axis) {
    case 0:
      //
      //  Axis is x, lowest is y, mid is z
      //
      g->mIdx0 = 1;
      g->mIdx1 = 2;
      g->mIdx2 = 0;
      break;
      
    case 1:
      //
      //  Axis is y, lowest is z, mid is x
      //
      g->mIdx0 = 2;
      g->mIdx1 = 0;
      g->mIdx2 = 1;
      break;
      
    case 2:
      //
      //  Axis is z, lowest is x, mid is y
      //
      g->mIdx0 = 0;
      g->mIdx1 = 1;
      g->mIdx2 = 2;
      break;
      
    default:
      fprintf(stderr, "INVALID AXIS %d\n", axis);
      g->mIdx0 = 0;
      g->mIdx1 = 1;
      g->mIdx2 = 2;
  };
  g->mR1Squared = args[6] * args[6];
  g->mR2Squared = args[7] * args[7];
}
//
//  TorusPointIn returns true of the point falls inside
//  the outer radius and outside the inner as well as in
//  range axially.
//
int TorusPointIn(Geom* g, Point3D* p, double tol)
{
  double dx0, dx1, rsq;
  double  tolsq = tol * tol;
  //  Start by testing axial position.
  //
  if ((p->m[g->mIdx2] < g->mMin.m[g->mIdx2]) ||
      (p->m[g->mIdx2] > g->mMax.m[g->mIdx2])) {
    return false;
  }
  //
  //  We're OK in the x direction. Try the radius.
  //
  dx0 = p->m[g->mIdx0] - g->mMin.m[g->mIdx0];
  dx1 = p->m[g->mIdx1] - g->mMin.m[g->mIdx1];
  rsq = dx0*dx0 + dx1*dx1;
  if  ((rsq < g->mR2Squared + tolsq) && (rsq > g->mR1Squared - tolsq)) {
    return true;
  }
  return false;
}

/*
 *  Add to a type array.
 *
int ICylinderAddTo(Geom* g, uint8_t* type, CD3Data* d)
{
  uint32_t iMin[3], iMax[3];
  double r = g->mParams[kICylRadius];
  double rSq = r*r, delsq = d->mDelta[0]*d->mDelta[0];
  uint32_t idx0, idx1, idx2;
  uint32_t ix0, iy0, iz0;
  //
  //  Clip to this Smoothable and see that there is something left!
  //
  double newMin[3], newMax[3];
  CD3ClipPt(d, g->mParams+kICylXMin, newMin);
  CD3ClipPt(d, g->mParams+kICylXMax, newMax);
  if (newMin[g->mAxis] >= newMax[g->mAxis]) {
    fprintf(stderr, "After clipping ICylinder has no volume.\n");
    return false;
  }
  //
  //  Now for the hard part. We have to figure out which points in the
  //  Smoothable are inside the cylinder and mark them.This is different
  //  for each of the three different axes.
  //
  switch (g->mAxis) {
    case 0:
      //
      //  Tube parallel to x axis.
      //
      //
      fprintf(stderr, "ICylinder: x axis not supported yet.\n");
      break;

    case 1:
      //
      //  y axis.
      //
      fprintf(stderr, "ICylinder: y axis not supported yet.\n");
      break;

    case 2:
      //
      //  Tube parallel to z axis.
      //
      assert(newMin[0] == newMax[0]);
      assert(newMin[1] == newMax[1]);
      assert(newMin[2] < newMax[2]);
      //
      //  Have to find limits in x and y directions.
      //  Convert a couple of points to find x and y ranges.
      //  Start by making lowest x and y and highest.
      //
      newMin[0] -= r;
      newMin[1] -= r;
      newMax[0] += r;
      newMax[1] += r;
      fprintf(stdout, "ICyl bound box (%lf,%lf,%lf)->(%lf,%lf,%lf)\n",
              newMin[0], newMin[1], newMin[2], newMax[0], newMax[1], newMax[2]);
      //
      //  Now map the max and min into indices.
      //
      if (!CD3Map(d, newMin, iMin)) {
        fprintf(stderr, "ICylinder:Min:Index conversion failure.\n");
        return false;
      }
      if (!CD3Map(d, newMax, iMax)) {
        fprintf(stderr, "ICylinder:Min:Index conversion failure.\n");
        return false;
      }
      ix0 = (iMax[0] + iMin[0])/2;
      iy0 = (iMax[1] + iMin[1])/2;
      //
      //  Now run through all the points in the parallelipiped and mark
      //  those that fall inside the tube.
      //  Work over the x axis first, then y, then z.
      //
      fprintf(stdout, "ICyl bound box (%d,%d,%d)->(%d,%d,%d)\n",
              iMin[0], iMin[1], iMin[2],iMax[0], iMax[1], iMax[2]);
      for (idx0 = iMin[0]; idx0 <= iMax[0]; idx0++) {
        for (idx1 = iMin[1]; idx1 <= iMax[1]; idx1++) {
          //
          //  Build coords of point and check if it is in cylinder.
          //
          double distsq = (idx0 - ix0)*(idx0 - ix0)*delsq +
                          (idx1 - iy0)*(idx1 - iy0)*delsq;
          if (distsq > rSq) {
            continue;
          }
          //
          //  Then do set every point at this x-y.
          //
          for (idx2 = iMin[2]; idx2 <= iMax[2]; idx2++) {
            type[CD3IndexAt(d, idx0, idx1, idx2)] = 0;
          }
        }
      }
      break;

    default:
      fprintf(stderr, "INVALID AXIS %d\n", g->mAxis);
      return false;
  }
  return true;
}
*/
