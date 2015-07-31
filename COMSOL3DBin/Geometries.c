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
      fprintf(ofp, "ICylinder from (%lfm %lf, %lf) to (%lf, %lf, %lf)\n",
              g->mMin.m[0], g->mMin.m[1], g->mMin.m[2],
              g->mMax.m[0], g->mMax.m[1], g->mMax.m[1]);
      fprintf(ofp, "Axis direction %d, radius %lf.\n",
              g->mAxis, sqrt(g->mR1Squared));
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
  g->mAxis = axis;
  g->mR1Squared = args[6] * args[6];
}
//
//  ICylinderPointIn returns true of the point falls inside
//  the cylinder itself.
//
int ICylinderPointIn(Geom* g, Point3D* p)
{
  double dx, dy, dz, rsq;
  double  tolsq = 0.01 * 0.01;
  switch (g->mAxis) {
    case 0:
      //
      //  Tube parallel to x axis.
      //  Start by testing x position.
      //
      if ((p->m[0] < g->mMin.m[0]) || (p->m[0] > g->mMax.m[0])) {
        return false;
      }
      //
      //  We're OK in the x direction. Try the radius.
      //
      dy = p->m[1] - g->mMin.m[1];
      dz = p->m[2] - g->mMin.m[2];
      rsq = dy*dy + dz*dz;
      if  (rsq < g->mR1Squared + tolsq) {
        return true;
      }
      break;

    case 1:
      //
      //  y axis.
      //  Start by testing y position.
      //
      if ((p->m[1] < g->mMin.m[1]) || (p->m[1] > g->mMax.m[1])) {
        return false;
      }
      //
      //  We're OK in the y direction. Try the radius.
      //
      dx = p->m[0] - g->mMin.m[0];
      dz = p->m[2] - g->mMin.m[2];
      rsq = dx*dx + dz*dz;
      if  (rsq < g->mR1Squared + tolsq) {
        return true;
      }
      break;

    case 2:
      //
      //  Tube parallel to z axis.
      //  Start by testing z position.
      //
      if ((p->m[2] <= g->mMin.m[2]) || (p->m[2] >= g->mMax.m[2])) {
        return false;
      }
      //
      //  We're OK in the z direction. Try the radius.
      //
      dx = p->m[0] - g->mMin.m[0];
      dy = p->m[1] - g->mMin.m[1];
      rsq = dx*dx + dy*dy;
      if  (rsq < g->mR1Squared + tolsq) {
        return true;
      }
      break;

    default:
      fprintf(stderr, "INVALID AXIS %d\n", g->mAxis);
      return false;
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
