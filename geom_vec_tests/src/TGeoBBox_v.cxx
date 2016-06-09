// @(#)root/geom:$Id: TGeoBBox_v.cxx 27731 2009-03-09 17:40:56Z brun $// Author: Andrei Gheata   24/10/01

// Contains() and DistromOutside/Out() implemented by Mihaela Gheata

/*************************************************************************
 * Copyright (C) 1995-2000, Rene Brun and Fons Rademakers.               *
 * All rights reserved.                                                  *
 *                                                                       *
 * For the licensing terms see $ROOTSYS/LICENSE.                         *
 * For the list of contributors see $ROOTSYS/README/CREDITS.             *
 *************************************************************************/

//--------------------------------------------------------------------------
// TGeoBBox_v - box class. All shape primitives inherit from this, their
//   constructor filling automatically the parameters of the box that bounds
//   the given shape. Defined by 6 parameters :
//      fDX, fDY, fDZ - half lengths on X, Y and Z axis
//      fOrigin[3]    - position of box origin
//
//--------------------------------------------------------------------------
//
//
//--- Building boxes
//  ==================
//  Normally a box has to be build only with 3 parameters : dx, dy, dz
// representing the half lengths on X, Y and Z axis. In this case, the origin
// of the box will match the one of its reference frame. The translation of the
// origin is used only by the constructors of all other shapes in order to
// define their own bounding boxes. Users should be aware that building a
// translated box that will represent a physical shape by itself will affect any
// further positioning of other shapes inside. Therefore in order to build a
// positioned box one should follow the recipe described in class TGeoNode.
//
// Creation of boxes
// 1.   TGeoBBox_v *box = new TGeoBBox_v("BOX", 20, 30, 40);
// Begin_Html
/*
<img src="gif/t_box.gif">
*/
// End_Html
//
// 2. A volume having a box shape can be built in one step:
//      TGeoVolume *vbox = gGeoManager->MakeBox("vbox", ptrMed, 20,30,40);
//
// Divisions of boxes.
//
//   Volumes having box shape can be divided with equal-length slices on
// X, Y or Z axis. The following options are supported:
// a) Dividing the full range of one axis in N slices
//      TGeoVolume *divx = vbox->Divide("SLICEX", 1, N);
//   - here 1 stands for the division axis (1-X, 2-Y, 3-Z)
// Begin_Html
/*
<img src="gif/t_boxdivX.gif">
*/
// End_Html
//
// b) Dividing in a limited range - general case.
//      TGeoVolume *divy = vbox->Divide("SLICEY",2,N,start,step);
//   - start = starting offset within (-fDY, fDY)
//   - step  = slicing step
//
// Begin_Html
/*
<img src="gif/t_boxdivstepZ.gif">
*/
// End_Html
//
// Both cases are supported by all shapes.
//   See also class TGeoShape for utility methods provided by any particular
// shape.
//_____________________________________________________________________________

#include <iostream>

#include "TGeoManager.h"
#include "TGeoMatrix.h"
#include "TGeoVolume.h"
#include "TVirtualGeoPainter.h"
#include "TGeoBBox_v.h"
#include "TVirtualPad.h"
#include "TBuffer3D.h"
#include "TBuffer3DTypes.h"
#include "TMath.h"
#include "TRandom.h"

using std::min;

//_____________________________________________________________________________
TGeoBBox_v::TGeoBBox_v() {
  // Default constructor
  SetShapeBit(TGeoShape::kGeoBox);
  fDX = fDY = fDZ = 0;
  fOrigin[0] = fOrigin[1] = fOrigin[2] = 0.0;
}

//_____________________________________________________________________________
TGeoBBox_v::TGeoBBox_v(double dx, double dy, double dz, double *origin) : TGeoBBox(dx, dy, dz, origin) {}

//_____________________________________________________________________________
TGeoBBox_v::TGeoBBox_v(const char *name, double dx, double dy, double dz, double *origin)
    : TGeoBBox(name, dx, dy, dz, origin) {}

//_____________________________________________________________________________
TGeoBBox_v::TGeoBBox_v(double *param) : TGeoBBox(param) {}

//_____________________________________________________________________________
TGeoBBox_v::~TGeoBBox_v() {
  // Destructor
}

//_____________________________________________________________________________
bool TGeoBBox_v::AreOverlapping(const TGeoBBox_v *box1, const TGeoMatrix *mat1, const TGeoBBox_v *box2,
                                  const TGeoMatrix *mat2) {
  // Check if 2 positioned boxes overlap.
  double master[3];
  double local[3];
  double ldir1[3], ldir2[3];
  const double *o1 = box1->GetOrigin();
  const double *o2 = box2->GetOrigin();
  // Convert center of first box to the local frame of second
  mat1->LocalToMaster(o1, master);
  mat2->MasterToLocal(master, local);
  if (TGeoBBox_v::Contains(local, box2->GetDX(), box2->GetDY(), box2->GetDZ(), o2))
    return true;
  double distsq = (local[0] - o2[0]) * (local[0] - o2[0]) + (local[1] - o2[1]) * (local[1] - o2[1]) +
                    (local[2] - o2[2]) * (local[2] - o2[2]);
  // Compute distance between box centers and compare with max value
  double rmaxsq = (box1->GetDX() + box2->GetDX()) * (box1->GetDX() + box2->GetDX()) +
                    (box1->GetDY() + box2->GetDY()) * (box1->GetDY() + box2->GetDY()) +
                    (box1->GetDZ() + box2->GetDZ()) * (box1->GetDZ() + box2->GetDZ());
  if (distsq > rmaxsq + TGeoShape::Tolerance())
    return false;
  // We are still not sure: shoot a ray from the center of "1" towards the
  // center of 2.
  double dir[3];
  mat1->LocalToMaster(o1, ldir1);
  mat2->LocalToMaster(o2, ldir2);
  distsq = sqrt(1. / distsq);
  dir[0] = (ldir2[0] - ldir1[0]) * distsq;
  dir[1] = (ldir2[1] - ldir1[1]) * distsq;
  dir[2] = (ldir2[2] - ldir1[2]) * distsq;
  mat1->MasterToLocalVect(dir, ldir1);
  mat2->MasterToLocalVect(dir, ldir2);
  // Distance to exit from o1
  double dist1 = TGeoBBox_v::DistFromInsideS(o1, ldir1, box1->GetDX(), box1->GetDY(), box1->GetDZ(), o1);
  // Distance to enter from o2
  double dist2 = TGeoBBox_v::DistFromOutsideS(local, ldir2, box2->GetDX(), box2->GetDY(), box2->GetDZ(), o2);
  if (dist1 > dist2)
    return true;
  return false;
}

//_____________________________________________________________________________
void TGeoBBox_v::ComputeNormal(const double *point, const double *dir, double *norm) const {
  // Computes normal to closest surface from POINT.
  memset(norm, 0, 3 * sizeof(double));
  double saf[3];
  int i;
  saf[0] = fabs(fabs(point[0] - fOrigin[0]) - fDX);
  saf[1] = fabs(fabs(point[1] - fOrigin[1]) - fDY);
  saf[2] = fabs(fabs(point[2] - fOrigin[2]) - fDZ);
  i = TMath::LocMin(3, saf);
  norm[i] = (dir[i] > 0) ? 1 : (-1);
}

// sw_____________________________________________________________________________
void TGeoBBox_v::ComputeNormal_l(const double *__restrict__ point, const double *__restrict__ dir,
                                 double *__restrict__ norm, int np) const {
  for (int i = 0; i < np; i++) //@EXPECTVEC
  {
    ComputeNormal(&point[3 * i], &dir[3 * i], &norm[3 * i]);
  }
}

// mb_____________________________________________________________________________
void TGeoBBox_v::ComputeNormal_v(const double *__restrict__ point, const double *__restrict__ dir,
                                 double *__restrict__ norm, int np) const {
  // Computes normal to closest surface from POINT.
  memset(norm, 0, 3 * np * sizeof(double));
  double saf[3];
  int min;
  for (int i = 0; i < np; i++) //@EXPECTVEC
  {
    saf[0] = fabs(fabs(point[3 * i] - fOrigin[0]) - fDX);
    saf[1] = fabs(fabs(point[3 * i + 1] - fOrigin[1]) - fDY);
    saf[2] = fabs(fabs(point[3 * i + 2] - fOrigin[2]) - fDZ);
    min = TMath::LocMin(3, saf);
    norm[3 * i + min] = (dir[3 * i + min] > 0) ? 1 : (-1);
  }
}

//_____________________________________________________________________________
bool TGeoBBox_v::CouldBeCrossed(const double *point, const double *dir) const {
  // Decides fast if the bounding box could be crossed by a vector.
  double mind = fDX;
  if (fDY < mind)
    mind = fDY;
  if (fDZ < mind)
    mind = fDZ;
  double dx = fOrigin[0] - point[0];
  double dy = fOrigin[1] - point[1];
  double dz = fOrigin[2] - point[2];
  double do2 = dx * dx + dy * dy + dz * dz;
  if (do2 <= (mind * mind))
    return true;
  double rmax2 = fDX * fDX + fDY * fDY + fDZ * fDZ;
  if (do2 <= rmax2)
    return true;
  // inside bounding sphere
  double doct = dx * dir[0] + dy * dir[1] + dz * dir[2];
  // leaving ray
  if (doct <= 0)
    return false;
  double dirnorm = dir[0] * dir[0] + dir[1] * dir[1] + dir[2] * dir[2];
  if ((doct * doct) >= (do2 - rmax2) * dirnorm)
    return true;
  return false;
}

void TGeoBBox_v::CouldBeCrossed_l(const double *point, const double *dir, bool *__restrict__ crossed,
                                  int np) const {
  // Decides fast if the bounding box could be crossed by a vector.
  for (unsigned int i = 0; i < np; ++i) {
    crossed[i] = CouldBeCrossed(&point[3 * i], &dir[3 * i]);
  }
}

void TGeoBBox_v::CouldBeCrossed_v(const StructOfCoord &point, const StructOfCoord &dir, bool *__restrict__ crossed,
                                  int np) const {
  // Decides fast if the bounding box could be crossed by a vector.
  double mind = min<double>(fDX, fDY);
  mind = min<double>(mind, fDZ);

  double rmax2 = fDX * fDX + fDY * fDY + fDZ * fDZ;

  for (int i = 0; i < np; i++) //@EXPECTVEC
  {
    double dx, dy, dz, do2;
    dx = fOrigin[0] - point.x[i];
    dy = fOrigin[1] - point.y[i];
    dz = fOrigin[2] - point.z[i];
    do2 = dx * dx + dy * dy + dz * dz;
    double doct = dx * dir.x[i] + dy * dir.y[i] + dz * dir.z[i];
    double dirnorm = dir.x[i] * dir.x[i] + dir.y[i] * dir.y[i] + dir.z[i] * dir.z[i];
    crossed[i] = ((do2 <= (mind * mind)) | (do2 <= rmax2) |
                  (doct > 0 & ((doct * doct) >= (do2 - rmax2) * dirnorm))); // to verify
    // !! use | instead of || .. as lazy operators are branching
  }
}

//_____________________________________________________________________________
void TGeoBBox_v::CouldBeCrossed_v(const double __restrict__ *point, const double __restrict__ *dir,
                                  bool *__restrict__ couldcrossed, int np) const {
  // Decides fast if the bounding box could be crossed by a vector.
  double mind = fDX;
  if (fDY < mind)
    mind = fDY;
  if (fDZ < mind)
    mind = fDZ;

  double dx, dy, dz, do2;
  double rmax2 = fDX * fDX + fDY * fDY + fDZ * fDZ;
  for (int i = 0; i < np; i++) //@EXPECTVEC
  {
    dx = fOrigin[0] - point[i * 3 + 0];
    dy = fOrigin[1] - point[i * 3 + 1];
    dz = fOrigin[2] - point[i * 3 + 2];

    do2 = dx * dx + dy * dy + dz * dz;
    // if (do2<=(mind*mind)) return true; // (1)

    // if (do2<=rmax2) return true; //(2)
    // inside bounding sphere
    double doct = dx * dir[3 * i + 0] + dy * dir[3 * i + 1] + dz * dir[3 * i + 2];
    // leaving ray
    // if (doct<=0) return false; //(3)
    double dirnorm =
        dir[3 * i + 0] * dir[3 * i + 0] + dir[3 * i + 1] * dir[3 * i + 1] + dir[3 * i + 2] * dir[3 * i + 2];
    // if ((doct*doct)>=(do2-rmax2)*dirnorm) return true; //(4)
    // return false; //(5)

    couldcrossed[i] = ((do2 <= (mind * mind)) | (do2 <= rmax2) |
                       (doct > 0 & ((doct * doct) >= (do2 - rmax2) * dirnorm))); // to verify
    // use | instead of || as lazy operators are branching
  }
}

//_____________________________________________________________________________
int TGeoBBox_v::DistancetoPrimitive(int px, int py) {
  // Compute closest distance from point px,py to each corner.
  const int numPoints = 8;
  return ShapeDistancetoPrimitive(numPoints, px, py);
}

//_____________________________________________________________________________
bool TGeoBBox_v::Contains(const double *point) const {
  // Test if point is inside this shape.
  if (fabs(point[2] - fOrigin[2]) > fDZ)
    return false;
  if (fabs(point[1] - fOrigin[1]) > fDY)
    return false;
  if (fabs(point[0] - fOrigin[0]) > fDX)
    return false;
  return true;
}

void TGeoBBox_v::Contains_l(const double *__restrict__ point, bool *__restrict__ isin, int np) const {
  //
  for (int i = 0; i < np; ++i) //@EXPECTVEC
  {
    isin[i] = this->TGeoBBox_v::Contains(&point[3 * i]);
  }
}

#define vector(elcount, type) __attribute__((vector_size((elcount) * sizeof(type)))) type
//_____________________________________________________________________________
void TGeoBBox_v::Contains_v(const double *__restrict__ pointi, bool *__restrict__ isin, int np) const {
// Test if point is inside this shape.
//  vector(32,double) *vx, *vy;
#ifndef __INTEL_COMPILER
  const double *point = (double *)__builtin_assume_aligned(pointi, 32);
#else
  const double *point = (const double *)pointi;
#endif

  //#pragma ivdep
  for (int i = 0; i < np; ++i) //@EXPECTVEC
  {
    double xx, yy, zz;
    xx = point[3 * i];
    yy = point[3 * i + 1];
    zz = point[3 * i + 2];
    isin[i] = (fabs(xx) < fDX) & (fabs(yy) < fDY) & (fabs(zz) < fDZ);
  }
}

#ifndef VEC_EXTENSIONS
void TGeoBBox_v::Contains_v(const StructOfCoord &__restrict__ pointi, bool *__restrict__ isin, int np) const {
// Test if point is inside this shape.
//  vector(32,double) *vx, *vy;
#ifndef __INTEL_COMPILER
  const double *x = (const double *)__builtin_assume_aligned(pointi.x, 16);
  const double *y = (const double *)__builtin_assume_aligned(pointi.y, 16);
  const double *z = (const double *)__builtin_assume_aligned(pointi.z, 16);
#else
  const double *x = (const double *)pointi.x;
  const double *y = (const double *)pointi.y;
  const double *z = (const double *)pointi.z;
#endif

  //#pragma ivdep
  for (int i = 0; i < np; ++i) //@EXPECTVEC
  {
    double xx, yy, zz;
    xx = x[i] - fOrigin[0];
    yy = y[i] - fOrigin[1];
    zz = z[i] - fOrigin[2];
    isin[i] = (fabs(xx) < fDX) & (fabs(yy) < fDY) & (fabs(zz) < fDZ);
  }
}
#endif

/*
//_____________________________________________________________________________
void TGeoBBox_v::Contains_v(const double * __restrict__ point, bool * __restrict__ isin, int np) const
{
 double px[np];
 double py[np];
 double pz[np];
 double fake[np];
 for(int i=0; i<np; ++i)
   {
     px[i]=point[3*i];
     py[i]=point[3*i+1];
     pz[i]=point[3*i+2];
     fake[i]=point[3*i+3];
   }

#pragma ivdep
 for(int i=0; i<np; ++i)  //@EXPECTVEC
   {
     double xx, yy, zz;
     xx=px[i]-fOrigin[0];
     yy=py[i]-fOrigin[1];
     zz=pz[i]-fOrigin[2];
     isin[i]=(fabs(xx)<fDX) & (fabs(yy)<fDY) & (fabs(zz)<fDZ);
   }
}
*/

//_____________________________________________________________________________
bool TGeoBBox_v::Contains(const double *point, double dx, double dy, double dz, const double *origin) {
  // Test if point is inside this shape.
  if (fabs(point[2] - origin[2]) > dz)
    return false;
  if (fabs(point[0] - origin[0]) > dx)
    return false;
  if (fabs(point[1] - origin[1]) > dy)
    return false;
  return true;
}

//__ static version of method _________________________________________________
double TGeoBBox_v::DistFromInsideS(const double *point, const double *dir, double dx, double dy, double dz,
                                     const double *origin, double /*stepmax*/) {
  // Computes distance from inside point to surface of the box.
  // Boundary safe algorithm.

  double s, smin, saf[6];
  double newpt[3];
  int i;

  newpt[0] = point[0] - origin[0];
  saf[0] = dx + newpt[0];
  saf[1] = dx - newpt[0];
  newpt[1] = point[1] - origin[1];
  saf[2] = dy + newpt[1];
  saf[3] = dy - newpt[1];
  newpt[2] = point[2] - origin[2];
  saf[4] = dz + newpt[2];
  saf[5] = dz - newpt[2];
  smin = TGeoShape::Big();

  double s1, s2;
  if (dir[0] != 0) {
    // sign will be zero for negative numbers
    s = (dir[0] > 0) ? (saf[1] / dir[0]) : (-saf[0] / dir[0]);
    if (s < smin)
      smin = s;
    if (s < 0)
      smin = 0.0;
  }

  if (dir[1] != 0) {
    s = (dir[1] > 0) ? (saf[3] / dir[1]) : (-saf[2] / dir[1]);
    if (s < smin)
      smin = s;
    if (s < 0)
      smin = 0.0;
  }

  if (dir[2] != 0) {
    s = (dir[2] > 0) ? (saf[5] / dir[2]) : (-saf[4] / dir[2]);
    if (s < smin)
      smin = s;
    if (s < 0)
      smin = 0.0;
  }

  return smin;
}

//_____________________________________________________________________________
double TGeoBBox_v::DistFromInside(double *point, double *dir, int iact, double step, double *safe) const {
  // Compute distance from inside point to surface of the box along some direction given by dir
  // Boundary safe algorithm.
  double s, smin, saf[6];
  double newpt[3];
  int i;
  for (i = 0; i < 3; i++)
    newpt[i] = point[i] - fOrigin[i];
  saf[0] = fDX + newpt[0];
  saf[1] = fDX - newpt[0];
  saf[2] = fDY + newpt[1];
  saf[3] = fDY - newpt[1];
  saf[4] = fDZ + newpt[2];
  saf[5] = fDZ - newpt[2];

  // tyring to vectorize the above already
  // idea: fDX, fDY, fDZ should be layed out in memory consecutively

  // double * fD = (double *) &fDX;
  // for (i=0; i<3; i++)
  //   {
  //     saf[2*i] = fD[i]+newpt[i];
  //     saf[2*i+1] = fD[i]-newpt[i];
  //   }
  if (iact < 3 && safe) {
    smin = saf[0];
    // compute safe distance
    for (i = 1; i < 6; i++)
      if (saf[i] < smin)
        smin = saf[i];
    *safe = smin;
    if (smin < 0)
      *safe = 0.0;
    if (iact == 0)
      return TGeoShape::Big();
    if (iact == 1 && step < *safe)
      return TGeoShape::Big();
  }

  // compute distance to surface
  smin = TGeoShape::Big();
  for (i = 0; i < 3; i++) {
    if (dir[i] != 0) {
      s = (dir[i] > 0) ? (saf[(i << 1) + 1] / dir[i]) : (-saf[i << 1] / dir[i]);
      if (s < 0)
        return 0.0;
      if (s < smin)
        smin = s;
    }
  }
  return smin;
}

//_____________________________________________________________________________
void TGeoBBox_v::DistFromInsideS_l(const double *__restrict__ point, const double *__restrict__ dir, double dx,
                                   double dy, double dz, const double *__restrict__ origin,
                                   const double *__restrict__ stepmax, double *__restrict__ distance, int npoints) {
#pragma simd
  // #pragma ivdep
  for (unsigned int k = 0; k < npoints; k++) //@EXPECTVEC
  {
    // maybe elemental function has to be declared declspec
    distance[k] = TGeoBBox_v::DistFromInsideS(&point[3 * k], &dir[3 * k], dx, dy, dz, origin, stepmax[k]);
  }
}

//_____________________________________________________________________________
void TGeoBBox_v::DistFromInsideS_l(const StructOfCoord &__restrict__ point, const StructOfCoord &__restrict__ dir,
                                   double dx, double dy, double dz, const double *__restrict__ origin,
                                   const double *__restrict__ stepmax, double *__restrict__ distance, int npoints) {
#pragma simd
  // #pragma ivdep
  for (unsigned int k = 0; k < npoints; k++) //@EXPECTVEC
  {
    // maybe elemental function has to be declared declspec
    double p[3] = {point.x[k], point.y[k], point.z[k]};
    double d[3] = {dir.x[k], dir.y[k], dir.z[k]};
    distance[k] = TGeoBBox_v::DistFromInsideS(p, d, dx, dy, dz, origin, stepmax[k]);
  }
}

//_____________________________________________________________________________
void TGeoBBox_v::DistFromInsideS_v(const double *__restrict__ point, const double *__restrict__ dir, double dx,
                                   double dy, double dz, const double *__restrict__ origin,
                                   const double *__restrict__ stepmax, double *__restrict__ distance, int npoints) {

  // this and the previous should be the same; here I have done manual inlining
  for (unsigned int k = 0; k < npoints; ++k) //@EXPECTVEC
  {
    double s, smin, saf[6];
    double newpt[3];
    int i;

    newpt[0] = point[3 * k + 0] - origin[0];
    saf[0] = dx + newpt[0];
    saf[1] = dx - newpt[0];
    newpt[1] = point[3 * k + 1] - origin[1];
    saf[2] = dy + newpt[1];
    saf[3] = dy - newpt[1];
    newpt[2] = point[3 * k + 2] - origin[2];
    saf[4] = dz + newpt[2];
    saf[5] = dz - newpt[2];

    smin = TGeoShape::Big();

    double s1, s2;
    int d = 3 * k;
    if (dir[d] != 0) {
      s = (dir[d] > 0) ? (saf[1] / dir[d]) : (-saf[0] / dir[d]);
    }
    if (s < smin)
      smin = s;
    if (s < 0)
      smin = 0.0;

    d = 3 * k + 1;
    if (dir[d] != 0) {
      s = (dir[d] > 0) ? (saf[3] / dir[d]) : (-saf[2] / dir[d]);
    }
    // if (s < smin) smin = (((int) s)>0) * s;
    if (s < smin)
      smin = s;
    if (s < 0)
      smin = 0.0;

    d = 3 * k + 2;
    if (dir[d] != 0) {
      s = (dir[d] > 0) ? (saf[5] / dir[d]) : (-saf[4] / dir[d]);
    }
    if (s < smin)
      smin = s; // smin = s;
    if (s < 0)
      smin = 0.0;

    // maybe elemental function has to be declared declspec
    distance[k] = smin;
  }
}

// SOA version____________________________________________________________________________
void TGeoBBox_v::DistFromInsideS_v(const StructOfCoord &__restrict__ point, const StructOfCoord &__restrict__ dir,
                                   double dx, double dy, double dz, const double *__restrict__ origin,
                                   const double *__restrict__ stepmax, double *__restrict__ distance, int npoints) {
#ifndef __INTEL_COMPILER
  const double *x = (const double *)__builtin_assume_aligned(point.x, 32);
  const double *y = (const double *)__builtin_assume_aligned(point.y, 32);
  const double *z = (const double *)__builtin_assume_aligned(point.z, 32);
  const double *dirx = (const double *)__builtin_assume_aligned(dir.x, 32);
  const double *diry = (const double *)__builtin_assume_aligned(dir.y, 32);
  const double *dirz = (const double *)__builtin_assume_aligned(dir.z, 32);
  double *dist = (double *)__builtin_assume_aligned(distance, 32);
#else
  const double *x = (const double *)point.x;
  const double *y = (const double *)point.y;
  const double *z = (const double *)point.z;
  const double *dirx = (const double *)dir.x;
  const double *diry = (const double *)dir.y;
  const double *dirz = (const double *)dir.z;
#endif

// this and the previous should be the same; here I have done manual inlining
#pragma vector aligned
  for (size_t k = 0; k < npoints; ++k) //@EXPECTVEC
  {
    double s, smin, saf[6];
    double newpt[3];

    newpt[0] = x[k] - origin[0];
    saf[0] = dx + newpt[0];
    saf[1] = dx - newpt[0];
    newpt[1] = y[k] - origin[1];
    saf[2] = dy + newpt[1];
    saf[3] = dy - newpt[1];
    newpt[2] = z[k] - origin[2];
    saf[4] = dz + newpt[2];
    saf[5] = dz - newpt[2];

    smin = TGeoShape::Big();
    double sx, sy, sz;
    double tiny = 1e-20;
    sx = (dirx[k] > 0) ? (saf[1] / (dirx[k] + tiny)) : (-saf[0] / (dirx[k] - tiny));
    sy = (diry[k] > 0) ? (saf[3] / (diry[k] + tiny)) : (-saf[2] / (diry[k] - tiny));
    sz = (dirz[k] > 0) ? (saf[5] / (dirz[k] + tiny)) : (-saf[4] / (dirz[k] - tiny));
    //      sx = (saf[1]/(std::fabs(dirx[k])+tiny));
    //      sy = (saf[3]/(std::fabs(diry[k])+tiny));
    //      sz = (saf[5]/(std::fabs(dirz[k])+tiny));

    smin = sx;
    smin = (sy < smin) ? sy : smin;
    smin = (sz < smin) ? sz : smin;
    distance[k] = (smin < 0) ? 0 : smin;
  }
}

//_____________________________________________________________________________
double TGeoBBox_v::DistFromOutside(const double *point, const double *dir, int iact, double step,
                                     double *safe) const

{
  // Compute distance from outside point to surface of the box.
  // Boundary safe algorithm.
  bool in = true;
  double saf[3];
  double par[3];
  double newpt[3];
  int i, j;
  for (i = 0; i < 3; i++)
    newpt[i] = point[i] - fOrigin[i];
  par[0] = fDX;
  par[1] = fDY;
  par[2] = fDZ;
  for (i = 0; i < 3; i++) {
    saf[i] = fabs(newpt[i]) - par[i];
    if (saf[i] >= step)
      return TGeoShape::Big();
    if (in && saf[i] > 0)
      in = false;
  }
  if (iact < 3 && safe) {
    // compute safe distance
    if (in) {
      *safe = 0.0;
    } else {
      *safe = saf[0];
      if (saf[1] > *safe)
        *safe = saf[1];
      if (saf[2] > *safe)
        *safe = saf[2];
    }
    if (iact == 0)
      return TGeoShape::Big();
    if (iact == 1 && step < *safe)
      return TGeoShape::Big();
  }
  // compute distance from point to box
  double coord, snxt = TGeoShape::Big();
  int ibreak = 0;
  // protection in case point is actually inside box
  if (in) {
    j = 0;
    double ss = saf[0];
    if (saf[1] > ss) {
      ss = saf[1];
      j = 1;
    }
    if (saf[2] > ss)
      j = 2;
    if (newpt[j] * dir[j] > 0)
      return TGeoShape::Big(); // in fact exiting
    return 0.0;
  }
  for (i = 0; i < 3; i++) {
    if (saf[i] < 0)
      continue;
    if (newpt[i] * dir[i] >= 0)
      continue;
    snxt = saf[i] / fabs(dir[i]);
    ibreak = 0;
    for (j = 0; j < 3; j++) {
      if (j == i)
        continue;
      coord = newpt[j] + snxt * dir[j];
      if (fabs(coord) > par[j]) {
        ibreak = 1;
        break;
      }
    }
    if (!ibreak)
      return snxt;
  }
  return TGeoShape::Big();
}

//_____________________________________________________________________________
void TGeoBBox_v::DistFromOutside_l(const double __restrict__ *point, const double __restrict__ *dir,
                                   double *__restrict__ distance, int np) const {
// trivial loop implementation calling "elemental function"
#pragma simd
  for (int i = 0; i < np; i++) //@EXPECTVEC
    distance[i] = TGeoBBox_v::DistFromOutside(&point[3 * i], &dir[3 * i]);
}

// static version_____________________________________________________________________________
double TGeoBBox_v::DistFromOutsideS(const double *point, const double *dir, double dx, double dy, double dz,
                                      const double *origin, double stepmax) {
  // Compute distance from outside point to surface of the box.
  // Boundary safe algorithm.
  bool in = true;
  double saf[3];
  double par[3];
  double newpt[3];
  int i, j;
  for (i = 0; i < 3; i++)
    newpt[i] = point[i] - origin[i];
  par[0] = dx;
  par[1] = dy;
  par[2] = dz;
  for (i = 0; i < 3; i++) {
    saf[i] = fabs(newpt[i]) - par[i];
    if (saf[i] >= stepmax)
      return TGeoShape::Big();
    if (in && saf[i] > 0)
      in = false;
  }
  // In case point is inside return ZERO
  if (in)
    return 0.0;
  double coord, snxt = TGeoShape::Big();
  int ibreak = 0;
  for (i = 0; i < 3; i++) {
    if (saf[i] < 0)
      continue;
    if (newpt[i] * dir[i] >= 0)
      continue;
    snxt = saf[i] / fabs(dir[i]);
    ibreak = 0;
    for (j = 0; j < 3; j++) {
      if (j == i)
        continue;
      coord = newpt[j] + snxt * dir[j];
      if (fabs(coord) > par[j]) {
        ibreak = 1;
        break;
      }
    }
    if (!ibreak)
      return snxt;
  }
  return TGeoShape::Big();
}

// SOA _l version of static method DistFromOutside
void TGeoBBox_v::DistFromOutsideS_l(const StructOfCoord &__restrict__ point, const StructOfCoord &__restrict__ dir,
                                    double dx, double dy, double dz, const double *__restrict__ origin,
                                    const double *__restrict__ stepmax, double *__restrict__ distance, int np) {
  for (unsigned int k = 0; k < np; ++k) // @EXPECTVEC
  {
    double p[3];
    double d[3];
    p[0] = point.x[k];
    p[1] = point.y[k];
    p[2] = point.z[k];
    d[0] = dir.x[k];
    d[1] = dir.y[k];
    d[2] = dir.z[k];
    distance[k] = TGeoBBox_v::DistFromOutsideS(p, d, dx, dy, dz, origin, stepmax[k]);
  }
}

// _l version of static method DistFromOutside
void TGeoBBox_v::DistFromOutsideS_l(const double *__restrict__ point, const double *__restrict__ dir, double dx,
                                    double dy, double dz, const double *__restrict__ origin,
                                    const double *__restrict__ stepmax, double *__restrict__ distance, int np) {
  for (unsigned int k = 0; k < np; ++k) // @EXPECTVEC
  {
    distance[k] = TGeoBBox_v::DistFromOutsideS(&point[3 * k], &dir[3 * k], dx, dy, dz, origin, stepmax[k]);
  }
}

#ifdef VEC_EXTENSIONS
#include "TGeoBBox_v_vecext_Vc.cxx"
#else // use normal implementation

// SOA version of static method DistFromOutside
void TGeoBBox_v::DistFromOutsideS_v(const StructOfCoord &__restrict__ point, const StructOfCoord &__restrict__ dir,
                                    double dx, double dy, double dz, const double *__restrict__ origin,
                                    const double *__restrict__ stepmax, double *__restrict__ distance, int np) {
// note: we take stepmax as a maxstep PER particle
// ALIGNMENTSTUFF HERE (TODO: should be generated automatically, MACRO? )
#ifndef __INTEL_COMPILER
  const double *p[3];
  p[0] = (const double *)__builtin_assume_aligned(point.x, 32); // 32 bit for AVX
  p[1] = (const double *)__builtin_assume_aligned(point.y, 32);
  p[2] = (const double *)__builtin_assume_aligned(point.z, 32);
  const double *d[3];
  d[0] = (const double *)__builtin_assume_aligned(dir.x, 32);
  d[1] = (const double *)__builtin_assume_aligned(dir.y, 32);
  d[2] = (const double *)__builtin_assume_aligned(dir.z, 32);
  double *dist = (double *)__builtin_assume_aligned(distance, 32);
  const double *stepm = (const double *)__builtin_assume_aligned(stepmax, 32);
#else
  const double *p[3];
  p[0] = (const double *)point.x;
  p[1] = (const double *)point.y;
  p[2] = (const double *)point.z;
  const double *d[3];
  d[0] = (const double *)dir.x;
  d[1] = (const double *)dir.y;
  d[2] = (const double *)dir.z;
  double *dist = (double *)distance;
#endif

  double par[3];
  par[0] = dx;
  par[1] = dy;
  par[2] = dz;
#pragma ivdep
  for (unsigned int k = 0; k < np; ++k) // @EXPECTVEC
  {
    bool in;
    double saf[3];
    double newpt[3];

    double factor = 1.;
    double infactor; // will be zero or one depending if we are inside/outside

    // for (unsigned int i=0; i<3; i++) {
    //	 newpt[i] = p[i][k] - origin[i];
    // saf[i] = fabs(newpt[i])-par[i];
    // factor* = (saf[i]>=stepmax[k]) ? TGeoShape::Big() : 1.;
    // if (in && saf[i]>0) in=false;
    //	 in = in & (saf[i]<0);
    //}
    // unrolled above block manually: ( it would be nice to have a template unrool loop and a lambda function ?? )

    newpt[0] = p[0][k] - origin[0];
    saf[0] = fabs(newpt[0]) - par[0];
    factor = (saf[0] >= stepmax[k]) ? TGeoShape::Big() : 1.; // this might be done at the end
    in = (saf[0] < 0);

    newpt[1] = p[1][k] - origin[1];
    saf[1] = fabs(newpt[1]) - par[1];
    factor *= (saf[1] >= stepmax[k]) ? TGeoShape::Big() : 1.; // this might be done at the end
    in = in & (saf[1] < 0);

    newpt[2] = p[2][k] - origin[2];
    saf[2] = fabs(newpt[2]) - par[2];
    factor *= (saf[2] >= stepmax[k]) ? TGeoShape::Big() : 1.; // this might be done at the end
    in = in & (saf[2] < 0);

    infactor = (double)!in;

    // NOW WE HAVE THE SAFETYS AND IF IN OUT

    double coord;
    double snxt[3] = {TGeoShape::Big(), TGeoShape::Big(), TGeoShape::Big()};

    /*
    int ibreak=0;
    for (i=0; i<3; i++) {
      if (saf[i]<0) continue;
      if (newpt[i]*dir[i] >= 0) continue;
      snxt = saf[i]/fabs(dir[i]);
      ibreak = 0;
      for (j=0; j<3; j++) {
        if (j==i) continue;
        coord=newpt[j]+snxt*dir[j];
        if (fabs(coord)>par[j]) {
          ibreak=1;
          break;
        }
      }
      if (!ibreak) return snxt;
      }
    */

    // THE FOLLOWING IS AN UNROLLED VERSION OF ABOVE CONSTRUCT WITHOUT EARLY RETURNS

    // i=0
    int hit0 = 0;
    if (saf[0] > 0 & newpt[0] * d[0][k] < 0) // if out and right direction
    {
      snxt[0] = saf[0] / fabs(d[0][k]);             // distance to y-z face
      double coord1 = newpt[1] + snxt[0] * d[1][k]; // calculate new y and z coordinate
      double coord2 = newpt[2] + snxt[0] * d[2][k];
      hit0 = (fabs(coord1) > par[1] | fabs(coord2) > par[2]) ? 0 : 1; // 0 means miss, 1 means hit
    }

    int hit1 = 0;
    if (saf[1] > 0 & newpt[1] * d[1][k] < 0) {
      snxt[1] = saf[1] / fabs(d[1][k]);
      double coord0 = newpt[0] + snxt[1] * d[0][k];
      double coord2 = newpt[2] + snxt[1] * d[2][k];
      hit1 = (fabs(coord0) > par[0] | fabs(coord2) > par[2]) ? 0 : 1;
    }

    int hit2 = 0;
    if (saf[2] > 0 & newpt[2] * d[2][k] < 0) {
      snxt[2] = saf[2] / fabs(d[2][k]);
      double coord0 = newpt[0] + snxt[2] * d[0][k];
      double coord1 = newpt[1] + snxt[2] * d[1][k];
      hit2 = (fabs(coord0) > par[0] | fabs(coord1) > par[1]) ? 0 : 1;
    }
    distance[k] = (hit0 | hit1 | hit2) ? factor * infactor * (hit0 * snxt[0] + hit1 * snxt[1] + hit2 * snxt[2])
                                       : infactor * TGeoShape::Big();
  }
}
#endif

//_____________________________________________________________________________
bool TGeoBBox_v::GetPointsOnFacet(int index, int npoints, double *array) const {
  // Fills array with n random points located on the surface of indexed facet.
  // The output array must be provided with a length of minimum 3*npoints. Returns
  // true if operation succeeded.
  // Possible index values:
  //    0 - all facets togeather
  //    1 to 6 - facet index from bottom to top Z
  if (index < 0 || index > 6)
    return false;
  double surf[6];
  double area = 0.;
  if (index == 0) {
    for (int isurf = 0; isurf < 6; isurf++) {
      surf[isurf] = TGeoBBox_v::GetFacetArea(isurf + 1);
      if (isurf > 0)
        surf[isurf] += surf[isurf - 1];
    }
    area = surf[5];
  }

  for (int i = 0; i < npoints; i++) {
    // Generate randomly a surface index if needed.
    double *point = &array[3 * i];
    int surfindex = index;
    if (surfindex == 0) {
      double val = area * gRandom->Rndm();
      surfindex = 2 + TMath::BinarySearch(6, surf, val);
      if (surfindex > 6)
        surfindex = 6;
    }
    switch (surfindex) {
    case 1:
      point[0] = -fDX + 2 * fDX * gRandom->Rndm();
      point[1] = -fDY + 2 * fDY * gRandom->Rndm();
      point[2] = -fDZ;
      break;
    case 2:
      point[0] = -fDX + 2 * fDX * gRandom->Rndm();
      point[1] = -fDY;
      point[2] = -fDZ + 2 * fDZ * gRandom->Rndm();
      break;
    case 3:
      point[0] = -fDX;
      point[1] = -fDY + 2 * fDY * gRandom->Rndm();
      point[2] = -fDZ + 2 * fDZ * gRandom->Rndm();
      break;
    case 4:
      point[0] = -fDX + 2 * fDX * gRandom->Rndm();
      point[1] = fDY;
      point[2] = -fDZ + 2 * fDZ * gRandom->Rndm();
      break;
    case 5:
      point[0] = fDX;
      point[1] = -fDY + 2 * fDY * gRandom->Rndm();
      point[2] = -fDZ + 2 * fDZ * gRandom->Rndm();
      break;
    case 6:
      point[0] = -fDX + 2 * fDX * gRandom->Rndm();
      point[1] = -fDY + 2 * fDY * gRandom->Rndm();
      point[2] = fDZ;
      break;
    }
  }
  return true;
}

//_____________________________________________________________________________
bool TGeoBBox_v::GetPointsOnSegments(int npoints, double *array) const {
  // Fills array with n random points located on the line segments of the shape mesh.
  // The output array must be provided with a length of minimum 3*npoints. Returns
  // true if operation is implemented.
  if (npoints < GetNmeshVertices()) {
    Error("GetPointsOnSegments", "You should require at least %d points", GetNmeshVertices());
    return false;
  }
  TBuffer3D &buff = (TBuffer3D &)GetBuffer3D(TBuffer3D::kRawSizes | TBuffer3D::kRaw, true);
  int npnts = buff.NbPnts();
  int nsegs = buff.NbSegs();
  // Copy buffered points  in the array
  memcpy(array, buff.fPnts, 3 * npnts * sizeof(double));
  int ipoints = npoints - npnts;
  int icrt = 3 * npnts;
  int nperseg = (int)(double(ipoints) / nsegs);
  double *p0, *p1;
  double x, y, z, dx, dy, dz;
  for (int i = 0; i < nsegs; i++) {
    p0 = &array[3 * buff.fSegs[3 * i + 1]];
    p1 = &array[3 * buff.fSegs[3 * i + 2]];
    if (i == (nsegs - 1))
      nperseg = ipoints;
    dx = (p1[0] - p0[0]) / (nperseg + 1);
    dy = (p1[1] - p0[1]) / (nperseg + 1);
    dz = (p1[2] - p0[2]) / (nperseg + 1);
    for (int j = 0; j < nperseg; j++) {
      x = p0[0] + (j + 1) * dx;
      y = p0[1] + (j + 1) * dy;
      z = p0[2] + (j + 1) * dz;
      array[icrt++] = x;
      array[icrt++] = y;
      array[icrt++] = z;
      ipoints--;
    }
  }
  return true;
}

//_____________________________________________________________________________
double TGeoBBox_v::Safety(const double *point, bool in) const {

  // Computes the closest distance from given point to this shape.

  /*
     double safe, safy, safz;
     if (in) {
        safe = fDX - fabs(point[0]-fOrigin[0]);
        safy = fDY - fabs(point[1]-fOrigin[1]);
        safz = fDZ - fabs(point[2]-fOrigin[2]);
        if (safy < safe) safe = safy;
        if (safz < safe) safe = safz;
     } else {
        safe = -fDX + fabs(point[0]-fOrigin[0]);
        safy = -fDY + fabs(point[1]-fOrigin[1]);
        safz = -fDZ + fabs(point[2]-fOrigin[2]);
        if (safy > safe) safe = safy;
        if (safz > safe) safe = safz;
     }
  */

  double safe, safy, safz;
  char insign = (in) ? 1 : -1;
  safe = insign * (fDX - fabs(point[0] - fOrigin[0]));
  safy = insign * (fDY - fabs(point[1] - fOrigin[1]));
  safz = insign * (fDZ - fabs(point[2] - fOrigin[2]));
  if (safy < safe)
    safe = safy;
  if (safz < safe)
    safe = safz;

  return safe;
}

//_____________________________________________________________________________
void TGeoBBox_v::Safety_l(const double *__restrict__ point, double *__restrict__ safety, const int n,
                          bool in) const {
  for (unsigned int i = 0; i < n; i++) {
    safety[i] = this->Safety(&point[3 * i], in);
  }
}

// prototype for vectorized min function
inline void getmin_v(const double *x, const double *y, double *z, int n) {
  //  return the component-wise min of two vectors
  for (unsigned int i = 0; i < n; i++) // this vectorizes
  {
    z[i] = (x[i] > y[i]) ? y[i] : x[i];
  }
}

inline double myabs(double x) {
  return (x >= 0) ? x : -x; // this is basically the same as ROOT
}

//_____________________________________________________________________________
void TGeoBBox_v::Safety_v(const double *__restrict__ point, bool in, double *__restrict__ safety,
                          const int n) const {
  double *fD = (double *)&fDX; // treat box sizes as vector instead of writing explicitly x,y,z ( this works because
                               // they are next to each other in memory

  // treat in as a simple sign
  double insign = (in - 0.5) * 2.;
#pragma ivdep
  for (unsigned int i = 0; i < n; i++) // @EXPECTVEC
  {
    double t = point[3 * i + 0] - fOrigin[0];
    safety[i] = insign * (fD[0] - myabs(t));
  }

#pragma ivdep
  for (unsigned int i = 0; i < n; i++) // @EXPECTVEC
  {
    double t = point[3 * i + 1] - fOrigin[1];
    double t2 = insign * (fD[1] - myabs(t));
    safety[i] = min<double>(safety[i], t2);
  }

#pragma ivdep
  for (unsigned int i = 0; i < n; i++) // @EXPECTVEC
  {
    double t = point[3 * i + 2] - fOrigin[2];
    double t2 = insign * (fD[2] - myabs(t));
    safety[i] = min<double>(safety[i], t2);
  }
}

// SOA version_____________________________________________________________________________
void TGeoBBox_v::Safety_v(const StructOfCoord &__restrict__ point, bool in, double *__restrict__ safety,
                          const int n) const {

#ifndef __INTEL_COMPILER
  const double *x = (const double *)__builtin_assume_aligned(point.x, 16);
  const double *y = (const double *)__builtin_assume_aligned(point.y, 16);
  const double *z = (const double *)__builtin_assume_aligned(point.z, 16);
#else
  const double *x = (const double *)point.x;
  const double *y = (const double *)point.y;
  const double *z = (const double *)point.z;
#endif

  double safe, safy, safz;
  char insign = (in) ? 1 : -1;
  for (unsigned int i = 0; i < n; i++) // @EXPECTVEC
  {
    safe = insign * (fDX - fabs(x[i] - fOrigin[0]));
    safy = insign * (fDY - fabs(y[i] - fOrigin[1]));
    safz = insign * (fDZ - fabs(z[i] - fOrigin[2]));
    if (safy < safe)
      safe = safy;
    if (safz < safe)
      safe = safz;
    safety[i] = safe;
  }
}

//_____________________________________________________________________________
void TGeoBBox_v::Safety_v(const double *__restrict__ point, double *__restrict__ safety,
                          const bool *__restrict__ in, const int n) const {
  double *fD = (double *)&fDX; // treat box sizes as vector instead of writing explicitly x,y,z ( this works because
                               // they are next to each other in memory )

#pragma ivdep
  for (unsigned int i = 0; i < n; i++) // @EXPECTVEC
  {
    double t = point[3 * i + 0] - fOrigin[0];
    safety[i] = 2. * (in[i] - 0.5) * (fD[0] - myabs(t));
  }

#pragma ivdep
  for (unsigned int i = 0; i < n; i++) // @EXPECTVEC
  {
    double t = point[3 * i + 1] - fOrigin[1];
    double t2 = 2. * (in[i] - 0.5) * (fD[1] - myabs(t));
    safety[i] = min<double>(safety[i], t2);
  }

#pragma ivdep
  for (unsigned int i = 0; i < n; i++) // @EXPECTVEC
  {
    double t = point[3 * i + 2] - fOrigin[2];
    double t2 = 2. * (in[i] - 0.5) * (fD[2] - myabs(t));
    safety[i] = min<double>(safety[i], t2);
  }
}
