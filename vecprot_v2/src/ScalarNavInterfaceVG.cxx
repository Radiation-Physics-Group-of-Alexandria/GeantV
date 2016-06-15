#include "ScalarNavInterfaceVG.h"

#include "backend/Backend.h"
#ifdef GEANT_NVCC
#include "navigation/SimpleNavigator.h"
#else
#include "navigation/ABBoxNavigator.h"
#endif
#include "volumes/PlacedVolume.h"
#include "base/Vector3D.h"
#include "base/Transformation3D.h"
#include "base/Global.h"
#include "management/GeoManager.h"

#ifndef GEANT_NVCC
#ifdef CROSSCHECK
#include "TGeoNavigator.h"
#include "TGeoNode.h"
#endif
#endif

#ifdef BUG_HUNT
#include "GeantPropagator.h"
#endif

namespace Geant {
inline namespace GEANT_IMPL_NAMESPACE {

using namespace VECGEOM_NAMESPACE;

//______________________________________________________________________________
GEANT_CUDA_BOTH_CODE
void ScalarNavInterfaceVG::NavFindNextBoundaryAndStep(int ntracks, const double *pstep, 
         const double *x, const double *y, const double *z,
         const double *dirx, const double *diry, const double *dirz, 
         const VolumePath_t **instate, VolumePath_t **outstate, 
         double *step, double *safe, bool *isonbdr) {
// Find the next boundary and state after propagating to the boundary. 
// Input:  ntracks - number of tracks
//         pstep - proposed step
//         x, y, z, dirx, diry, dirz - particle position and direction
//         instate - input particle navigation state
//         safe - estimated safety value for the input point
//         isonbdr - starting point is on a boundary
// Output: outstate - navigation state after propagation. If boundary further
//           than proposed step, outstate has to match instate
//         step - propagation step for which the state is sampled
//         safety - calculated safety value for the input point
//         isonbdr - propagated point is on a boundary

  typedef Vector3D<Precision> Vector3D_t;
#ifdef GEANT_NVCC
  SimpleNavigator nav;
#else
  ABBoxNavigator nav;
#endif // GEANT_NVCC

  for (int itr = 0; itr < ntracks; ++itr) {
    // Check if current safety allows for the proposed step
    if (safe[itr] > pstep[itr]) {
      step[itr] = pstep[itr];
      isonbdr[itr] = false;
      *outstate[itr] = *instate[itr];
      continue;
    }
    
    nav.FindNextBoundaryAndStep(Vector3D_t(x[itr], y[itr], z[itr]),
                                Vector3D_t(dirx[itr], diry[itr], dirz[itr]), 
                                *instate[itr], *outstate[itr], 
                                Math::Min<double>(1.E20, pstep[itr]), step[itr]);
    step[itr] = Math::Max<double>(2. * gTolerance, step[itr] + 2. * gTolerance);
    safe[itr] = (isonbdr[itr]) ? 0 : nav.GetSafety(Vector3D_t(x[itr], y[itr], z[itr]), *instate[itr]);
    safe[itr] = Math::Max<double>(safe[itr], 0);
    // onboundary with respect to new point
    isonbdr[itr] = outstate[itr]->IsOnBoundary();
    
    //#### To add small step detection and correction - see ScalarNavInterfaceTGeo ####//

#ifndef GEANT_NVCC
#ifdef CROSSCHECK
    //************
    // CROSS CHECK USING TGEO
    //************
    TGeoNavigator *rootnav = gGeoManager->GetCurrentNavigator();
    rootnav->ResetState();
    rootnav->SetCurrentPoint(x[itr], y[itr], z[itr]);
    rootnav->SetCurrentDirection(dirx[itr], diry[itr], dirz[itr]);
    TGeoBranchArray *tmp = instate[itr]->ToTGeoBranchArray();
    tmp->UpdateNavigator(rootnav);
    delete tmp;
    rootnav->FindNextBoundaryAndStep(Math::Min<double>(1.E20, pstep[itr]), !isonbdr[itr]);
    double stepcmp = Math::Max<double>(2 * gTolerance, rootnav->GetStep());
    double safecmp = rootnav->GetSafeDistance();
    if (Math::Abs(step[itr] - stepcmp) > 1E-6) {
      Geant::Print("","## PSTEP %lf VECGEOMSTEP %lf ROOTSTEP %lf", pstep[itr], step[itr], stepcmp);
      Geant::Print("","## PSTEP %lf ONBOUND %d VECGEOMSAFETY %lf ROOTSAFETY %lf BRUTEFORCEROOT %lf", pstep[itr], isonbdr[itr],
             safe[itr], rootnav->Safety());

      // check nextpath
      tmp = outstate[itr]->ToTGeoBranchArray();
      tmp->InitFromNavigator(rootnav);
      Geant::Print("","## VECGEOMNEXTNODE %p ROOTNEXTNODE %p", outstate[itr]->GetCurrentNode(), tmp->GetCurrentNode());
      Geant::Print("","## VECGEOMBOUNDARY %d ROOTBOUNDARY %d", outstate[itr]->IsOnBoundary(), rootnav->IsOnBoundary());

      Geant::Print("","INCONSISTENT STEP");
      nav.InspectEnvironmentForPointAndDirection(Vector3D_t(x[itr], y[itr], z[itr]) /*global pos*/,
                                                 Vector3D_t(dirx[itr], diry[itr], dirz[itr]), *instate[itr]);
    }
#endif // CROSSCHECK
#endif 

#ifdef VERBOSE
    Geant::Print("","navfindbound on %p track %d with pstep %lf yields step %lf and safety %lf\n", this, itr, pstep[itr], step[itr],
           safe[itr]);
#endif // VERBOSE
  }
}

//______________________________________________________________________________
GEANT_CUDA_BOTH_CODE
void ScalarNavInterfaceVG::NavFindNextBoundaryAndStep(const double &pstep, 
         const double &x, const double &y, const double &z,
         const double &dirx, const double &diry, const double &dirz, 
         const VolumePath_t *instate, VolumePath_t *outstate, 
         double &step, double &safe, bool &isonbdr) {
// Find the next boundary and state after propagating to the boundary. 
// Input:  pstep - proposed step
//         x, y, z, dirx, diry, dirz - particle position and direction
//         instate - input particle navigation state
//         safe - estimated safety value for the input point
//         isonbdr - starting point is on a boundary
// Output: outstate - navigation state after propagation. If boundary further
//           than proposed step, outstate has to match instate
//         step - propagation step for which the state is sampled
//         safety - calculated safety value for the input point
//         isonbdr - propagated point is on a boundary

  typedef Vector3D<Precision> Vector3D_t;
#ifdef GEANT_NVCC
  SimpleNavigator nav;
#else
  ABBoxNavigator nav;
#endif // GEANT_NVCC

  // Check if current safety allows for the proposed step
  if (safe > pstep) {
    step = pstep;
    isonbdr = false;
    *outstate = *instate;
    return;
  }

  nav.FindNextBoundaryAndStep(Vector3D_t(x, y, z),
                              Vector3D_t(dirx, diry, dirz), 
                              *instate, *outstate, 
                              Math::Min<double>(1.E20, pstep), step);
  step = Math::Max<double>(2. * gTolerance, step + 2. * gTolerance);
  safe = (isonbdr) ? 0 : nav.GetSafety(Vector3D_t(x, y, z), *instate);
  safe = Math::Max<double>(safe, 0);
  // onboundary with respect to new point
  isonbdr = outstate->IsOnBoundary();

  //#### To add small step detection and correction - see ScalarNavInterfaceTGeo ####//

#ifdef CROSSCHECK
  //************
  // CROSS CHECK USING TGEO
  //************
  TGeoNavigator *rootnav = gGeoManager->GetCurrentNavigator();
  rootnav->ResetState();
  rootnav->SetCurrentPoint(x, y, z);
  rootnav->SetCurrentDirection(dirx, diry, dirz);
  TGeoBranchArray *tmp = instate->ToTGeoBranchArray();
  tmp->UpdateNavigator(rootnav);
  delete tmp;
  rootnav->FindNextBoundaryAndStep(Math::Min<double>(1.E20, pstep), !isonbdr);
  double stepcmp = Math::Max<double>(2 * gTolerance, rootnav->GetStep());
  double safecmp = rootnav->GetSafeDistance();
  if (Math::Abs(step - stepcmp) > 1E-6) {
    Geant::Print("","## PSTEP %lf VECGEOMSTEP %lf ROOTSTEP %lf", pstep, step, stepcmp);
    Geant::Print("","## PSTEP %lf ONBOUND %d VECGEOMSAFETY %lf ROOTSAFETY %lf BRUTEFORCEROOT %lf", pstep, isonbdr,
           safe, rootnav->Safety());

    // check nextpath
    tmp = outstate->ToTGeoBranchArray();
    tmp->InitFromNavigator(rootnav);
    Geant::Print("","## VECGEOMNEXTNODE %p ROOTNEXTNODE %p", outstate->GetCurrentNode(), tmp->GetCurrentNode());
    Geant::Print("","## VECGEOMBOUNDARY %d ROOTBOUNDARY %d", outstate->IsOnBoundary(), rootnav->IsOnBoundary());

    Geant::Print("","INCONSISTENT STEP");
    nav.InspectEnvironmentForPointAndDirection(Vector3D_t(x, y, z) /*global pos*/,
                                               Vector3D_t(dirx, diry, dirz), *instate);
  }
#endif // CROSSCHECK

#ifdef VERBOSE
  Geant::Print("","navfindbound on %p with pstep %lf yields step %lf and safety %lf\n", this, pstep, step,
         safe);
#endif // VERBOSE
}

//______________________________________________________________________________
void ScalarNavInterfaceVG::NavIsSameLocation(int ntracks,
       const double *x, const double *y, const double *z,
       const double */*dirx*/, const double */*diry*/, const double */*dirz*/,
       const VolumePath_t **start, VolumePath_t **end, bool *same, VolumePath_t *tmpstate) {
// 
// Checks if the navigation states corresponding to positions (x,y,z) are the
// same as the ones pointed by start. Update new states in end.
// Input:  ntracks - number of tracks to be checked
//         x,y,z   - arrays of positions
//         start   - starting navigation paths to compare with
// Output: end     - navigation paths corresponding to given positions
//         same    - flags showing if the end and start positions are matching
  
  typedef Vector3D<Precision> Vector3D_t;

  SimpleNavigator nav;
  for (int itr = 0; itr < ntracks; ++itr) {

// cross check with answer from ROOT
#ifndef GEANT_NVCC
#ifdef CROSSCHECK
    TGeoBranchArray *sb = start[itr]->ToTGeoBranchArray();
    TGeoBranchArray *eb = end[itr]->ToTGeoBranchArray();
#endif
#endif

    // TODO: not using the direction yet here !!
    bool samepath = nav.HasSamePath(
      Vector3D_t(x[itr], y[itr], z[itr]), *start[itr], *tmpstate);
    if (!samepath) tmpstate->CopyTo(end[itr]);
#ifndef GEANT_NVCC
#ifdef CROSSCHECK
    TGeoNavigator *nav = gGeoManager->GetCurrentNavigator();
    nav->ResetState();
    nav->SetLastSafetyForPoint(0, 0, 0, 0);
    nav->SetCurrentPoint(x[itr], y[itr], z[itr]);
    sb->UpdateNavigator(nav);
    bool rootsame = nav->IsSameLocation(x[itr], y[itr], z[itr], true);
    if (rootsame != samepath) {
      Geant::Print("","INCONSISTENT ANSWER ROOT(%d) VECGEOM(%d)", rootsame, samepath);
      std::cout << Vector3D_t(x[itr], y[itr], z[itr]) << "\n";
      Geant::Print("","old state");
      sb->Print();
      nav->ResetState();
      nav->SetLastSafetyForPoint(0, 0, 0, 0);
      nav->SetCurrentPoint(x[itr], y[itr], z[itr]);
//      nav->SetCurrentDirection(xdir[itr], ydir[itr], zdir[itr]);
      sb->UpdateNavigator(nav);
      nav->InspectState();
      bool rootsame = nav->IsSameLocation(x[itr], y[itr], z[itr], true);
      nav->InspectState();
      eb->InitFromNavigator(nav);
      Geant::Print("","new state");
      eb->Print();
      Geant::Print("","VERSUS VECGEOM OLD AND NEW");
      start[itr]->printVolumePath();
      end[itr]->printVolumePath();
    }
    delete sb;
    delete eb;
#endif // CROSSCHECK
#endif 
    same[itr] = samepath;
  }
}

//______________________________________________________________________________
void ScalarNavInterfaceVG::NavIsSameLocation(const double &x, const double &y, const double &z,
       const double &/*dirx*/, const double &/*diry*/, const double &/*dirz*/,
       const VolumePath_t *start, VolumePath_t *end, bool &same, VolumePath_t *tmpstate) {
// 
// Checks if the navigation states corresponding to positions (x,y,z) are the
// same as the ones pointed by start. Update new states in end.
// Input:  ntracks - number of tracks to be checked
//         x,y,z   - arrays of positions
//         start   - starting navigation paths to compare with
// Output: end     - navigation paths corresponding to given positions
//         same    - flags showing if the end and start positions are matching
  
  typedef Vector3D<Precision> Vector3D_t;

  SimpleNavigator nav;

// cross check with answer from ROOT
#ifdef CROSSCHECK
  TGeoBranchArray *sb = start->ToTGeoBranchArray();
  TGeoBranchArray *eb = end->ToTGeoBranchArray();
#endif

  // TODO: not using the direction yet here !!
  bool samepath = nav.HasSamePath(
    Vector3D_t(x, y, z), *start, *tmpstate);
  if (!samepath) tmpstate->CopyTo(end);

#ifdef CROSSCHECK
  TGeoNavigator *nav = gGeoManager->GetCurrentNavigator();
  nav->ResetState();
  nav->SetLastSafetyForPoint(0, 0, 0, 0);
  nav->SetCurrentPoint(x, y, z);
  sb->UpdateNavigator(nav);
  bool rootsame = nav->IsSameLocation(x, y, z, true);
  if (rootsame != samepath) {
    Geant::Print("","INCONSISTENT ANSWER ROOT(%d) VECGEOM(%d)", rootsame, samepath);
    std::cout << Vector3D_t(x, y, z) << "\n";
    Geant::Print("","old state");
    sb->Print();
    nav->ResetState();
    nav->SetLastSafetyForPoint(0, 0, 0, 0);
    nav->SetCurrentPoint(x, y, z);
//      nav->SetCurrentDirection(xdir, ydir, zdir);
    sb->UpdateNavigator(nav);
    nav->InspectState();
    bool rootsame = nav->IsSameLocation(x, y, z, true);
    nav->InspectState();
    eb->InitFromNavigator(nav);
    Geant::Print("","new state");
    eb->Print();
    Geant::Print("","VERSUS VECGEOM OLD AND NEW");
    start->printVolumePath();
    end->printVolumePath();
  }
  delete sb;
  delete eb;
#endif // CROSSCHECK
  same = samepath;
}

} // GEANT_IMPL_NAMESPACE
} // Geant
