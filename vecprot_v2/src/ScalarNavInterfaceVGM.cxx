#include "ScalarNavInterfaceVGM.h"

#include "backend/Backend.h"
#include "navigation/VNavigator.h"
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

#ifdef CROSSCHECK
#include "TGeoNavigator.h"
#include "TGeoNode.h"
#endif

#ifdef BUG_HUNT
#include "GeantPropagator.h"
#endif

namespace Geant {
inline namespace GEANT_IMPL_NAMESPACE {

using namespace VECGEOM_NAMESPACE;

//______________________________________________________________________________
GEANT_CUDA_BOTH_CODE
void ScalarNavInterfaceVGM::NavFindNextBoundaryAndStep(int ntracks, const double *pstep, 
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
    // If the basket is mixed volumes/navigators may be different
    VNavigator const * newnav = instate[itr]->Top()->GetLogicalVolume()->GetNavigator();
    // Check if current safety allows for the proposed step
    if (safe[itr] > pstep[itr]) {
      step[itr] = pstep[itr];
      isonbdr[itr] = false;
      *outstate[itr] = *instate[itr];
      continue;
    }

    step[itr] = newnav->ComputeStepAndSafetyAndPropagatedState(Vector3D_t(x[itr], y[itr], z[itr]),
                               Vector3D_t(dirx[itr], diry[itr], dirz[itr]),
                               Math::Min<double>(1.E20, pstep[itr]), *instate[itr], *outstate[itr] /* the paths */, !isonbdr[itr], safe[itr]);
    step[itr] = Math::Max<double>(2. * gTolerance, step[itr] + 2. * gTolerance);
    safe[itr] = Math::Max<double>(safe[itr], 0);
    // onboundary with respect to new point
    isonbdr[itr] = outstate[itr]->IsOnBoundary();
    
    //#### To add small step detection and correction - see ScalarNavInterfaceTGeo ####//

#if defined(CROSSCHECK) && !defined(GEANT_NVCC)
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

#ifdef VERBOSE
    Geant::Print("","navfindbound on %p track %d with pstep %lf yields step %lf and safety %lf\n", this, itr, pstep[itr], step[itr],
           safe[itr]);
#endif // VERBOSE
  }
}

//______________________________________________________________________________
GEANT_CUDA_BOTH_CODE
void ScalarNavInterfaceVGM::NavFindNextBoundaryAndStep(GeantTrack &track) {

  typedef Vector3D<Precision> Vector3D_t;
#ifdef GEANT_NVCC
  SimpleNavigator nav;
#else
  ABBoxNavigator nav;
#endif // GEANT_NVCC

  // Retrieve navigator for the track
  VNavigator const * newnav = track.fPath->Top()->GetLogicalVolume()->GetNavigator();
  // Check if current safety allows for the proposed step
  if (track.fSafety > track.fPstep) {
    track.fStep = track.fPstep;
    track.fBoundary = false;
    *track.fNextpath = *track.fPath;
    return;
  }

  track.fStep = newnav->ComputeStepAndSafetyAndPropagatedState(Vector3D_t(track.fXpos, track.fYpos, track.fZpos),
                             Vector3D_t(track.fXdir, track.fYdir, track.fZdir),
                             Math::Min<double>(1.E20, track.fPstep),
                             *track.fPath, *track.fNextpath, !track.fBoundary, track.fSafety);
  track.fStep = Math::Max<double>(2. * gTolerance, track.fStep + 2. * gTolerance);
  track.fSafety = Math::Max<double>(track.fSafety, 0);
  // onboundary with respect to new point
  track.fBoundary = track.fNextpath->IsOnBoundary();

  //#### To add small step detection and correction - see ScalarNavInterfaceTGeo ####//

#ifdef CROSSCHECK
  //************
  // CROSS CHECK USING TGEO
  //************
  TGeoNavigator *rootnav = gGeoManager->GetCurrentNavigator();
  rootnav->ResetState();
  rootnav->SetCurrentPoint(track.fXpos, track.fYpos, track.fZpos);
  rootnav->SetCurrentDirection(track.fXdir, track.fYdir, track.fZdir);
  TGeoBranchArray *tmp = track.fPath->ToTGeoBranchArray();
  tmp->UpdateNavigator(rootnav);
  delete tmp;
  rootnav->FindNextBoundaryAndStep(Math::Min<double>(1.E20, track.fPstep), !track.fBoundary);
  double stepcmp = Math::Max<double>(2 * gTolerance, rootnav->GetStep());
  double safecmp = rootnav->GetSafeDistance();
  if (Math::Abs(track.fStep - stepcmp) > 1E-6) {
    Geant::Print("","## PSTEP %lf VECGEOMSTEP %lf ROOTSTEP %lf", track.fPstep, track.fStep, stepcmp);
    Geant::Print("","## PSTEP %lf ONBOUND %d VECGEOMSAFETY %lf ROOTSAFETY %lf BRUTEFORCEROOT %lf", track.fPstep, track.fBoundary,
           track.fSafety, rootnav->Safety());

    // check nextpath
    tmp = track.fNextpath->ToTGeoBranchArray();
    tmp->InitFromNavigator(rootnav);
    Geant::Print("","## VECGEOMNEXTNODE %p ROOTNEXTNODE %p", track.fNextpath->GetCurrentNode(), tmp->GetCurrentNode());
    Geant::Print("","## VECGEOMBOUNDARY %d ROOTBOUNDARY %d", track.fNextpath->IsOnBoundary(), rootnav->IsOnBoundary());

    Geant::Print("","INCONSISTENT STEP");
    nav.InspectEnvironmentForPointAndDirection(Vector3D_t(track.fXpos, track.fYpos, track.fZpos) /*global pos*/,
                                               Vector3D_t(track.fXdir, track.fYdir, track.fZdir),
                                               *track.fPath);
  }
#endif // CROSSCHECK

#ifdef VERBOSE
  Geant::Print("","navfindbound on %p with pstep %lf yields step %lf and safety %lf\n", 
               this, track.fPtrack.fStep, track.fStep, track.fSafety);
#endif // VERBOSE
}

//______________________________________________________________________________
void ScalarNavInterfaceVGM::NavIsSameLocation(int ntracks,
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
  

  //#### NOT USING YET THE NEW NAVIGATORS ####//

  typedef Vector3D<Precision> Vector3D_t;

  SimpleNavigator nav;
  for (int itr = 0; itr < ntracks; ++itr) {

// cross check with answer from ROOT

    // TODO: not using the direction yet here !!
    bool samepath = nav.HasSamePath(
      Vector3D_t(x[itr], y[itr], z[itr]), *start[itr], *tmpstate);
    if (!samepath) tmpstate->CopyTo(end[itr]);

#if defined(CROSSCHECK) && !defined(GEANT_NVCC)
    TGeoBranchArray *sb = start[itr]->ToTGeoBranchArray();
    TGeoBranchArray *eb = end[itr]->ToTGeoBranchArray();
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
#endif // CROSSCHECK && !GEANT_NVCC
    same[itr] = samepath;
  }
}

//______________________________________________________________________________
void ScalarNavInterfaceVGM::NavIsSameLocation(GeantTrack &track, bool &same, VolumePath_t *tmpstate) {

  //#### NOT USING YET THE NEW NAVIGATORS ####//

  typedef Vector3D<Precision> Vector3D_t;

  SimpleNavigator nav;
// cross check with answer from ROOT
#ifdef CROSSCHECK
  TGeoBranchArray *sb = track.fPath->ToTGeoBranchArray();
  TGeoBranchArray *eb = track.fNextpath->ToTGeoBranchArray();
#endif

  // TODO: not using the direction yet here !!
  bool samepath = nav.HasSamePath(
    Vector3D_t(track.fXpos, track.fYpos, track.fZpos), *track.fPath, *tmpstate);
  if (!samepath) tmpstate->CopyTo(track.fNextpath);

#ifdef CROSSCHECK
  TGeoNavigator *nav = gGeoManager->GetCurrentNavigator();
  nav->ResetState();
  nav->SetLastSafetyForPoint(0, 0, 0, 0);
  nav->SetCurrentPoint(track.fXpos, track.fYpos, track.fZpos);
  sb->UpdateNavigator(nav);
  bool rootsame = nav->IsSameLocation(track.fXpos, track.fYpos, track.fZpos, true);
  if (rootsame != samepath) {
    Geant::Print("","INCONSISTENT ANSWER ROOT(%d) VECGEOM(%d)", rootsame, samepath);
    std::cout << Vector3D_t(track.fXpos, track.fYpos, track.fZpos) << "\n";
    Geant::Print("","old state");
    sb->Print();
    nav->ResetState();
    nav->SetLastSafetyForPoint(0, 0, 0, 0);
    nav->SetCurrentPoint(track.fXpos, track.fYpos, track.fZpos);
//      nav->SetCurrentDirection(xdir, ydir, zdir);
    sb->UpdateNavigator(nav);
    nav->InspectState();
    bool rootsame = nav->IsSameLocation(track.fXpos, track.fYpos, track.fZpos, true);
    nav->InspectState();
    eb->InitFromNavigator(nav);
    Geant::Print("","new state");
    eb->Print();
    Geant::Print("","VERSUS VECGEOM OLD AND NEW");
    track.fPath->printVolumePath();
    track.fNextpath->printVolumePath();
  }
  delete sb;
  delete eb;
#endif // CROSSCHECK
  same = samepath;
}

} // GEANT_IMPL_NAMESPACE
} // Geant
