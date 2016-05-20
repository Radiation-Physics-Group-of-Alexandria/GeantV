/*
 Author: M. Bandieramonte
 */
#include "VecPhysCompton.h"

#include "GeantTaskData.h"
#include "GeantVApplication.h"
#include "WorkloadManager.h"
#include "globals.h"
#include "GeantTrack.h"

#ifdef USE_VECGEOM_NAVIGATOR
#include "navigation/NavigationState.h"
#else
#include "TGeoBranchArray.h"
#endif

#ifdef USE_ROOT
#include "Rtypes.h"
ClassImp( VecPhysCompton )
#endif


VecPhysCompton::VecPhysCompton() : PhysicsProcess() {}


void VecPhysCompton::Initialize() {
    std::cout << "VecPhysCompton::Initialize : Start" << std::endl;
    std::cout << "VecPhysCompton::Initialize : End" << std::endl;
}

void VecPhysCompton::Eloss( Material_t *mat, int ntracks, GeantTrack_v &tracks, int &nout, GeantTaskData *td ) {
    
    //  std::cout << "VecPhysCompton::Eloss : Start" << std::endl;
    
    for ( int i = 0; i < ntracks; i++ ) {
        tracks.fEdepV[i] = 10;
    }
    //  std::cout << "VecPhysCompton::Eloss : --- End ---" << std::endl;
}


void VecPhysCompton::ComputeIntLen( Material_t *mat, int ntracks, GeantTrack_v &tracks, double * /*lengths*/,
                          GeantTaskData *td ) {
    //  std::cout << "VecPhysCompton::ComputeIntLen : Start : ntracks=" << ntracks << std::endl;
    for ( int i = 0; i < ntracks; i++ ) {
        //    std::cout << " VecPhysCompton::ComputeIntLen : proposedStepLengths[" << i << "]=" << 10
        //             << " ; position=("  << tracks.fXposV[i] << "," << tracks.fYposV[i] << "," << tracks.fZposV[i]
        //             << ")" << std::endl;  // Debug
        tracks.fPstepV[i] = 10;
        tracks.fEindexV[i] = 1000;  // Forced to be always treated as continuous & discrete.
    }
    //  std::cout << "VecPhysCompton::ComputeIntLen : --- End ---" << std::endl;
}


void VecPhysCompton::PostStepTypeOfIntrActSampling( Material_t *mat, int ntracks, GeantTrack_v &tracks,
                                          GeantTaskData *td ) {
    //  std::cout << "VecPhysCompton::PostStepTypeOfIntrActSampling : Start & End" << std::endl;
}


void VecPhysCompton::PostStepFinalStateSampling( Material_t *mat, int ntracks, GeantTrack_v &tracks, int &nout,
                                       GeantTaskData *td ) {
    //  std::cout << "VecPhysCompton::PostStepTypeOfIntrActSampling : Start" << std::endl;
    //  std::cout << "VecPhysCompton::PostStepTypeOfIntrActSampling : --- End ---" << std::endl;
}

