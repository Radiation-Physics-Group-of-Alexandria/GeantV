#ifndef GEANT_VECPHYSORCHESTRATOR
#define GEANT_VECPHYSORCHESTRATOR

#include "Geant/Config.h"
#include "base/Global.h"
#include "Geant/Typedefs.h"

#include "GeantFwd.h"
#include "GUTrack.h"
#include "GeantTrack.h"

#include "ComptonKleinNishina.h"


//class GUTrack_v;
class VecPhysOrchestrator{
public:
    
    using GeantTrack = Geant::GeantTrack;
    using GeantTrack_v = Geant::GeantTrack_v;
    
    VecPhysOrchestrator();
    VecPhysOrchestrator(int processId){fProcessId=processId;};
    VecPhysOrchestrator(int processId, double energyLimit, int maxnumtracks);
    ~VecPhysOrchestrator();

    void Allocator(int size);
    void Deallocator();
    void GUTrackAllocator(GUTrack_v &gutrack_v, int size);
    void GUTrackDeallocator(GUTrack_v &gutrack_v);
    
    void PerformInteraction();
    void FilterPrimaryTracks(GeantTrack_v &gTrackV, int numtracks);
    void FilterTracksForTabPhys(GeantTrack_v &gTrackV, GeantTrack_v &gTabulatedPhysicsTracks, int numtracksIn, int &numtracksOut, int  *parentTrackIndices);
    int  WriteBackTracks(GeantTrack_v& gTrackV, int tid);
    void SetGeantTrack(GeantTrack &left, GeantTrack_v &right, int ip);
    int ApplyPostStepProcess(GeantTrack_v &gTrackV, int numtracks, int tid);
    void ConvertEnergiesToVecPhys();
    void ConvertEnergiesFromVecPhys();
    void DebugTracksEnergies(GeantTrack_v &gTrackV, int numtracks, GUTrack_v &primaries, bool checkIdentity);
    void DebugPrimaryAndSecondaryTracks();
    void CheckEnergyConservation(GeantTrack_v &gTrackV, int numtracks, GUTrack_v &primaries, GUTrack_v &secondaries);
    void CheckDirectionUnitVector(GeantTrack_v &gTrackV, GUTrack_v &primaries, GUTrack_v &secondaries);
    void RotateNewTrack(double oldXdir, double oldYdir, double oldZdir, GeantTrack_v &track, int index);
    int     fComptonTotTracks; //temporary
private:
    int     fProcessId;       // Index of the physics process as stored in tabulated physics (TPartIndex)
    double  fEnergyLimit;     // Tracking cut in kinetic energy (In GeantV units [GeV])
    
    vecphys::ComptonKleinNishina*  fVComptonModel;  // Compton Vector physics model
    
    // Intermediate containers
    int     *fTargetElements;       // Z-of target atoms
    int     *fParentTrackIndices;   // indices of the parent tracks in GeantTrack_v
    GUTrack_v *fPrimaryTracks;      // GUTrack_v with the proper primary tracks for Compton
    GUTrack_v *fSecondaryTracks;    // GUTrack_v with the corresponding secondary tracks
    GUTrack_v *fTabulatedPhysicsTracks; // GUTrack_v with the tracks for TabulatedPhysics
    
    
};

#endif
