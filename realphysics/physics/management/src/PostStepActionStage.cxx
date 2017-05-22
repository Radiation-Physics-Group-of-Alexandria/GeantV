
#include "PostStepActionStage.h"

// from geantV
#include "GeantPropagator.h"
#include "GeantTaskData.h"
#include "GeantTrack.h"
#include "Handler.h"

// from realphysics
#include "Material.h"
#include "MaterialCuts.h"
#include "Region.h"

#include "PhysicsProcess.h"
#include "PhysicsManagerPerParticle.h"
#include "LightTrack.h"

// handler(s)
#include "PostStepActionHandler.h"


namespace geantphysics {

PostStepActionStage::PostStepActionStage(Geant::GeantPropagator *prop)
: SimulationStage(Geant::kPostStepActionStage, prop) {}


// base class will delete the created handlers
PostStepActionStage::~PostStepActionStage() {}


int PostStepActionStage::CreateHandlers() {
  int threshold = fPropagator->fConfig->fNperBasket;
  // create the only one handler
  AddHandler(new PostStepActionHandler(threshold, fPropagator));
  // only one handler is created
  return 1;
}


// Selects tracks that have any processes, any post step processes i.e. discrete part and that limited the step
Geant::Handler* PostStepActionStage::Select(Geant::GeantTrack *track, Geant::GeantTaskData * /*td*/) {
  if (track->fStatus==Geant::TrackStatus_t::kPhysics && track->fEindex==1000) {
    // these tracks should always have psorcesses active in the given region moreover should always have discrete
    // processes that limited the step (fEindex==1000)
    //
    // reset number of interaction length left
    track->fNintLen = -1;
    return fHandlers[0];
    /*
    //
    // here we will get the MaterialCuts from the LogicalVolume later
    int   matIndx              = track->GetMaterial()->GetIndex();
    int   regIndx              = const_cast<vecgeom::LogicalVolume*>(track->GetVolume())->GetRegion()->GetIndex();
    const MaterialCuts *matCut =  MaterialCuts::GetMaterialCut(regIndx,matIndx);
    // get the internal code of the particle
    int   particleCode         = track->GVcode();
    const Particle *particle   = Particle::GetParticleByInternalCode(particleCode);
    // get the PhysicsManagerPerParticle for this particle: will be nullptr if the particle has no any PhysicsProcess-es
    PhysicsManagerPerParticle *pManager = particle->GetPhysicsManagerPerParticlePerRegion(matCut->GetRegionIndex());
    // here we could drop the
    if (pManager && pManager->GetListPostStepCandidateProcesses().size()>0) {
      // give back the only one handler for all other particles that will compute the along-step-action(s)
      return fHandlers[0];
    }
    */
  }
  // not physics or not discrete part of limited the step
  return nullptr;
}



}  // namespace geantphysics
