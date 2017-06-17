
#include "MSCProcess.h"

#include "SystemOfUnits.h"

#include "Electron.h"
#include "Positron.h"

#include "Region.h"

#include "MSCModel.h"
#include "EMModelManager.h"

// from geantV
#include "GeantTaskData.h"
#include "GeantTrack.h"

#include "Geant/NavigationInterface.h"

namespace geantphysics {

MSCProcess::MSCProcess(const std::string &name) : EMPhysicsProcess(name), fGeomMinLimit(0.05*geant::nm) {
  fGeomMinLimit2 = fGeomMinLimit*fGeomMinLimit;
  // process type is kElectromagnetic in the base EMPhysicsProcess calss so set it to kMSC
  SetType(ProcessType::kMSC);
  // set to be a pure continuous process
  SetIsContinuous(true);
  SetIsDiscrete(true); // just for the framework: msc has no discrete part (at the moment)
  // fill the list of particles that this process can be used to i.e. gamma particle
  AddToListParticlesAlloedToAssigned(Electron::Definition());
  AddToListParticlesAlloedToAssigned(Positron::Definition());
}

MSCProcess::~MSCProcess() {}

// called at the PrePropagationStage(in the Handler)
void MSCProcess::AlongStepLimitationLength(Geant::GeantTrack *gtrack, Geant::GeantTaskData *td) const {
  // init all lengths to the current minimum physics step length (that is the true length)
  //
  bool isOnBoundaryPostStp = gtrack->fBoundary;
  gtrack->fBoundary        = gtrack->fIsOnBoundaryPreStp;
  // get the phyics step limit due to all (other than msc) physics processes (that was used in the geometry stage)
  double minPhysicsStepLength   = gtrack->fPstep;
  gtrack->fTheTrueStepLenght    = minPhysicsStepLength;
  gtrack->fTheTransportDistance = minPhysicsStepLength;
  gtrack->fTheZPathLenght       = minPhysicsStepLength;
  gtrack->SetDisplacement(0.,0.,0.);
  gtrack->SetNewDirectionMsc(0.,0.,1.);
  // select msc model
  double ekin        = gtrack->fE-gtrack->fMass;
  int    regIndx     = const_cast<vecgeom::LogicalVolume*>(gtrack->GetVolume())->GetRegion()->GetIndex();
  MSCModel *mscModel = static_cast<MSCModel*>(GetModelManager()->SelectModel(ekin,regIndx));
  // check if:
  // - current true step length is > limit
  // - current kinetic energy is above the minimum
  if (mscModel && minPhysicsStepLength>GetGeomMinLimit()) {
    // will check min/max usage limits, possible additional step limit, update the
    // gtrack->fTheTrueStepLenght and will convert the true->to->geometic length (that will be in gtrack->fTheZPathLenght)
    mscModel->StepLimit(gtrack, td);
    // check if msc limited the step: set that continuous process was the winer
    if (gtrack->fTheTrueStepLenght<minPhysicsStepLength) {
      gtrack->fEindex = -1; // indicate that continuous process was the winer
    }
  }
  // update gtrack->fPstep to be the geometric step length (this is what propagation needs):
  // straight line discrete along the original direction
  // protection: geometric legth must always be <= than the true physics step length
  gtrack->fTheZPathLenght = std::min(gtrack->fTheZPathLenght, minPhysicsStepLength);
  gtrack->fPstep          = gtrack->fTheZPathLenght;
  // check if the geometrical physics step (true is change now by MSC to geometric) become shorter than snext and limit
  // snext to this geometrical physics step (track will be propagated to snext distance) and set the post-step point
  // boundary flag to false (we pulled back)
  if (gtrack->fSnext>gtrack->fPstep) {
    gtrack->fSnext    = gtrack->fPstep;
    gtrack->fBoundary = false;
  } else { // write back the original post-step point boundary flag and do nothing
    gtrack->fBoundary = isOnBoundaryPostStp;
  }
}

// called at the PostPropagationStage(in the Handler)
void MSCProcess::AlongStepDoIt(Geant::GeantTrack *gtrack, Geant::GeantTaskData *td) const {
  // get the step length (the geometric one)
  double geometricStepLength = gtrack->fStep;
  double truePathLength      = geometricStepLength;
  // select msc model
  double ekin        = gtrack->fE-gtrack->fMass;
  int    regIndx     = const_cast<vecgeom::LogicalVolume*>(gtrack->GetVolume())->GetRegion()->GetIndex();
  MSCModel *mscModel = static_cast<MSCModel*>(GetModelManager()->SelectModel(ekin,regIndx));
  // check if:
  // - current true step length is > limit
  // - current kinetic energy is above the minimum usea
  if (mscModel && gtrack->fTheTrueStepLenght>GetGeomMinLimit()) {
    truePathLength = gtrack->fTheTrueStepLenght;
    // gtrack->fTheTrueStepLenght be updated during the conversion
    mscModel->ConvertGeometricToTrueLength(gtrack, td);
    // protection againts wrong true-geometic-true gonversion
    truePathLength = std::min(truePathLength,gtrack->fTheTrueStepLenght);
    // optimization: scattring is not sampled if the particle is reanged out in this step or short step
    if (gtrack->fRange>truePathLength && truePathLength>GetGeomMinLimit()) {
      // sample scattering: might have been done during the step limit phase
      // NOTE: in G4 the SampleScattering method won't use the possible shrinked truePathLength!!! but make it correct
      gtrack->fTheTrueStepLenght = truePathLength;
      bool hasNewDir = mscModel->SampleScattering(gtrack, td);
      // compute displacement vector length
      double dl =   gtrack->fTheDisplacementVectorX*gtrack->fTheDisplacementVectorX
                  + gtrack->fTheDisplacementVectorY*gtrack->fTheDisplacementVectorY
                  + gtrack->fTheDisplacementVectorZ*gtrack->fTheDisplacementVectorZ;
      if (dl>fGeomMinLimit2 && !gtrack->fBoundary) {
        // displace the post-step point
        dl = std::sqrt(dl);
        double dir[3]={gtrack->fTheDisplacementVectorX/dl, gtrack->fTheDisplacementVectorY/dl, gtrack->fTheDisplacementVectorZ/dl};
        ScalarNavInterface::DisplaceTrack(*gtrack,dir,dl,GetGeomMinLimit());
      }
      // apply msc agular deflection
//      if (!gtrack->fBoundary && hasNewDir) {
      if (hasNewDir) {
        gtrack->fXdir = gtrack->fTheNewDirectionX;
        gtrack->fYdir = gtrack->fTheNewDirectionY;
        gtrack->fZdir = gtrack->fTheNewDirectionZ;
      }
    }
  }
  // update step length to store the true step length
  gtrack->fStep = truePathLength;
}



}  // namespace geantphysics
