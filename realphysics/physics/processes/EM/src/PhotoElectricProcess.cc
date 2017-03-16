//
//  J. Apostolakis & S.Y. Jun,  March 2017
//
#include "PhotoElectricProcess.h"
#include "PhotoElectronSauterGavrila.h"
#include "GammaModelWrapper.h"

namespace geantphysics {

PhotoElectricProcess::PhotoElectricProcess(const std::string &name)
: EMPhysicsProcess(name) {
  using vecphys::VECPHYS_IMPL_NAMESPACE::PhotoElectronSauterGavrila;

  // set process type to be an EM process.  Note: i.e. not relevant for E-loss
  SetType(ProcessType::kElectromagnetic);
  // Set it as a discrete process
  SetIsContinuous(false);
  SetIsDiscrete(true);
  // fill the list of particles that this process can be used to (note: models need to set either to be for e- or e+)
  AddToListParticlesAlloedToAssigned(Gamma::Definition());

  // PhotoElectricModel(Random_t *states = 0, int threadId = -1);    
  auto vecphysPhotoElectricSG= new PhotoElectronSauterGavrila(nullptr);
  double lowE_inKeV =    1.0;
  double highE_inMeV = 100.0;
  
  vecphysPhotoElectricSG->SetLowEnergyLimit (  lowE_inKeV * CLHEP::keV);  // geant::keV);
  vecphysPhotoElectricSG->SetHighEnergyLimit( highE_inMeV * CLHEP::MeV);  // geant::MeV);
  
  EMModel *photoElectricModel=
     new GammaModelWrapper<PhotoElectronSauterGavrila>(
        "PhotoElectronSauterGavrila model (Wrapped VecPhys model)",
        vecphysPhotoElectricSG);
  photoElectricModel->SetLowEnergyUsageLimit (  lowE_inKeV * geant::keV);
  photoElectricModel->SetHighEnergyUsageLimit( highE_inMeV * geant::MeV);
  
  AddModel( photoElectricModel );
}

} // namespace geantphysics
