
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
  auto photoElectricModelSG= new PhotoElectronSauterGavrila(nullptr);

  EMModel *comptonModel=
     new GammaModelWrapper<PhotoElectronSauterGavrila>(
        "PhotoElectronSautherGavrila model (Wrapped VecPhys model)",
        photoElectricModelSG);
  AddModel( comptonModel );
}

} // namespace geantphysics
