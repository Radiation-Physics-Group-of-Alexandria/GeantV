
#include "PhotoElectricProcess.h"
#include "ComptonKleinNishina.h"
#include "GammaModelWrapper.h"

namespace geantphysics {

PhotoElectricProcess::PhotoElectricProcess(const std::string &name)
: EMPhysicsProcess(name) {
  // set process type to be an EM process.  Note: i.e. not relevant for E-loss
  SetType(ProcessType::kElectromagnetic);
  // Set it as a discrete process
  SetIsContinuous(false);
  SetIsDiscrete(true);
  // fill the list of particles that this process can be used to (note: models need to set either to be for e- or e+)
  AddToListParticlesAlloedToAssigned(Gamma::Definition());

  // ComptonKleinNishina(Random_t *states = 0, int threadId = -1);    
  auto comptonKNvec= new vecphys::ComptonKleinNishina();

  EMModel *comptonModel=
     new GammaModelWrapper<vecphys::ComptonKleinNishina>(
        "Compton Klein Nishina model (Wrapped VecPhys model)",
        comptonKNvec);
  AddModel( comptonModel );
}

} // namespace geantphysics
