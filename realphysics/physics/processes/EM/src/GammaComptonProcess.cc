#include "GammaComptonProcess.h"

namespace geantphysics {

GammaComptonProcess::GammaComptonProcess(const std::string &name)
: EMPhysicsProcess(name) {
  // set process type to be a discrete EM process
  SetType(ProcessType::kElectromagnetic);
  // set to be a discrete process
  SetIsDiscrete(true);
  // fill the list of particles that this process can be used to
  AddToListParticlesAlloedToAssigned(Gamma::Definition());
}

} // namespace geantphysics
