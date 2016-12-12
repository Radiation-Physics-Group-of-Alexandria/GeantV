#include "GammaComptonProcess.h"

namespace geantphysics {

GammaComptonProcess::GammaComptonProcess(const std::string &name)
: EMPhysicsProcess(name) {
  // set process type to be an energy loss process (note: loss tables will be built automatically)
  SetType(ProcessType::kElectromagnetic); //mb: compton
  // set to be a continuous-discrete process
  //SetIsContinuous(true);
  SetIsDiscrete(true);
  // fill the list of particles that this process can be used to (note: models need to set either to be for e- or e+)
  //AddToListParticlesAlloedToAssigned(Electron::Definition());//mb: commented out
  //AddToListParticlesAlloedToAssigned(Positron::Definition());//mb: commented out
  AddToListParticlesAlloedToAssigned(Gamma::Definition());//mb: added
}

} // namespace geantphysics
