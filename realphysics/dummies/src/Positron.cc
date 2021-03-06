#include "Positron.h"

#include "PhysicalConstants.h"
#include "SystemOfUnits.h"

namespace geantphysics {

Positron* Positron::Definition() {
  static Positron instance("e+",   -11, 23, geant::kElectronMassC2,  1.0*geant::eplus);
  return &instance;
}

} // namespace geantphysics
