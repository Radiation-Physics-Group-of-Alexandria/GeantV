#include "PionMinus.h"

#include "PhysicalConstants.h"
#include "SystemOfUnits.h"

namespace geantphysics {

PionMinus* PionMinus::Definition() {
  static PionMinus instance("pi-", -211, 11, 0.13957*geant::GeV, -1.0*geant::eplus); // mass value taken from Geant4 10.3
  return &instance;
}

} // namespace geantphysics
