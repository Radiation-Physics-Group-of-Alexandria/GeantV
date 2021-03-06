#include "KaonZero.h"

#include "PhysicalConstants.h"
#include "SystemOfUnits.h"

namespace geantphysics {

KaonZero* KaonZero::Definition() {
  static KaonZero instance("K0", 311, 17, 0.497614*geant::GeV, 0); // mass value taken from Geant4 10.3
  return &instance;
}

} // namespace geantphysics
