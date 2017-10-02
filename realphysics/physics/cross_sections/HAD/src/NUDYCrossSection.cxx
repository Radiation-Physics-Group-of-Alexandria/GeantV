//============================================================
// @file NUDYCrossSection.cxx
// @author Abhijit Bhattacharyya
// @brief This file makes and interface of NUDY to GV
//============================================================

#include <iostream>
#include "NUDYCrossSection.h"


using namespace geantphysics;
using namespace Nudy;
using namespace NudyPhysics;

// namespace geantphysics {

NUDYCrossSection::NUDYCrossSection(): fisoName("Fe"), ftargetN(56), ftargetZ(26),
				      fprojEnergy("4.0") {}

NUDYCrossSection::NUDYCrossSection( const std::string isoName, const int Nval, const int Zval,
				    const double projectileKE, const HadronicProcessType IntType ) :
  fisoName(isoName), ftargetN(Nval), ftargetZ(Zval), fprojEnergy(projectileKE) {}

NUDYCrossSection::~NUDYCrossSection() {}

bool NUDYCrossSection::IsApplicable(const int projCode, cost double projKE ) {
  bool isOK = false;
  if (projCode == 2112 && projectileKE < 20.0 * geant::MeV) isOK = true;  // neutron 3 || 2112 ? and KE < 20 MeV
  return isOK;
}

double GetIsotopeCrossSection( const int projCode, const double projectileKE, const int Z, const int N) {

  double csValue = 0.0;

  std::string isoName = GetisoName();    //  geantphysics::Material::GetName();
  // NudyInterface::GetXS(projCode, projectileKE, isoName, Z, N);
  

  return csValue;
}

// }  // namespace
