/**********************************************
 * @file NudyInterface.cxx
 * @author Abhijit Bhattacharyya
 * @brief This interface links to the GV through NUDYCrossSection of realphysics/physics/HAD/CrossSection
 **********************************************/
#include <iostream>
#include "NudyInterface.h"

using namespace geantphysics;
using namespace Nudy;
using namespace NudyPhysics;

NudyInterface::NudyInterface() :
 fProjCode(2112), fProjKE(1.0*geant::MeV), fTemperature(293.60608), fIsoName("Fe"), ftZ(26), ftN(56), fReactType("Total"){};

NudyInterface::NudyInterface(
  const int projCode, const double projKE, double temp, const std::string isoName, const int tZ, const int tN, const std::string reactType)
  // : fProjCode(projCode), fProjKE(projKE), fIsoName(isoName), ftZ(tZ), ftN(tN) {
  {
    std:;string newReact;
    if (!reactType.length()) {
      newReact = "Total";
    }else {
      newReact = reactType;
    }
    SetProjectileCode (projCode);
    SetProjectileKE (projKE);
    SetTemp (temp);
    SetIsotopeName (isoName);
    SetZ (tZ);
    SetN (tN);
    SetReactionType(newReact);
  };

NudyInterface::~NudyInterface() {}

double NudyInterface::GetNudyXS( int projCode, double projKE, double temp, std::string isoName, int tZ, int tN, std::string reactType ) {
  NudyPhysics::NudyXSProcess xsProc;

  xsProc.GetXS( projCode, projKE, temp, isoName, tZ, tN, reactType) ;
  return 0;
}
