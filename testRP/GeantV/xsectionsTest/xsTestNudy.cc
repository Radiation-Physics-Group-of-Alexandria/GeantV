#include "GlauberGribovTotalXsc.h"
#include "GlauberGribovInelasticXsc.h"
#include "GlauberGribovElasticXsc.h"
#include "Particle.h"
#include <iostream>
#include <fstream>
#include <sstream>
#include <string>
#include "Proton.h"
#include "Neutron.h"
#include "KaonMinus.h"
#include "Electron.h"

#include "SystemOfUnits.h"
#include "PhysicalConstants.h"

#include "NudyInterface.h"

using geantphysics::Particle;
using geantphysics::Proton;
using geantphysics::Neutron;
using geantphysics::KaonMinus;

using namespace Nudy;
using NudyPhysics::NudyInterface;


int main(int /*argc*/, char** /*argv*/) {

  NudyPhysics::NudyInterface nudyxs;


  std::ofstream writef("sumXsc.dat", std::ios::out ) ;

  writef.setf( std::ios::scientific, std::ios::floatfield );

  double nxST; // total
  double nxSE; // elastic
  double nxSI; // inelastic
  double nxSF; // fission
  double nxSH; // thermal
  int projectileCode = 2112; // examples
  std::string eleName = "Pu";
  std::string reactType = "Fission"; // "Elastic, Enelastic, total, fission, thermal......"
  int Zvalue = 94;
  int Massvalue = 241;
  int Nvalue = Massvalue - Zvalue;
  double temperature = 293.60608;
  /*  This is showing energy as 0.004 and not as 4.0E=6 and so
  for the time being I am testing using raw number.
  double EnergyValue = 4.0 * geant::MeV;
  */
  double EnergyValue = 1.0e+6;

// @brief Here we provide code for projectile say 2112 for neutron, energy of projectile say 1.0 MeV
// @brief Then we provide temperature which is required for fission etc.
// @brief Then we provide element name like Pu , U etc. for which cross section is required
// @brief Z for element name
// @brief N for neutron number for the element.

nxSF = nudyxs.GetNudyXS(projectileCode, EnergyValue, temperature, eleName, Zvalue, Nvalue, "Fission");
std::cout << "Fission CrossSection (NuT) at En : " << EnergyValue << " eV = " <<  nxSF << std::endl;

// commented for testing one by one
/*
  nxST = nudyxs.GetNudyXS(projectileCode, EnergyValue, temperature, eleName, Zvalue, Nvalue, "Total");
  std::cout << "Total - " << nxST << std::endl;

  nxSE = nudyxs.GetNudyXS(projectileCode, EnergyValue, temperature, eleName, Zvalue, Nvalue, "Elastic");
  std::cout << "Elastic - " << nxSE << std::endl;

  nxSI = nudyxs.GetNudyXS(projectileCode, EnergyValue, temperature, eleName, Zvalue, Nvalue, "Inelastic");
  std::cout << "Inelastic - " << nxSI << std::endl;

  nxSH = nudyxs.GetNudyXS(projectileCode, EnergyValue, temperature, eleName, Zvalue, Nvalue, "Thermal");
  std::cout << "Thermal - " << nxSH << std::endl;
*/


/*
  double kinEnergy = 10.* geant::MeV;
  double maxEnergy = 1000 * geant::GeV;
  int particlePDG = -321;
  //int particlePDG = 2212;
  int Z = 82;
  // number of nucleons
  int N =207;
  double txs, ixs, exs;

  while (kinEnergy < maxEnergy)
    {

    txs = totxs.GetIsotopeCrossSection(geantphysics::Particle::GetParticleByPDGCode(particlePDG)->GetInternalCode(),
			       kinEnergy, geantphysics::Particle::GetParticleByPDGCode(particlePDG)->GetPDGMass(), Z, N) / geant::millibarn;

    ixs = inexs.GetIsotopeCrossSection(geantphysics::Particle::GetParticleByPDGCode(particlePDG)->GetInternalCode(),
			       kinEnergy, geantphysics::Particle::GetParticleByPDGCode(particlePDG)->GetPDGMass(), Z, N) / geant::millibarn;

    exs = elaxs.GetIsotopeCrossSection(geantphysics::Particle::GetParticleByPDGCode(particlePDG)->GetInternalCode(),
			       kinEnergy, geantphysics::Particle::GetParticleByPDGCode(particlePDG)->GetPDGMass(), Z, N) / geant::millibarn;

    std::cout <<"energyKin " << kinEnergy <<  " total " << txs << " elastic " << exs << " inelastic " << ixs <<  std::endl;

    writef << kinEnergy/geant::MeV <<"  \t" << ixs << "  \t"
     << exs << "  \t"<< 0 << "  \t"<< txs << std::endl;

    kinEnergy *= 1.001;
  }
*/
  return 0;
}
