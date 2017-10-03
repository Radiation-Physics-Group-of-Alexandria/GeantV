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

  double nXS;
  std::string eleName = "Fe";

  nXS = nudyxs.GetXS(2112, 1.0, eleName, 26, 56);


  std::cout << nXS << std::endl;



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
