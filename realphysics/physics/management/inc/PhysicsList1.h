#include "SystemOfUnits.h"

#include "PhysicsList.h"
#include "PhysicsProcess.h"

#include "Particle.h"

#include "Electron.h"
#include "Positron.h"
#include "Gamma.h"

#include "ElectronIonizationProcess.h"
#include "MollerBhabhaIonizationModel.h"
#include "ElectronBremsstrahlungProcess.h"
#include "SeltzerBergerBremsModel.h"
#include "RelativisticBremsModel.h"

#include "ComptonProcess.h"
#include "PhotoElectricProcess.h"
#include "ConversionProcess.h"

namespace geantphysics {

class PhysicsList1 : public PhysicsList {

public:
  PhysicsList1(const std::string &name) : PhysicsList(name) {}

  virtual void Initialize() {
    // get the partcile table and loop over it
    std::vector<Particle*> pTable = Particle::GetTheParticleTable();
    for (unsigned int i=0; i<pTable.size(); ++i) {
      Particle *particle = pTable[i];
// #if 0      
      if (particle==Electron::Definition()) {
        //std::cout<<"  ELECTRON" <<std::endl;
        //
        // create ionization process for e- with 1 model:
       //
        EMPhysicsProcess *eIoniProc = new ElectronIonizationProcess("eIoniFore-");
        // create the Moller-Bhabha model for ionization i.e. for e- + e- -> e- + e- intercation
        EMModel          *eMBModel  = new MollerBhabhaIonizationModel(true);
        // set min/max energies of the model
        eMBModel->SetLowEnergyUsageLimit ( 1.0*geant::keV);
        eMBModel->SetHighEnergyUsageLimit(100.0*geant::TeV);
        // add the model to the process
        eIoniProc->AddModel(eMBModel);
        //
        // add the process to the e- particle
        AddProcessToParticle(particle, eIoniProc);
        //
        // create bremsstrahlung process for e- with 2 models:
        //
        EMPhysicsProcess *eBremProc = new ElectronBremsstrahlungProcess("eBremFore-");
        // create a SeltzerBergerBremsModel for e-
        EMModel          *eSBModel  = new SeltzerBergerBremsModel(true);
        // set min/max energies of the model
        eSBModel->SetLowEnergyUsageLimit (1.0*geant::keV);
        eSBModel->SetHighEnergyUsageLimit(1.0*geant::GeV);
        // how to inactivate this model in a given region i.e. region with index 1
        // active regions for a model are set based on their process active regions + user requested inactive regions
        //eSBModel->AddToUserRequestedInActiveRegions(1);
        //
        // add this model to the process
        eBremProc->AddModel(eSBModel);
        //
        // create a RelativisticBremsModel for e-
        EMModel          *eRelBModel = new RelativisticBremsModel();
        // set min/max energies of the model
        eRelBModel->SetLowEnergyUsageLimit ( 1.0*geant::GeV);
        eRelBModel->SetHighEnergyUsageLimit(100.0*geant::TeV);
        // add this model to the process
        eBremProc->AddModel(eRelBModel);
        //
        // add the process to the e- particle
        AddProcessToParticle(particle, eBremProc);
      }
      if (particle==Positron::Definition()) {
        //std::cout<<"  Positron" <<std::endl;
        //
        // create ionization process for e+ with 1 model:
        //
        EMPhysicsProcess *eIoniProc = new ElectronIonizationProcess("eIoniFore+");
        // create the Moller-Bhabha model for ionization i.e. for e+ + e- -> e+ + e- intercation
        EMModel          *eMBModel  = new MollerBhabhaIonizationModel(false);
        // set min/max energies of the model
        eMBModel->SetLowEnergyUsageLimit ( 1.0*geant::keV);
        eMBModel->SetHighEnergyUsageLimit(10.0*geant::TeV);
        // add the model to the process
        eIoniProc->AddModel(eMBModel);
        // add the process to the e+ particle
        AddProcessToParticle(particle, eIoniProc);
        //
        // create bremsstrahlung process for e+ with 2 models:
        //
        EMPhysicsProcess *eBremProc = new ElectronBremsstrahlungProcess("eBremFore+");
        // create a SeltzerBergerBremsModel for e-
        EMModel          *eSBModel  = new SeltzerBergerBremsModel(false);
        // set min/max energies of the model
        eSBModel->SetLowEnergyUsageLimit (1.0*geant::keV);
        eSBModel->SetHighEnergyUsageLimit(1.0*geant::GeV);
        // how to inactivate this model in a given region i.e. region with index 1
        // active regions for a model are set based on their process active regions + user requested inactive regions
        //eSBModel->AddToUserRequestedInActiveRegions(1);
        //
        // add this model to the process
        eBremProc->AddModel(eSBModel);
        //
        // create a RelativisticBremsModel for e+
        EMModel          *eRelBModel = new RelativisticBremsModel();
        // set min/max energies of the model
        eRelBModel->SetLowEnergyUsageLimit ( 1.0*geant::GeV);
        eRelBModel->SetHighEnergyUsageLimit(10.0*geant::TeV);
        // add this model to the process
        eBremProc->AddModel(eRelBModel);
        //
        // add the process to the e+ particle
        AddProcessToParticle(particle, eBremProc);
      }
// #endif
      if (particle==Gamma::Definition()) {
// #if 0      
        std::cout<<"**** PhysicsList1:  Gamma - defining processes **************** " <<std::endl;

        EMPhysicsProcess *PhotoElectricProc= new PhotoElectricProcess("Photo Electric Process (Wrapped VecPhys)");
        // EMModel          *photoElecticModel;
        // set min/max energies of the model
        // photoElectricModel->SetLowEnergyUsageLimit (   1.0*geant::keV);
        // photoElectricModel->SetHighEnergyUsageLimit( 100.0*geant::MeV);
        std::cout<< " Adding Photo-Electric to gamma." << std::endl;
        AddProcessToParticle(particle, PhotoElectricProc);
        std::cout<< " Adding Photo-Electric to gamma - done." << std::endl;
        //
        EMPhysicsProcess *ComptonProcess= new class ComptonProcess("Compton Process (Wrapped VecPhys)");
        // EMModel          *comptonModel;
        // set min/max energies of the model
        // comptonModel->SetLowEnergyUsageLimit ( 1.0*geant::keV);
        // comptonModel->SetHighEnergyUsageLimit(10.0*geant::GeV);
        std::cout<< " Adding Compton to gamma."<<std::endl;
        AddProcessToParticle(particle, ComptonProcess);
        std::cout<< " Adding Compton to gamma - done."<<std::endl;
// #endif
#if 0        
        EMPhysicsProcess *ConversionProcess= new class ConversionProcess("Conversion Process (Wrapped VecPhys)");
        std::cout<< " Adding Conversion to gamma."<<std::endl;
        AddProcessToParticle(particle, ConversionProcess);
        std::cout<< " Adding Conversion to gamma - done."<<std::endl;        
#endif
      }        
    }
  }
};

}
