//
// ********************************************************************
// * License and Disclaimer                                           *
// *                                                                  *
// * The  Geant4 software  is  copyright of the Copyright Holders  of *
// * the Geant4 Collaboration.  It is provided  under  the terms  and *
// * conditions of the Geant4 Software License,  included in the file *
// * LICENSE and available at  http://cern.ch/geant4/license .  These *
// * include a list of copyright holders.                             *
// *                                                                  *
// * Neither the authors of this software system, nor their employing *
// * institutes,nor the agencies providing financial support for this *
// * work  make  any representation or  warranty, express or implied, *
// * regarding  this  software system or assume any liability for its *
// * use.  Please see the license in the file  LICENSE  and URL above *
// * for the full disclaimer and the limitation of liability.         *
// *                                                                  *
// * This  code  implementation is the result of  the  scientific and *
// * technical work of the GEANT4 collaboration.                      *
// * By using,  copying,  modifying or  distributing the software (or *
// * any work based  on the software)  you  agree  to acknowledge its *
// * use  in  resulting  scientific  publications,  and indicate your *
// * acceptance of all terms of the Geant4 Software license.          *
// ********************************************************************
//
//
// $Id$
//
//

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

#include "PhysicsList.hh"

#include "G4ProcessManager.hh"

#include "G4BosonConstructor.hh"
#include "G4LeptonConstructor.hh"
#include "G4MesonConstructor.hh"
#include "G4BosonConstructor.hh"
#include "G4BaryonConstructor.hh"
#include "G4IonConstructor.hh"
#include "G4SystemOfUnits.hh"

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

PhysicsList::PhysicsList():  G4VUserPhysicsList()
			  , theNeutrons(0)
			  , theBertiniNeutron(0)
			  , theFTFPNeutron(0)
                          // , theLEPNeutron(0)
			  , thePiK(0)
			  , theBertiniPiK(0)
			  , theFTFPPiK(0)
			  , thePro(0)
			  , theBertiniPro(0)
			  , theFTFPPro(0)
			  , QuasiElastic(false)
{
  defaultCutValue = 1.0*mm;
  SetVerboseLevel(1);
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

PhysicsList::~PhysicsList()
{
  delete theNeutrons;
  delete theBertiniNeutron;
  delete theFTFPNeutron;
  // delete theLEPNeutron;

  delete thePiK;
  delete theBertiniPiK;
  delete theFTFPPiK;

  delete thePro;
  delete theBertiniPro;
  delete theFTFPPro;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void PhysicsList::ConstructParticle()
{
  // In this method, static member functions should be called
  // for all particles which you want to use.
  // This ensures that objects of these particle types will be
  // created in the program.

  G4BosonConstructor  pBosonConstructor;
  pBosonConstructor.ConstructParticle();

  G4LeptonConstructor pLeptonConstructor;
  pLeptonConstructor.ConstructParticle();

  G4MesonConstructor pMesonConstructor;
  pMesonConstructor.ConstructParticle();

  G4BaryonConstructor pBaryonConstructor;
  pBaryonConstructor.ConstructParticle();

  G4IonConstructor pIonConstructor;
  pIonConstructor.ConstructParticle();

}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void PhysicsList::ConstructProcess()
{
  AddTransportation();
  ConstructEM();
  ConstructDecay();

  //  HadronPhysicsFTFP_BERT_WP();
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

#include "G4PhysicsListHelper.hh"

#include "G4ComptonScattering.hh"
#include "G4GammaConversion.hh"
#include "G4PhotoElectricEffect.hh"

#include "G4eMultipleScattering.hh"
#include "G4eIonisation.hh"
#include "G4eBremsstrahlung.hh"
#include "G4eplusAnnihilation.hh"

#include "G4MuMultipleScattering.hh"
#include "G4MuIonisation.hh"
#include "G4MuBremsstrahlung.hh"
#include "G4MuPairProduction.hh"

#include "G4hMultipleScattering.hh"
#include "G4hIonisation.hh"
#include "G4hBremsstrahlung.hh"
#include "G4hPairProduction.hh"

#include "G4ionIonisation.hh"

#include "TabulatedProcess.hh"
#include "VectorizedProcess.hh"
#include "TabulatedDataManager.hh"
#include "TPartIndex.h"

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void PhysicsList::ConstructEM()
{
  G4PhysicsListHelper* ph = G4PhysicsListHelper::GetPhysicsListHelper();

	auto theParticleIterator = GetParticleIterator();
  theParticleIterator->reset();
  while( (*theParticleIterator)() ){
    G4ParticleDefinition* particle = theParticleIterator->value();
    G4String particleName = particle->GetParticleName();

    if (particleName == "gamma") {
      // gamma
      ph->RegisterProcess(new G4PhotoElectricEffect, particle);
      ph->RegisterProcess(new G4ComptonScattering,   particle);
      ph->RegisterProcess(new G4GammaConversion,     particle);

    } else if (particleName == "e-") {
      //electron

      //choices for different physics list for Ionisation and Bremsstrahlung
      char* plname = getenv("PHYSLIST");

      if ( plname && strcmp(plname,"TabulatedPhysics")==0) {
	G4eIonisation* eIoniProc = new G4eIonisation();
	TabulatedProcess* eIoniWrapperProc =
	  new TabulatedProcess(eIoniProc->GetProcessName(),
			       eIoniProc->GetProcessType(),
			       kIonisation);
	eIoniWrapperProc->SetProcessSubType(eIoniProc->GetProcessSubType());
	eIoniWrapperProc->RegisterProcess(eIoniProc);
	ph->RegisterProcess(eIoniWrapperProc, particle);

	G4eBremsstrahlung* eBremProc = new G4eBremsstrahlung();
	TabulatedProcess* eBremWrapperProc =
	  new TabulatedProcess(eBremProc->GetProcessName(),
			       eBremProc->GetProcessType(),
			       kBrehms);
	eBremWrapperProc->SetProcessSubType(eBremProc->GetProcessSubType());
	eBremWrapperProc->RegisterProcess(eBremProc);
	ph->RegisterProcess(eBremWrapperProc, particle);

      }
      else if ( plname && strcmp(plname,"VectorizedPhysics")==0) {
	//use vectorized physics
	G4eIonisation* eIoniProc = new G4eIonisation();
	VectorizedProcess* eIoniWrapperProc =
	  new VectorizedProcess(eIoniProc->GetProcessName(),
				eIoniProc->GetProcessType());
	eIoniWrapperProc->SetProcessSubType(eIoniProc->GetProcessSubType());
	eIoniWrapperProc->RegisterProcess(eIoniProc);
	ph->RegisterProcess(eIoniWrapperProc, particle);

	G4eBremsstrahlung* eBremProc = new G4eBremsstrahlung();
	VectorizedProcess* eBremWrapperProc =
	  new VectorizedProcess(eBremProc->GetProcessName(),
				eBremProc->GetProcessType());
	eBremWrapperProc->SetProcessSubType(eBremProc->GetProcessSubType());
	eBremWrapperProc->RegisterProcess(eBremProc);
	ph->RegisterProcess(eBremWrapperProc, particle);
      }
      else {
	//user standard electron processes
	ph->RegisterProcess(new G4eMultipleScattering, particle);
	ph->RegisterProcess(new G4eIonisation,         particle);
	ph->RegisterProcess(new G4eBremsstrahlung,     particle);
      }
    } else if (particleName == "e+") {
      //positron
      ph->RegisterProcess(new G4eMultipleScattering, particle);
      ph->RegisterProcess(new G4eIonisation,         particle);
      ph->RegisterProcess(new G4eBremsstrahlung,     particle);
      ph->RegisterProcess(new G4eplusAnnihilation,   particle);

    } else if( particleName == "mu+" ||
               particleName == "mu-"    ) {
      //muon
      ph->RegisterProcess(new G4MuMultipleScattering, particle);
      ph->RegisterProcess(new G4MuIonisation,         particle);
      ph->RegisterProcess(new G4MuBremsstrahlung,     particle);
      ph->RegisterProcess(new G4MuPairProduction,     particle);

    } else if( particleName == "proton" ||
               particleName == "pi-" ||
               particleName == "pi+"    ) {
      //proton
      ph->RegisterProcess(new G4hMultipleScattering, particle);
      ph->RegisterProcess(new G4hIonisation,         particle);
      ph->RegisterProcess(new G4hBremsstrahlung,     particle);
      ph->RegisterProcess(new G4hPairProduction,     particle);

    } else if( particleName == "alpha" ||
               particleName == "He3" )     {
      //alpha
      ph->RegisterProcess(new G4hMultipleScattering, particle);
      ph->RegisterProcess(new G4ionIonisation,       particle);

    } else if( particleName == "GenericIon" ) {
      //Ions
      ph->RegisterProcess(new G4hMultipleScattering, particle);
      ph->RegisterProcess(new G4ionIonisation,       particle);

      } else if ((!particle->IsShortLived()) &&
               (particle->GetPDGCharge() != 0.0) &&
               (particle->GetParticleName() != "chargedgeantino")) {
      //all others charged particles except geantino
      ph->RegisterProcess(new G4hMultipleScattering, particle);
      ph->RegisterProcess(new G4hIonisation,         particle);
    }
  }
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

#include "G4Decay.hh"

void PhysicsList::ConstructDecay()
{
  G4PhysicsListHelper* ph = G4PhysicsListHelper::GetPhysicsListHelper();

  // Add Decay Process
  G4Decay* theDecayProcess = new G4Decay();

	auto theParticleIterator = GetParticleIterator();
  theParticleIterator->reset();
  while( (*theParticleIterator)() ){
    G4ParticleDefinition* particle = theParticleIterator->value();
    if (theDecayProcess->IsApplicable(*particle)) {
      ph->RegisterProcess(theDecayProcess, particle);
    }
  }
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void PhysicsList::SetCuts()
{
  if (verboseLevel >0){
    G4cout << "PhysicsList::SetCuts:";
    G4cout << "CutLength : " << G4BestUnit(defaultCutValue,"Length") << G4endl;
  }

  // set cut values for gamma at first and for e- second and next for e+,
  // because some processes for e+/e- need cut values for gamma
  //
  SetCutValue(defaultCutValue, "gamma");
  SetCutValue(defaultCutValue, "e-");
  SetCutValue(defaultCutValue, "e+");
  SetCutValue(defaultCutValue, "proton");

  if (verboseLevel>0) DumpCutValuesTable();
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void PhysicsList::HadronPhysicsFTFP_BERT_WP()
{
  theNeutrons=new G4NeutronBuilder_WP;
  theFTFPNeutron=new G4FTFPNeutronBuilder(QuasiElastic);
  theNeutrons->RegisterMe(theFTFPNeutron);
  theNeutrons->RegisterMe(theBertiniNeutron=new G4BertiniNeutronBuilder);
  theBertiniNeutron->SetMinEnergy(0.0*GeV);
  theBertiniNeutron->SetMaxEnergy(5*GeV);
  // theNeutrons->RegisterMe(theLEPNeutron=new G4LEPNeutronBuilder);
  // theLEPNeutron->SetMinInelasticEnergy(0.0*eV);   // no inelastic from LEP
  // theLEPNeutron->SetMaxInelasticEnergy(0.0*eV);

  thePro=new G4ProtonBuilder_WP;
  theFTFPPro=new G4FTFPProtonBuilder(QuasiElastic);
  thePro->RegisterMe(theFTFPPro);
  thePro->RegisterMe(theBertiniPro=new G4BertiniProtonBuilder);
  theBertiniPro->SetMaxEnergy(5*GeV);

  thePiK=new G4PiKBuilder_WP;
  theFTFPPiK=new G4FTFPPiKBuilder(QuasiElastic);
  thePiK->RegisterMe(theFTFPPiK);
  thePiK->RegisterMe(theBertiniPiK=new G4BertiniPiKBuilder);
  theBertiniPiK->SetMaxEnergy(5*GeV);

  theNeutrons->Build();
  thePro->Build();
  thePiK->Build();
}
