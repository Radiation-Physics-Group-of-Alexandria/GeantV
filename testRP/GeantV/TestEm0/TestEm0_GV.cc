
/**
 * @brief  TestEm0: simple GeantV-Real-physics application to test integarted quatities.
 * @author M Novak
 * @date   April 2017
 *
 * A simple GeantV-Real-physics application to test integarted physics quantities (e.g.
 * atomic and macroscopic cross sections, mean free path, restricted and unrestricted
 * stopping powers and the corresponding ranges) for a given particle, target material,
 * particle kinetic energy and secondary production threshold. EM processes, that has
 * discrete part and assigned to the given particle in the given physics-list, are
 * extracted from the physics-list and the integrated quantities are computed by calling
 * EMProcess/EMModel interface methods. Quantites, that are available in tables (built
 * at the physics initialisation), also extracted (i.e. interpolated) from these tables
 * by using the corresponding interface methods. Note, that atomic cross sections are
 * computed only for materials that has a single elemnt.
 *
 * Run ./TestEm0_GV --help    for more details!
 * The corresponding Geant4 test application is TestEm0_G4.
 */

#include <iostream>
#include <vector>

#include <getopt.h>
#include <err.h>

#include "Material.h"
#include "Element.h"
#include "MaterialCuts.h"

// vecgeom includes
#include "volumes/LogicalVolume.h"
#include "volumes/Box.h"

#include "Region.h"
#include "PhysicsListManager.h"
#include "PhysicsList.h"
#include "PhysicsParameters.h"

// include the user defined physics list
#include "UserPhysicsList.h"

#include "ELossTableManager.h"


#include "PhysicsManagerPerParticle.h"
#include "PhysicsProcess.h"

#include "Particle.h"
#include "Electron.h"
#include "Positron.h"
#include "Gamma.h"

#include "EMPhysicsProcess.h"
#include "EMModelManager.h"
#include "EMModel.h"


using geantphysics::Material;
using geantphysics::Element;
using geantphysics::MaterialCuts;  // this is just to print the table

using geantphysics::PhysicsListManager;
using geantphysics::PhysicsList;
using geantphysics::PhysicsParameters;


using geantphysics::PhysicsManagerPerParticle;
using geantphysics::PhysicsProcess;

using geantphysics::Particle;
using geantphysics::Electron;
using geantphysics::Positron;
using geantphysics::Gamma;

using geantphysics::EMPhysicsProcess;

using geantphysics::ELossTableManager;

using geantphysics::EMModelManager;
using geantphysics::EMModel;

//
// default values of the input parameters
static std::string   particleName("e-");            // primary particle is electron
static std::string   materialName("NIST_MAT_Pb");   // material is lead
static double        primaryEnergy     = 0.1;       // primary particle energy in [GeV]
static double        prodCutValue      = 0.1;       // by default in length and internal units i.e. [cm]
static bool          isProdCutInLength = true;      // is the production cut value given in length ?

static struct option options[] = {
  {"particle-name    (possible particle names: e-, e+, gamma)                 - default: e-"          , required_argument, 0, 'p'},
  {"material-name    (with a NIST_MAT_ prefix; see more in material doc.)     - default: NIST_MAT_Pb" , required_argument, 0, 'm'},
  {"primary-energy   (in internal energy units i.e. [GeV])                    - default: 0.1"         , required_argument, 0, 'E'},
  {"cut-vale         (secondary production threshold value for all particles) - default: 0.1"         , required_argument, 0, 'c'},
  {"cut-in-energy    (is the production cut value given in energy ? )         - default: NO"          , no_argument, 0, 'e'},
  {"help"                                                                                             , no_argument, 0, 'h'},
  {0, 0, 0, 0}
};
void help();

//===========================================================================================//
int main(int argc, char *argv[]) {
  //
  // Get input parameters
  while (true) {
    int c, optidx = 0;
    c = getopt_long(argc, argv, "ehp:m:E:c:", options, &optidx);
    if (c == -1)
      break;
    switch (c) {
    case 0:
      c = options[optidx].val;
    /* fall through */
    case 'p':
       particleName = optarg;
        break;
    case 'm':
       materialName = optarg;
       break;
    case 'E':
      primaryEnergy = (double)strtof(optarg, NULL);
      if (primaryEnergy<=0)
        errx(1, "primary particle energy must be positive");
      break;
    case 'c':
      prodCutValue = (double)strtof(optarg, NULL);
      if (prodCutValue<=0)
        errx(1, "production cut value must be positive");
      break;
    case 'e':
      isProdCutInLength = false;
      break;
    case 'h':
       help();
       return 0;
       break;
    default:
      help();
      errx(1, "unknown option %c", c);
    }
  }

//============================= User defined input data =====================================//
//
// Create target material: which is supposed to be a NIST Material
Material *matDetector = Material::NISTMaterial(materialName);
//
// Set particle (Electron,Positron,Gamma)
Particle *particle = nullptr;
if (particleName=="e-") {
  particle = Electron::Definition();
} else if (particleName=="e+") {
  particle = Positron::Definition();
} else if (particleName=="gamma") {
  particle = Gamma::Definition();
} else {
  std::cerr<< " ***** ERROR: unknown particle name: " << particleName << std::endl;
  help();
  return 0;
}
//
// Set particle kinetic energy
double kineticEnergy = primaryEnergy;

//
// Set production cuts if needed
bool   iscutinlength   = isProdCutInLength;
double gcut            = prodCutValue;
double emcut           = prodCutValue;
double epcut           = prodCutValue;
//===========================================================================================//



//============= Initialization i.e. building up and init the physics ========================//
//
// Create a dummy vecgeom::geometry with only one volume i.e. world; create a region and set production cuts
//
// create a vecgeom::LogicalVolume
vecgeom::UnplacedBox worldParams = vecgeom::UnplacedBox(1.,1.,1.);
vecgeom::LogicalVolume worldl(&worldParams);
// create one region and assigne to the logical volume
vecgeom::Region *aRegion = new vecgeom::Region("ARegion",iscutinlength, gcut, emcut, epcut);
worldl.SetRegion(aRegion);
// create a dummy vecgeom::media for this volume:
// Only the vecgeom::material index (that will be 0) is used later to identify our material during the init.
worldl.SetTrackingMediumPtr(new vecgeom::Medium((matDetector->GetName()).c_str(),new vecgeom::Material(),new double()));
vecgeom::GeoManager::Instance().SetWorld(worldl.Place());
vecgeom::GeoManager::Instance().CloseGeometry();

// print the material table
//std::cerr<< Material::GetTheMaterialTable();

// We set the min value of the secondary production thresholds to 100eV according to Geant4/TestEm0
// NOTE: these must be set before creating all the MaterialCuts-s
PhysicsParameters::SetMinAllowedGammaCutEnergy   (100.*geant::eV);
PhysicsParameters::SetMinAllowedElectronCutEnergy(100.*geant::eV);
PhysicsParameters::SetMinAllowedPositronCutEnergy(100.*geant::eV);

//
// Create all MaterialCuts
//
MaterialCuts::CreateAll();
// print all MaterialCuts
// std::cout<<MaterialCuts::GetTheMaterialCutsTable()<<std::endl;




//
// Set number of regions in the PhysicsListManager: during normal physics init., it is done
// automatically in the PhysicsProcessHandler::Initialize()
PhysicsListManager::Instance().SetNumberOfRegions(vecgeom::Region::GetNumberOfRegions());


//
// Create one user physics list.
//
//   THIS IS VERY SIMILAR To Geant4 STYLE
//   The same as above but if we have only one physics list and the active region vector is not provided then that one
//   phyics list will be used in all regions

    // this is what the user needs to do:
    // 1. create the physics list
    PhysicsList *thePhysicsList = new userapplication::UserPhysicsList("UserDefined-PhysicsList");
    // The user can get the PhysicsParameters form the physics list and change default values:
    //  - unlike in case of usual simulation, now we request to compute the CSDA range table
    thePhysicsList->GetPhysicsParameters()->SetIsComputeCSDARange(true);
    // 2. register the physics list:
    PhysicsListManager::Instance().RegisterPhysicsList(thePhysicsList);

//
// print out the PhysicsParameters obejct: we ha only one physics list so we have only one physics parameter object
 std::cout<<PhysicsParameters::GetThePhysicsParametersTable()[0];

// Build all registered physics lists: during normal physics init., it is done
// automatically in the PhysicsProcessHandler::Initialize()
PhysicsListManager::Instance().BuildPhysicsLists();
//===========================================================================================//



//======================== Getting the required information =================================//
//
//  Now get the list of EM processes assigned to the particle and use them.
//
std::cout<<std::endl<<std::endl;
std::cout<< "   ================================================================================ \n"
         << "   ================================    RESULTS    ================================= \n"
         << "   ================================================================================ \n";

// first get the MaterialCuts: we have only one
const MaterialCuts *matCut = MaterialCuts::GetMaterialCut(aRegion->GetIndex(),matDetector->GetIndex());
std::cout<< "  "<< matCut->GetMaterial() << std::endl;
std::cout<< "   -------------------------------------------------------------------------------- "<<std::endl;
std::cout<< "   MaterialCuts: \n" << matCut;
std::cout<< "   -------------------------------------------------------------------------------- "<<std::endl;
std::cout<< "   Particle       =  " << particle->GetName() << std::endl;
std::cout<< "   -------------------------------------------------------------------------------- "<<std::endl;
std::cout<< "   Kinetic energy =  " << kineticEnergy/geant::MeV << "  [MeV] " << std::endl;
std::cout<< "   -------------------------------------------------------------------------------- "<<std::endl;
std::cout<< "   Physics list   =  " << thePhysicsList->GetName() << std::endl;
std::cout<< "   -------------------------------------------------------------------------------- "<<std::endl;


//
// Get all processes, assigned to this partcile and active in the given region, that has discrete part.
PhysicsManagerPerParticle *thePhysManager = particle->GetPhysicsManagerPerParticlePerRegion(matCut->GetRegionIndex());
// check if the particle has any processes assined to (in the current region): if not than the thePhysManager = nullptr
size_t numPostStepCandidateProcesses = 0;
std::vector<PhysicsProcess*> thePostStepCandProcVect;
if (thePhysManager) {
  thePostStepCandProcVect       = thePhysManager->GetListPostStepCandidateProcesses();
  numPostStepCandidateProcesses = thePostStepCandProcVect.size();
}
std::cout<< "   The particle has " << numPostStepCandidateProcesses
         << " processes with discrete part assinged to." << std::endl;
std::cout<< "   -------------------------------------------------------------------------------- "<<std::endl;
// terminate if there is no any post-step candidate processes assigned to the particle in the given region.
if (!numPostStepCandidateProcesses) {
  std::cout << "   ================================================================================ "<< std::endl << std::endl;
  return 0;
}


double  compTotalMacXsec = 0.0;  // computed total macroscopic cross section i.e. summed up per-process mac. x-secs
double  getTotalMacXsec  = 0.0;  // same but interpolated from the lambda table that is built at initialization
double  compTotalRestrictedDEDX   = 0.0; // the corresponding can be obtained from table as well
double  compTotalUnRestrictedDEDX = 0.0; // the corresponding can be obtained from table as well (if CSDA range is TRUE)
double  compTotalAtomicXsec       = 0.0; // computed total atomic cross section (only if material has 1 element)
std::vector<std::string> processNameVect; // process names
std::vector<double>  compMacXsecPerProcessVect; // computed macroscopic cross sections per-process
std::vector<double>  getMacXsecPerProcessVect;  // same but interpolated from the lambda table that is built at init.
std::vector<double>  compRestrictedDEDXVect;    // computed restricted dE/dx per process (only if there is 1 ElossP.)
std::vector<double>  compUnRestrictedDEDXVect;  // computed full dE/dx per process (only if there is 1 ElossP.)
std::vector<double>  compAtomicXsectionVect;    // computed atomic cross sections (only if material has 1 element)


// quick check if the particle has kEnergyLoss process (that is cont-discrete so...)
bool isHasEnergyLossProcess = false;
for (size_t i=0; i<thePostStepCandProcVect.size(); ++i) {
  if (thePostStepCandProcVect[i]->GetType()==geantphysics::ProcessType::kEnergyLoss) {
    isHasEnergyLossProcess = true;
    break;
  }
}

bool isSingleElementMaterial = false;
if (matCut->GetMaterial()->GetNumberOfElements()==1) {
  isSingleElementMaterial = true;
}

getTotalMacXsec = thePhysManager->GetInvTotalLambda(matCut, kineticEnergy);
for (size_t i=0; i<thePostStepCandProcVect.size(); ++i) {
  PhysicsProcess *proc = thePostStepCandProcVect[i];
  processNameVect.push_back(proc->GetName());
  compMacXsecPerProcessVect.push_back(proc->ComputeMacroscopicXSection(matCut, kineticEnergy, particle));
  compTotalMacXsec += compMacXsecPerProcessVect[i];
  getMacXsecPerProcessVect.push_back(thePhysManager->GetMacroscopicXSectionForProcess(matCut, kineticEnergy, i));
  EMPhysicsProcess *emProc = nullptr;
  if (proc->GetType()==geantphysics::ProcessType::kElectromagnetic
      || proc->GetType()==geantphysics::ProcessType::kEnergyLoss) {
    emProc = static_cast<EMPhysicsProcess*>(proc);
  }
  if (isHasEnergyLossProcess) {
    if (emProc) {
      compRestrictedDEDXVect.push_back(emProc->ComputeDEDX(matCut, kineticEnergy, particle));
      compUnRestrictedDEDXVect.push_back(emProc->ComputeDEDX(matCut, kineticEnergy, particle, true));
    } else {
      compRestrictedDEDXVect.push_back(0.0);
      compUnRestrictedDEDXVect.push_back(0.0);
    }
    compTotalRestrictedDEDX   += compRestrictedDEDXVect[i];
    compTotalUnRestrictedDEDX += compUnRestrictedDEDXVect[i];
  }
  if (isSingleElementMaterial) {
    if (emProc) {
      EMPhysicsProcess     *emProc         = static_cast<EMPhysicsProcess*>(proc);
      EMModelManager       *emModelManager = emProc->GetModelManager();
      EMModel              *emModel        = emModelManager->SelectModel(kineticEnergy, matCut->GetRegionIndex());
      const Element        *elem           = (matCut->GetMaterial()->GetElementVector())[0];
      compAtomicXsectionVect.push_back(emModel->ComputeXSectionPerAtom(elem, matCut, kineticEnergy, particle));
    } else {
      compAtomicXsectionVect.push_back(0.0);
    }
    compTotalAtomicXsec += compAtomicXsectionVect[i];
  }
}

processNameVect.push_back("total");
compMacXsecPerProcessVect.push_back(compTotalMacXsec);
getMacXsecPerProcessVect.push_back(getTotalMacXsec);
compRestrictedDEDXVect.push_back(compTotalRestrictedDEDX);
compUnRestrictedDEDXVect.push_back(compTotalUnRestrictedDEDX);
compAtomicXsectionVect.push_back(compTotalAtomicXsec);


std::cout<< "   process name                :";
for (size_t i=0; i<processNameVect.size();++i) {
  std::cout<< std::setw(16) << std::right << processNameVect[i] << std::setw(10) << " ";
}
std::cout<<std::endl; std::cout<<std::endl;

if (isSingleElementMaterial) {
  std::cout<< "   cross section per atom      :";
  for (size_t i=0; i<processNameVect.size();++i) {
    std::cout<< std::setw(14) << std::scientific << std::right << compAtomicXsectionVect[i]/(geant::barn)
             << std::setw(14) << std::left << "     [barn]";
  }
  std::cout<<std::endl; std::cout<<std::endl;
}

std::cout<< "   compCrossSectionPerVolume   :";
for (size_t i=0; i<processNameVect.size();++i) {
  std::cout<< std::setw(14) << std::scientific << std::right << compMacXsecPerProcessVect[i]/(1./geant::cm)
           << std::setw(14) << std::left << "     [1/cm]";
}
std::cout<<std::endl;

std::cout<< "   cross section per volume    :";
for (size_t i=0; i<processNameVect.size();++i) {
  std::cout<< std::setw(14) << std::scientific << std::right << getMacXsecPerProcessVect[i]/(1./geant::cm)
           << std::setw(14) << std::left << "     [1/cm]";
}
std::cout<<std::endl;

double density = matCut->GetMaterial()->GetDensity()/(geant::g/geant::cm3); // density in [g/cm3]
std::cout<< "   cross section per mass      :";
for (size_t i=0; i<processNameVect.size();++i) {
  std::cout<< std::setw(14) << std::scientific << std::right << getMacXsecPerProcessVect[i]/density/(1./geant::cm)
           << std::setw(14) << std::left << "    [cm2/g]";
}
std::cout<<std::endl; std::cout<<std::endl;


std::cout<< "   mean free path (length)     :";
for (size_t i=0; i<processNameVect.size();++i) {
  double lambda = PhysicsProcess::GetAVeryLargeValue();
  if (getMacXsecPerProcessVect[i]>0.) {
    lambda = 1./getMacXsecPerProcessVect[i]/geant::cm;
  }
  std::cout<< std::setw(14) << std::scientific << std::right << lambda
           << std::setw(14) << std::left << "       [cm]";
}
std::cout<<std::endl;

std::cout<< "   mean free path (g/cm2)      :";
for (size_t i=0; i<processNameVect.size();++i) {
  double lambda = PhysicsProcess::GetAVeryLargeValue();
  if (getMacXsecPerProcessVect[i]>0.) {
    lambda = getMacXsecPerProcessVect[i]/density/(1./geant::cm);
    lambda = 1./lambda;
  }
  std::cout<< std::setw(14) << std::scientific << std::right << lambda
           << std::setw(14) << std::left << "    [g/cm2]";
}
std::cout<<std::endl; std::cout<<std::endl;


if (isHasEnergyLossProcess) {
  std::cout<< "   restricted dE/dx  (MeV/cm)  :";
  for (size_t i=0; i<processNameVect.size();++i) {
    std::cout<< std::setw(14) << std::scientific << std::right << compRestrictedDEDXVect[i]/(geant::MeV/geant::cm)
             << std::setw(14) << std::left << "   [MeV/cm]";
  }
  std::cout<<std::endl;

  double getRestrictedDEDX = ELossTableManager::Instance().GetRestrictedDEDX(matCut, particle,kineticEnergy)/(geant::MeV/geant::cm);
  std::cout<< "   restricted dE/dx (MeV/cm)   :";
  std::cout<< std::setw(29) << std::setfill('.') << std::right << " ( interpolated from"
           << std::left << " the ELossTable )"
           << std::setw(24) << std::scientific << std::right << getRestrictedDEDX << std::setfill(' ')
           << std::setw(14) << std::left << "   [MeV/cm]";
  std::cout<<std::endl; std::cout<<std::endl;


  std::cout<< "   restricted dE/dx (MeVcm2/g) :";
  for (size_t i=0; i<processNameVect.size();++i) {
    std::cout<< std::setw(14) << std::scientific << std::right << compRestrictedDEDXVect[i]/(geant::MeV/geant::cm)/density
             << std::setw(14) << std::left << " [MeVcm2/g]";
  }
  std::cout<<std::endl;

  std::cout<< "   restricted dE/dx (MeVcm2/g) :";
  std::cout<< std::setw(29) << std::setfill('.') <<std::right << " ( interpolated from"
           << std::left << " the ELossTable )"
           << std::setw(24) << std::scientific << std::right << getRestrictedDEDX/density << std::setfill(' ')
           << std::setw(14) << std::left << " [MeVcm2/g]";
  std::cout<<std::endl; std::cout<<std::endl;


  std::cout<< "   unrestricted dE/dx (MeV/cm) :";
  for (size_t i=0; i<processNameVect.size();++i) {
    std::cout<< std::setw(14) << std::scientific << std::right << compUnRestrictedDEDXVect[i]/(geant::MeV/geant::cm)
             << std::setw(14) << std::left << "   [MeV/cm]";
  }
  std::cout<<std::endl;

  std::cout<< "   unrestricted dE/dx(MeVcm2/g):";
  for (size_t i=0; i<processNameVect.size();++i) {
    std::cout<< std::setw(14) << std::scientific << std::right << compUnRestrictedDEDXVect[i]/(geant::MeV/geant::cm)/density
             << std::setw(14) << std::left << " [MeVcm2/g]";
  }
  std::cout<<std::endl; std::cout<<std::endl;

  double getRestrictedRange = ELossTableManager::Instance().GetRestrictedRange(matCut, particle,kineticEnergy)/(geant::mm);
  std::cout<< "   range from restricted dE/dx :";
  std::cout<< std::setw(29) << std::setfill('.') << std::right << " ( interpolated from"
           << std::left << " the ELossTable )"
           << std::setw(24) << std::scientific << std::right << getRestrictedRange << std::setfill(' ')
           << std::setw(14) << std::left << "       [mm]";
  std::cout<<std::endl;

  double getUnRestrictedRange = ELossTableManager::Instance().GetRange(matCut, particle,kineticEnergy)/(geant::mm);
  std::cout<< "   range from full dE/dx       :";
  std::cout<< std::setw(29) << std::setfill('.') << std::right << " ( interpolated from"
           << std::left << " the ELossTable )"
           << std::setw(24) << std::scientific << std::right << getUnRestrictedRange << std::setfill(' ')
           << std::setw(14) << std::left << "       [mm]";
  std::cout<<std::endl; std::cout<<std::endl;

}

std::cout << "   ================================================================================ "<< std::endl << std::endl;

return 0;
}

void help() {
  std::cout<<"\n "<<std::setw(120)<<std::setfill('=')<<""<<std::setfill(' ')<<std::endl;
  std::cout<<"  TestEm0 like GeantV application for testing integrated physics quantities using a given user physics-list"
           << std::endl;
  std::cout<<"\n  Usage: TestEm0_GV [OPTIONS] \n"<<std::endl;
  for (int i = 0; options[i].name != NULL; i++) {
    printf("\t-%c  --%s\n", options[i].val, options[i].name);
  }
  std::cout<<"\n "<<std::setw(120)<<std::setfill('=')<<""<<std::setfill(' ')<<std::endl;
}
