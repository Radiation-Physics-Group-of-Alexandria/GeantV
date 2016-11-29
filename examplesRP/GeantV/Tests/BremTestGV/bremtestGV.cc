#include <iostream>
#include <vector>
#include <cmath>
#include <ctime>

// this is NOT the GeantV but our local GeantTaskData
#include "PhysicsData.h"
#include "LightTrack.h"
//#include "GeantTaskData.h"


#include "PhysicsListManager.h"

#include "PhysicsList.h"
#include "PhysicsList1.h"
//#include "PhysicsList2.h"

#include "Material.h"
#include "MaterialCuts.h"
#include "Region.h"


// just to delete and clean everything ragarding the physics at the end of run
//#include "ELossTableManager.h"
//#include "ELossTableRegister.h"


#include "PhysicsManagerPerParticle.h"
#include "Electron.h"


#include "EMPhysicsProcess.h"
#include "EMModelManager.h"
#include "EMModel.h"

#include "SeltzerBergerBremsModel.h"
#include "RelativisticBremsModel.h"

//#include "GenLogSpacedGrid.h"

//#include "Random.h"
#include "PhysicsProcess.h"

#include "Hist.h"

//#include "GeantTaskData.h"
int main() {

using geantphysics::PhysicsListManager;
using geantphysics::PhysicsList;
using geantphysics::PhysicsList1;
//using geant::PhysicsList2;

using geantphysics::Material;
using geantphysics::MaterialCuts;  // this is just to print the table
using geantphysics::Region;

using geantphysics::PhysicsManagerPerParticle;
using geantphysics::Electron;

using geantphysics::EMPhysicsProcess;
using geantphysics::EMModelManager;
using geantphysics::EMModel;

using geantphysics::SeltzerBergerBremsModel;
using geantphysics::RelativisticBremsModel;

//using geant::GenLogSpacedGrid;
//using geant::Random;
using geantphysics::PhysicsProcess;
//using geant::Hist;

using geantphysics::LightTrack;
using geantphysics::PhysicsData;

using geant::Hist;

//
// Create materials and add them to some regions: (MaterialCuts will be created automatically)
//
// create a material
Material *matAl = Material::NISTMaterial("NIST_MAT_Al");
//Material *matAu = Material::NISTMaterial("NIST_MAT_Au");
//Material *matPbWO4 = Material::NISTMaterial("NIST_MAT_PbWO4");

// print the material table
//std::cerr<< Material::GetTheMaterialTable();


// create some dummy regions; add some materials;
/*
bool   iscutinlength = true;
double gcut          = 1.0*geant::mm;
double emcut         = 1.0*geant::mm;
double epcut         = 1.0*geant::mm;
*/
bool   iscutinlength = false;
double gcut          = 2.0*geant::keV;
double emcut         = 2.0*geant::keV;
double epcut         = 2.0*geant::keV;

Region *reg0    = new Region("REGION_0",iscutinlength, gcut, emcut, epcut);
reg0->AddMaterial(matAl);
//reg0->AddMaterial(matAu);
//reg0->AddMaterial(matPbWO4);
////Region *reg1    = new Region("REGION_1",iscutinlength, gcut*1.01, emcut*10, epcut*10);
///reg1->AddMaterial(matAu);

// create all MaterialCuts
MaterialCuts::CreateAll();

std::cout<<MaterialCuts::GetTheMaterialCutsTable()<<std::endl;

// set number of regions in the PhysicsListManager: before we start to use it must be set.
PhysicsListManager::Instance().SetNumberOfRegions(Region::GetTheRegionTable().size());
// it will be done automatically later in the PhysicsListManager

//
//  THIS IS VERY SIMILAR To Geant4 STYLE
//  The same as above but if we have only one physics list and the active region vector is not provided then that one
//  phyics list will be used in all regions
//

    // this is what the user needs to be done
    PhysicsList *physList1 = new PhysicsList1("Physics-List-1");
    PhysicsListManager::Instance().RegisterPhysicsList(physList1);

    // this we will do after the user has done its part
    PhysicsListManager::Instance().BuildPhysicsLists();

//    PhysicsListManager::Instance().PrintAll();







    // get one MaterialCuts from the table
    const MaterialCuts *matCut = MaterialCuts::GetMaterialCut(0);
//    double aaGet  = Electron::Definition()->GetPhysicsManagerPerParticlePerRegion(matCut->GetRegionIndex())->GetInvTotalLambda(matCut,1.0*geant::GeV);

    // the models that we are going for
    SeltzerBergerBremsModel *bremSB   = nullptr; // SeltzerBergerBremsModel
    double eMinSB, eMaxSB; // min/max energy usage limits
    RelativisticBremsModel  *bremRel  = nullptr; // RelativisticBremsModel
    double eMinRel, eMaxRel; // min/max energy usage limits

    const std::vector<PhysicsProcess*> processVector = Electron::Definition()->GetPhysicsManagerPerParticlePerRegion(matCut->GetRegionIndex())->GetListProcesses();
    for (int i=0; i<processVector.size(); ++i) {
      std::cout <<" e- process name = " << processVector[i]->GetName() << std::endl;
      if (processVector[i]->GetName()=="eBremFore-") {
         EMPhysicsProcess *emProcess          = static_cast<EMPhysicsProcess*>(processVector[i]);
         const EMModelManager *emModelMgr     = emProcess->GetModelManager();
         const std::vector<EMModel*> emModels = emModelMgr->GetModelListInRegion(matCut->GetRegionIndex());
         for (int j=0; j<emModels.size(); ++j) {
           std::cout<< " e- models in EMProcess = "<<emProcess->GetName()<<"  "<< j <<" ==> "<<emModels[j]->GetName()<<std::endl;
           std::cout<< "   ==> used between E_min = "<<emModels[j]->GetLowEnergyUsageLimit()/geant::keV
                    << " [keV] and E_max = " << emModels[j]->GetHighEnergyUsageLimit()/geant::GeV << " [GeV]"
                    <<std::endl;
          if (emModels[j]->GetName()=="eSeltzerBergerBrems") {
            bremSB  = static_cast<SeltzerBergerBremsModel*>(emModels[j]);
            eMinSB  = emModels[j]->GetLowEnergyUsageLimit();
            eMaxSB  = emModels[j]->GetHighEnergyUsageLimit();
          } else if (emModels[j]->GetName()=="eRelativisticBrems") {
            bremRel = static_cast<RelativisticBremsModel*>(emModels[j]);
            eMinRel = emModels[j]->GetLowEnergyUsageLimit();
            eMaxRel = emModels[j]->GetHighEnergyUsageLimit();
          }
         }
      }
    }



    double ekin   = 10.0*geant::MeV;  // e- kinetic energy
    double dirx   = 0.0;   // direction
    double diry   = 0.0;
    double dirz   = 1.0;
    int    gvcode = Electron::Definition()->GetInternalCode(); // internal code of e-






    Geant::GeantTaskData *td = new Geant::GeantTaskData(1,1,1);
    PhysicsData *phd = new PhysicsData();
    td->fPhysicsData = phd;
    // set up a the primary light track for brem.
    LightTrack primaryLT;
    //LightTrack secondaryLT;
    std::vector<LightTrack> secLt;  // dummy because we fill secondaries into GeantTaskData::PhysicsData


    double gamProdCut = matCut->GetProductionCutsInEnergy()[0];
    Hist *h = new Hist(std::log10(gamProdCut/ekin), 0.01, 100);
//    Hist *h = new Hist(-6., 0.01, 100);

    long int NUMRAND = 10000000;
    clock_t  start_time = clock();
    for (long int i=0;i<NUMRAND;++i) {
       // we will use members:
       //  fMaterialCutCoupleIndex <==>  // current MaterialCuts index
       //  fKinE         <==>  fE-fMass  // kinetic energy;  will be set to the new kinetic energy
       //  fGVcode       <==>  fGVcode   // internal particle code
       //  fIntLen       <==>  fIntLen   // pre-step lambda for accounting energy loss i.e. to see if it is a delta inter.
       //  fXdir         <==>  fXdir     // direction vector x comp. will be set to the new direction x comp.
       //  fYdir         <==>  fYdir     // direction vector y comp. will be set to the new direction y comp.
       //  fZdir         <==>  fZdir     // direction vector z comp. will be set to the new direction z comp.
       primaryLT.SetMaterialCutCoupleIndex(matCut->GetIndex());
       primaryLT.SetKinE(ekin);
       primaryLT.SetGVcode(gvcode);
//       primaryLT.SetTrackIndex(0); // not important now
       primaryLT.SetDirX(dirx);
       primaryLT.SetDirY(diry);
       primaryLT.SetDirZ(dirz);
//       primaryLT.SetTotalMFP(1.0); // not important now
       //
       // clean the number of secondary tracks used (in PhysicsData)
       td->fPhysicsData->SetNumUsedSecondaries(0);
       //
       // invoke the interaction
       int numSecs = bremSB->SampleSecondaries(primaryLT,secLt,td);
       //int numSecs = bremRel->SampleSecondaries(primaryLT,secLt,td);
        
        if (numSecs>0) {
         LightTrack &secondaryLT = ((td->fPhysicsData->GetListOfSecondaries())[0]);
         double egamma = secondaryLT.GetKinE();
         if (egamma<gamProdCut) {
           std::cerr<<"  ***  gammae = "<<egamma << " < gamProdCut = "<<gamProdCut <<std::endl;
           exit(-1);
         }
         h->Fill(std::log10(egamma/ekin),1.0);
       }
     }
     clock_t end_time = clock();
     std::cerr<<" --- Time = "<<(end_time-start_time)/(double(CLOCKS_PER_SEC))<<std::endl;


     double norm = 0.25/((double)NUMRAND);
     for (int i=0; i< h->GetNumBins(); ++i) {
      std::cout<<h->GetX()[i]+0.5*h->GetDelta()<<"  "<<std::setprecision(8)<<h->GetY()[i]*norm<<std::endl;
     }


return 0;
}
