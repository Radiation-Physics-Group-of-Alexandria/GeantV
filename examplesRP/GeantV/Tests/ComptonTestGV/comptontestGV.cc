#include <iostream>
#include <vector>
#include <cmath>
#include <ctime>

// this is NOT the GeantV but our local GeantTaskData
#include "PhysicsData.h"
#include "LightTrack.h"

#include "PhysicsListManager.h"

#include "PhysicsList.h"
#include "PhysicsList2.h"

#include "Material.h"
#include "MaterialCuts.h"
#include "Region.h"

#include "PhysicsManagerPerParticle.h"
#include "Electron.h"


#include "EMPhysicsProcess.h"
#include "EMModelManager.h"
#include "EMModel.h"

#include "SeltzerBergerBremsModel.h"
#include "RelativisticBremsModel.h"
#include "KleinNishinaComptonModel.h"

#include "PhysicsProcess.h"
#include "Hist.h"

#include "TH1F.h"
#include "TFile.h"


int main() {
    
    using geantphysics::PhysicsListManager;
    using geantphysics::PhysicsList;
    using geantphysics::PhysicsList2;
    
    using geantphysics::Material;
    using geantphysics::MaterialCuts;  // this is just to print the table
    using geantphysics::Region;
    
    using geantphysics::PhysicsManagerPerParticle;
    using geantphysics::Electron;
    using geantphysics::Gamma;
    
    using geantphysics::EMPhysicsProcess;
    using geantphysics::EMModelManager;
    using geantphysics::EMModel;
    
    using geantphysics::SeltzerBergerBremsModel;
    using geantphysics::RelativisticBremsModel;
    using geantphysics::KleinNishinaComptonModel;

    using geantphysics::PhysicsProcess;

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
    double gcut          = 100.0 * geant::eV;//2.0*geant::keV;
    double emcut         = 100.0 * geant::eV;//2.0*geant::keV;
    double epcut         = 100.0 * geant::eV;//2.0*geant::keV;
    
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
    PhysicsList *physList2 = new PhysicsList2("Physics-List-2");
    PhysicsListManager::Instance().RegisterPhysicsList(physList2);
    
    // this we will do after the user has done its part
    PhysicsListManager::Instance().BuildPhysicsLists();
    std::cout<<"Physics list created\n";
    
    //    PhysicsListManager::Instance().PrintAll();
    
    
    // get one MaterialCuts from the table
    const MaterialCuts *matCut = MaterialCuts::GetMaterialCut(0);
    //    double aaGet  = Electron::Definition()->GetPhysicsManagerPerParticlePerRegion(matCut->GetRegionIndex())->GetInvTotalLambda(matCut,1.0*geant::GeV);
    
    // the models that we are going for
    KleinNishinaComptonModel *comptonKN   = nullptr; // KleinNishinaComptonModel
    double eMinKN, eMaxKN; // min/max energy usage limits

    
    const std::vector<PhysicsProcess*> processVector = Gamma::Definition()->GetPhysicsManagerPerParticlePerRegion(matCut->GetRegionIndex())->GetListProcesses();
    
    for (int i=0; i<processVector.size(); ++i) {
        std::cout <<" Process name = " << processVector[i]->GetName() << std::endl;
        if (processVector[i]->GetName()=="gCompton") {
            EMPhysicsProcess *emProcess          = static_cast<EMPhysicsProcess*>(processVector[i]);
            const EMModelManager *emModelMgr     = emProcess->GetModelManager();
            const std::vector<EMModel*> emModels = emModelMgr->GetModelListInRegion(matCut->GetRegionIndex());
            for (int j=0; j<emModels.size(); ++j) {
                std::cout<< " gamma models in EMProcess = "<<emProcess->GetName()<<"  "<< j <<" ==> "<<emModels[j]->GetName()<<std::endl;
                std::cout<< "   ==> used between E_min = "<<emModels[j]->GetLowEnergyUsageLimit()/geant::keV
                << " [keV] and E_max = " << emModels[j]->GetHighEnergyUsageLimit()/geant::GeV << " [GeV]"
                <<std::endl;
                if (emModels[j]->GetName()=="KleinNishinaCompton") {
                    
                    comptonKN = static_cast<KleinNishinaComptonModel*>(emModels[j]);
                    eMinKN = emModels[j]->GetLowEnergyUsageLimit();
                    eMaxKN = emModels[j]->GetHighEnergyUsageLimit();
                    std::cout<<"**********\n\neminKN: "<<eMinKN<<"\nemaxKN: "<<eMaxKN<<std::endl;
                    std::cout<<comptonKN->GetName()<<std::endl;
                }
            }
        }
    }
    
    double ekin   = 10.0*geant::GeV;  // gamma kinetic energy
    double dirx   = 0.0;              // direction
    double diry   = 0.0;
    double dirz   = 1.0;
    int    gvcode = Gamma::Definition()->GetInternalCode(); // internal code of gamma
    std::cout<<"***********CODE: "<<gvcode<<std::endl;
    
    Geant::GeantTaskData *td = new Geant::GeantTaskData(1,1,1);
    PhysicsData *phd = new PhysicsData();
    td->fPhysicsData = phd;
    // set up a the primary light track for brem.
    LightTrack primaryLT;
    //LightTrack secondaryLT;
    std::vector<LightTrack> secLt;  // dummy because we fill secondaries into GeantTaskData::PhysicsData
    
    double eProdCut = matCut->GetProductionCutsInEnergy()[0];
    Hist *h = new Hist((eProdCut/ekin), 1, 100);
    Hist *hcosTheta = new Hist(-1, 1, 100);
    
    
    TFile* fHistFile = new TFile("gV_new.root", "RECREATE");
    TH1F* fEnergyOut1 = new TH1F("EnergyOut1", "EnergyOut1", 100,  (eProdCut/ekin), 1);
    TH1F* fEnergyOut2 = new TH1F("EnergyOut2", "EnergyOut2", 100,  (eProdCut/ekin), 1);
    TH1F* fAngleOut1 = new TH1F("AngleOut1", "AngleOut1", 100, -1., 1.0);
    TH1F* fAngleOut2 = new TH1F("AngleOut2", "AngleOut2", 100, -1., 1.0);
    
    
    
    
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
        int numSecs = comptonKN->SampleSecondaries(primaryLT,secLt,td);
        
        
        
        if (numSecs>0) {
            LightTrack &secondaryLT = ((td->fPhysicsData->GetListOfSecondaries())[0]);
            double eelectron = secondaryLT.GetKinE();
            double egamma = primaryLT.GetKinE();

            double gcostheta = primaryLT.GetDirZ();
            double ecostheta = secondaryLT.GetDirZ();
            
            if (eelectron<eProdCut) {
                std::cerr<<"  ***  electronenergy = "<<eelectron << " < eProdCut = "<<eProdCut <<std::endl;
                exit(-1);
            }
            //h->Fill(std::log10(eelectron/ekin),1.0);
            
            if(std::abs(egamma)>2*1.0000e-12)
            h->Fill((egamma/ekin),1.0);
            //h->Fill(std::log10(eelectron/ekin),1.0);

            hcosTheta->Fill(gcostheta,1.0);
            //hcosTheta->Fill(gcostheta,1.0);
            
            if(std::abs(egamma)>2*1.0000e-12)
            fEnergyOut1->Fill(egamma/ekin);
            fEnergyOut2->Fill(eelectron/ekin);
            fAngleOut1->Fill(gcostheta);
            fAngleOut2->Fill(ecostheta);

        }
    }
    clock_t end_time = clock();
    std::cerr<<" --- Time = "<<(end_time-start_time)/(double(CLOCKS_PER_SEC))<<std::endl;
    
    
    double norm = 1./((double)NUMRAND);
    for (int i=0; i< h->GetNumBins(); ++i) {
        std::cout<<h->GetX()[i]+0.5*h->GetDelta()<<"  "<<std::setprecision(8)<<h->GetY()[i]*norm<<std::endl;
    }
    
    FILE * fp = fopen("gv.ascii", "w");
    for (int index = 0; index< h->GetNumBins(); index++)
        fprintf(fp, "%.8f %.8f\n", h->GetX()[index]+0.5*h->GetDelta(), h->GetY()[index]*norm);
        fclose(fp);
    

    
    //for (int i=0; i< hcosTheta->GetNumBins(); ++i) {
    //    std::cout<<hcosTheta->GetX()[i]+0.5*h->GetDelta()<<"  "<<std::setprecision(8)<<hcosTheta->GetY()[i]*norm<<std::endl;
    //}
    
    
    fHistFile->Write();
    fHistFile->Close();
    
    //delete fEnergyOut1;
    //delete fEnergyOut2;
    //delete fAngleOut1;
    //delete fAngleOut2;
    return 0;
}
