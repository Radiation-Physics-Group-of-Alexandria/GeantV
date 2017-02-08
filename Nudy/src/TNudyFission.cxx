//SPB
// From this class, we are getting the neutron fission cross section at temprature of 293.606087K.
//This class provides the energy and angle of prompt neutrons produced in fission process.
//Also provides fission fragment mass, charge, no. of prompt and delayed neutron, fission heat etc.
#include "TNudyCore.h"
#include "TNudyEndfRecoPoint.h"
#include "TNudyFission.h"
#include "TNudyEndfMat.h"
#include "TNudyElement.h"
#include <fstream>
#include <iostream>

#ifdef USE_ROOT
#include "TRandom3.h"
#endif
#include "TCanvas.h"
#include "TFile.h"
using namespace std;
#ifdef USE_ROOT
ClassImp(TNudyFission)
#endif

TNudyFission::TNudyFission(){}
TNudyFission::TNudyFission(int ElemId, const char* rootfileName)
{
    ielemId =  ElemId;// for element number like say for H=1, U=2, the way 
                      //we number the elements for material having composite elements
    fRnd     = new TRandom3(0);
    const char* irENDF=rootfileName;//"/home/shiba/geant/rootendf/O16.root";
    rp = new TNudyEndfRecoPoint(ielemId, irENDF);
    mtValues = rp->MtValues[ielemId].size();
    sigmaPartial=0;
    sigmaTotal=0;
    sigmaFission = 0;
}
void TNudyFission::nFissionXsec(double nuEn, TNudyElement *targetElement)
{
     kineticE = nuEn;
     div_t divr;
     charge = targetElement->GetatomicNumber();
     mass   = targetElement->GetatomicMass();
     double sumXsec=0;
     //fout.open("/home/shiba/output/fissioneutron.txt",ios::out);
     for(unsigned int crsp = 0; crsp < mtValues; crsp++) {
        MT  = rp->MtValues[ielemId][crsp];
        MF4 = rp->GetMt4(ielemId, MT);
        MF5 = rp->GetMt5(ielemId, MT);
        MF6 = rp->GetMt6(ielemId, MT);
        LCT = rp->GetCos4Lct(ielemId, MT);
        //if(MT== 18){// || MT == 454 || MT == 459){
        switch (MT) {
               case 18:{
        
            sigmaFission = rp->GetSigmaPartial(ielemId, crsp, kineticE);
            sigmaTotal   = rp->GetSigmaTotal(ielemId, kineticE);
            procNameFis  = "nFission";
            nameF        =  "Fission fragements and others";
            double zaf   = rp->GetFisYield(ielemId, kineticE);
      	          divr         = div(zaf, 10000);
                  fissionfragmentsCharge = divr.quot;//fission fragment charge
                  divr         = div(divr.rem, 10);
            fissionFragmentmass1 = divr.quot;
            massFissionFragments.push_back(fissionFragmentmass1);
            chargeFissionFragments.push_back(fissionfragmentsCharge);
            chargeFissionFragments.push_back(charge - fissionfragmentsCharge);
            
            double nut = rp->GetNuTotal(ielemId, kineticE);
            // cout<<"B4 neutron fission part:\t"<<nu<<"\t"<<nut<<endl;
            double kineticRand = kineticE;// * fRnd->Uniform(1);
                      fissHeat = rp->GetFissHeat(ielemId, kineticRand);
                           nup = rp->GetNuPrompt(ielemId, kineticRand);
                           nud = rp->GetNuDelayed(ielemId, kineticRand);
            //cout<<"sigmaFission: "<<sigmaFission<<endl;
            //cout<<" B4 MC process  "<< nup <<"  "<< nud <<"  "<< nut <<endl;
            int nu                              = (int)nut;
            if(fRnd->Uniform(1) < nut - nu) nu = nu + 1;
              //cout<<"After MC process:\t"<<nup<<"  "<<nud<<"  "<<nut<<" "<<nu<<endl;
                for(int nui = 0; nui < nu; nui++) {
                    secEnergyLab = rp->GetEnergy5(ielemId, MT, kineticE);// energy of neutron produced from fission
                    cosLab       = 2 * fRnd->Uniform(1) - 1;//isotropic angle of neutron produced from fission
                    //cout<<nui<<"  "<<nu<<" prompt neutron energy : "<< secEnergyLab<<endl; 
                    secEnergyLabF.push_back(secEnergyLab);
                    cosLabF.push_back(cosLab);
                 }
                 fissionFragmentmass2 = mass - fissionFragmentmass1 - nu;
                 massFissionFragments.push_back(fissionFragmentmass2);
                //if(fissionFragmentmass1 == 0.0 || fissionFragmentmass2 == 0.0)cout<<mass<<"    "<<fissionFragmentmass1<<"    "<<fissionFragmentmass1<<endl; 
           }break;
        }
     }
    totalkineticEnergyFF = rp->GetkineticEnergyFF(ielemId,kineticE);  
    kineticEnergyFragment1 = totalkineticEnergyFF*fissionFragmentmass2/(fissionFragmentmass1 + fissionFragmentmass2);  
    kineticEnergyFragment2 = totalkineticEnergyFF - kineticEnergyFragment1;
    //cout<<"total KE FFs: "<<totalkineticEnergyFF<< " KE fragm1:"<< kineticEnergyFragment1<< " KE fragm2: "<<kineticEnergyFragment2<<endl;
     kineticEnergyFissionFragments.push_back(kineticEnergyFragment1);
     kineticEnergyFissionFragments.push_back(kineticEnergyFragment2);
}
//______________________________________________________________________________     
std::string TNudyFission::GetFissionProcessName()
{ 
     return procNameFis;
}
//______________________________________________________________________________
double TNudyFission::GetFissionXsec()
{ 
     return sigmaFission;
} 
//______________________________________________________________________________
std::string TNudyFission::GetParticleFission()
{ 
     return nameF;
}
//______________________________________________________________________________
std::vector<double> TNudyFission::GetKineticEnergyNeutron()
{ 
     return secEnergyLabF;
}
//______________________________________________________________________________
std::vector<double> TNudyFission::GetcosAngNeutron()
{ 
     return cosLabF;
}
//______________________________________________________________________________
std::vector<int> TNudyFission::GetFissionFragmentsmass()
{ 
     return massFissionFragments;
}
//______________________________________________________________________________
std::vector<int> TNudyFission::GetFissionFragmentscharge()
{ 
     return chargeFissionFragments;
}
//______________________________________________________________________________
std::vector<double> TNudyFission::GetKineticEnergyFissionFragments()
{
     return kineticEnergyFissionFragments;
}
//______________________________________________________________________________
double TNudyFission::GetPromptneutron()
{ 
     return nup;
}
//______________________________________________________________________________
void TNudyFission::processReset()
{
massFissionFragments.clear();
chargeFissionFragments.clear();
kineticEnergyFissionFragments.clear();
secEnergyLabF.clear();
cosLabF.clear();
}
//______________________________________________________________________________
TNudyFission::~TNudyFission()
{
  delete rp;
}

