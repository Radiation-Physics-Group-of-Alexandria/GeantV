// Various interaction process of neutron: Selection of secondary particle type, its energy and angle
//The final state product residues mass and charges are assigned considering the incident particle
//as a neutron. The charge and mass of the residue will change depending upon incident particle as
//a charge particle or photon.
#include <iostream>
//#include "TNudyEndfDoppler.h"
//#include "TNudyEndfAng.h"
//#include "TNudyEndfEnergy.h"
//#include "TNudyEndfEnergyAng.h"
//#include "TNudyEndfFissionYield.h"
#include <fstream>
#include "ElementProp.h"
#include "TNudyCore.h"
#include "TNudyEndfRecoPoint.h"
#include "TNudyFission.h"
#include "TNudyEndfMat.h"
#include "TNudyElement.h"
#ifdef USE_ROOT
#include "TRandom3.h"
#endif
#include "TCanvas.h"
#include "TFile.h"
using namespace std;
#ifdef USE_ROOT
ClassImp(TNudyFission)
#endif

TNudyFission::TNudyFission(){
   
}

TNudyFission::TNudyFission(int ElemId, const char* rootfileName)
{
    ielemId =  ElemId;// for element number like say for H=1, U=2, the way we number the elements for material having composite elements
    
    elemId=ElemId;
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
     double sumXsec=0;
     
     //fout.open("/home/shiba/output/fissioneutron.txt",ios::out);
     for(unsigned int crsp = 0; crsp < mtValues; crsp++) {
        
        MT  = rp->MtValues[ielemId][crsp];
        MF4 = rp->GetMt4(ielemId, MT);
        MF5 = rp->GetMt5(ielemId, MT);
        MF6 = rp->GetMt6(ielemId, MT);
        LCT = rp->GetCos4Lct(ielemId, MT);
        if(MT== 18){// || MT == 454 || MT == 459){
            sigmaFission = rp->GetSigmaPartial(ielemId, crsp, kineticE);
            sigmaTotal   = rp->GetSigmaTotal(ielemId, kineticE);
            procNameFis  = "nFission";
            nameF        =  "Fission fragements and others";
            double zaf   = rp->GetFisYield(ielemId, kineticE);
      	    divr         = div(zaf, 1000);
            residueZFs     = divr.quot;//fission fragment mass
            residueAFs = divr.rem;//fission fragment charge
            double nut = rp->GetNuTotal(ielemId, kineticE);
            // cout<<"B4 neutron fission part:\t"<<nu<<"\t"<<nut<<endl;
            double kineticRand = kineticE;// * fRnd->Uniform(1);
            double fissHeat = rp->GetFissHeat(ielemId, kineticRand);
                        nup = rp->GetNuPrompt(ielemId, kineticRand);
            double nud = rp->GetNuDelayed(elemId, kineticRand);
           
           cout<<"sigmaFission: "<<sigmaFission<<endl;
           
           //fout<<" sample heat "<< fissHeat/1E6 <<"  "<< nut <<"  "<< nup <<"  "<< nud <<"  "<<
            int nu                              = (int)nut;
            if(fRnd->Uniform(1) < nut - nu) nu = nu + 1;
              if(fRnd->Uniform(1) < nut - nu) nu = nu + 1;
               //cout<<"After neutron fission part:\t"<<nu<<"\t"<<nut<<endl;
                 for(int nui = 0; nui < nu; nui++) {
                    secEnergyLabF = rp->GetEnergy5(ielemId, MT, kineticE);// energy of neutron produced from fission
                    cosLabF       = 2 * fRnd->Uniform(1) - 1;//isotropic angle of neutron produced from fission
                    }
                }
          }
  }
//________________________________________________________________________________________     

std::string TNudyFission::GetFissionProcessName()
{ 
     return procNameFis;
}


double TNudyFission::GetFissionXsec()
{ 
       return sigmaFission;
} 


std::string TNudyFission::GetParticleFission()
{ 
     return nameF;
}

//________________________________________________________________________________________
double TNudyFission::GetKiEnFission(){ 
       return secEnergyLabF;
  }
//________________________________________________________________________________________
double TNudyFission::GetcosAngFission()
{ 
      return cosLabF;
}

double TNudyFission::GetFissionFragmentmass()
{ 
      return residueAFs;
}

double TNudyFission::GetPromptneutron()
{ 
      return nup;
}


//_________________________________________________________________________________________
//LCT Flag to specify the frame of reference used
//LCT=1, the data are given in the LAB system
//LCT=2, the data are given in the CM system
//MF=File number
//MT=Reaction number due to various reaction type
void TNudyFission::GetSecParameter(TNudyElement *targetElement, TNudyEndfRecoPoint *recoPoint)
{

   double massT= targetElement->GetatomicMass();
  // cout<<"Mass of Target:::::::::::::::\t"<<massT<<endl;
   if (MF4 == 4 && MF5 == 5) {
        //cout<<"angle in CM:\t"<<recoPoint->GetCos4(ielemId, MT, kineticE)<<endl;
        cosCM        = recoPoint->GetCos4(ielemId, MT, kineticE);
        secEnergyCM  = recoPoint->GetEnergy5(ielemId, MT, kineticE);
       //cout<< "CM cos(angle) and energy:\t"<<cosCM<<"\t"<<secEnergyCM<<endl;
        secEnergyLab = TNudyCore::Instance()->cmToLabInelasticE(secEnergyCM, kineticE, cosCM, massT);
       if (LCT == 2) {
          cosLab = TNudyCore::Instance()->cmToLabInelasticCosT(secEnergyLab, secEnergyCM, kineticE, cosCM,
                                                           massT);
          }else if (LCT == 1) {
           cosLab       = cosCM;
           secEnergyLab = secEnergyCM;
          }
    }else if (MF4 == 4 && MF5 == -1) {
    cosCM      = recoPoint->GetCos4(ielemId, MT, kineticE);
    double fm1 = 1;
    double fm2 = massT;
    double fm4 = residueAFs;
    double fm3 = massT - residueAFs;
  //cout<<"Mass of particles:"<<fm1<<"\t"<<fm2<<"\t"<<fm3<<"\t"<<fm4<<"\t"<<MT<<endl;
  // double fqval =  recoPoint->GetQValue(ielemId, MT) ;
    double fsi  = fm1 + fm2;
    double fdi  = fm1 - fm2;
    double fsf  = fm3 + fm4;
    double fdf  = fm3 - fm4;
    double fs   = fsi * fsi + 2 * fm2 * kineticE;
    double fki2 = (fs - fsi * fsi) * (fs - fdi * fdi) / (4 * fs);
    double fkf2 = (fs - fsf * fsf) * (fs - fdf * fdf) / (4 * fs);
    double fz2  = sqrt(fki2 + fm2 * fm2);
    double fz3  = sqrt(fkf2 + fm3 * fm3);
    // double fz4 = sqrt(fkf2 + fm4 * fm4);
    double fe3 = (fz2 * fz3 - sqrt(fki2 * fkf2) / fm2) - fm3;
    // double fe4 = (fz2 * fz4 - sqrt(fki2 * fkf2)/fm2) - fm4;
    // double fp12 = kineticE * kineticE + 2 * fm1 * kineticE;
    // double fp32 = fe3 * fe3 + 2 * fm3 * fe3;
    // double fp42 = fe4 * fe4 + 2 * fm4 * fe4;
    // double fcos3 = ((kineticE + fsi)*(fe3 + fm3)- fz3 * sqrt(fs))/sqrt(fp12*fp32);
    // double fcos4 = ((kineticE + fsi)*(fe4 + fm4)- fz4 * sqrt(fs))/sqrt(fp12*fp42);
   //if(cosCM>1)std::cout<<"cos CM = "<< cosCM << std::endl;
    // secEnergyCM = recoPoint->GetEnergy5(ielemId, MT, kineticE);
    secEnergyLab = fe3;
    if (LCT == 2) {
      cosLab = TNudyCore::Instance()->cmToLabInelasticCosT(secEnergyLab, secEnergyCM, kineticE, cosCM,
                                                           massT);
     } else if (LCT == 1) {
       cosLab = cosCM;
     }
   }else if (MF4 == 99 && MF5 == 5) {
    double fm1 = 1;
    double fm2 = massT;
    double fm4 = residueAFs;
    double fm3 = massT - residueAFs;
    // double fqval =  recoPoint->GetQValue(ielemId, MT) ;
    double fsi  = fm1 + fm2;
    double fdi  = fm1 - fm2;
    double fsf  = fm3 + fm4;
    double fdf  = fm3 - fm4;
    double fs   = fsi * fsi + 2 * fm2 * kineticE;
    double fki2 = (fs - fsi * fsi) * (fs - fdi * fdi) / (4 * fs);
    double fkf2 = (fs - fsf * fsf) * (fs - fdf * fdf) / (4 * fs);
    double fz2  = sqrt(fki2 + fm2 * fm2);
    double fz3  = sqrt(fkf2 + fm3 * fm3);
    // double fz4 = sqrt(fkf2 + fm4 * fm4);
    double fe3 = (fz2 * fz3 - sqrt(fki2 * fkf2) / fm2) - fm3;
    // double fe4 = (fz2 * fz4 - sqrt(fki2 * fkf2)/fm2) - fm4;
    double fp12 = kineticE * kineticE + 2 * fm1 * kineticE;
    double fp32 = fe3 * fe3 + 2 * fm3 * fe3;
    // double fp42 = fe4 * fe4 + 2 * fm4 * fe4;
    cosLab = ((kineticE + fsi) * (fe3 + fm3) - fz3 * sqrt(fs)) / sqrt(fp12 * fp32);
    // double fcos4 = ((kineticE + fsi)*(fe4 + fm4)- fz4 * sqrt(fs))/sqrt(fp12*fp42);
    secEnergyLab = fe3;
  } else if (MF4 == 99 && MF5 == -1 && MF6 == 6) {
    int law = recoPoint->GetLaw6(ielemId, MT);
    // std::cout<<"law "<< law << std::endl;
    switch (law) {
    case 2: {
      cosCM         = recoPoint->GetCos64(ielemId, MT, kineticE);
      double fm1    = 1;
      double fm2    = massT;
      double fm3    = massT - residueAFs;
      double fm4    = residueAFs;
      double fqval  = recoPoint->GetQValue(ielemId, MT);
      double fEt    = kineticE + fqval;
      double fm1m2  = fm1 + fm2;
      double fm3m4  = fm3 + fm4;
      double fA1    = fm1 * fm4 * (kineticE / fEt) / (fm1m2 * fm3m4);
      double fB1    = fm1 * fm3 * (kineticE / fEt) / (fm1m2 * fm3m4);
      double fC1    = fm2 * fm3 * (1 + fm1 * fqval * fm4 / (fEt * fm1)) / (fm1m2 * fm3m4);
      double fD1    = fm2 * fm4 * (1 + fm1 * fqval * fm4 / (fEt * fm1)) / (fm1m2 * fm3m4);
      secEnergyLab  = fEt * (fB1 + fD1 + 2 * sqrt(fA1 * fC1) * cosCM);
      double sinLab = fD1 * cosCM * fEt / secEnergyLab;
      cosLab        = sqrt(1 - sinLab * sinLab);
      // std::cout<< cosLab <<"  "<<secEnergyLab << std::endl;
    } break;
    case 3:
      break;
    case 1:
    case 6:
      cosCM = recoPoint->GetCos6(ielemId, MT, kineticE);
      // std::cout<<MT<<"\tcos:\t "<< cosCM << std::endl;
      secEnergyCM = recoPoint->GetEnergy6(ielemId, MT, kineticE);
      // std::cout<<"E "<< secEnergyCM<< std::endl;
      secEnergyLab = TNudyCore::Instance()->cmToLabInelasticE(secEnergyCM, kineticE, cosCM, massT);
      cosLab       = TNudyCore::Instance()->cmToLabInelasticCosT(secEnergyLab, secEnergyCM, kineticE, cosCM,
                                                           massT);
      break;
    case 7:
      cosLab = recoPoint->GetCos6(ielemId, MT, kineticE);
      // std::cout<< cosLab << std::endl;
      secEnergyLab = recoPoint->GetEnergy6(ielemId, MT, kineticE);
      // std::cout<< secEnergyLab << std::endl;
      break;
    }
  }
// } 
  
  
}
void TNudyFission::FillHisto(double icosLab, double isecEnergyLab)
{
 // h->Fill(icosLab, isecEnergyLab / 1E9);
  //h1->Fill(icosLab);
  //h2->Fill(isecEnergyLab / 1E9);
  // x[ecounter] = isecEnergyLab/1E9 ;
  // y[ecounter] = icosLab ;
  // if(events<=1000000)ecounter ++ ;
}
//
//_________________________________________________________________________________________________________________
TNudyFission::~TNudyFission()
{
 ProcName.clear();
 xec.clear();

delete fRnd;
delete rp;

 }

