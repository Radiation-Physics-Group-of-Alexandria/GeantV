//SPB
//Neutron Capture process
//The energy and angle of out going gamma is not yet set. The default values are set zero
//From this class, we are getting the neutron capture cross-section at temprature of 293.606087K.
#include "TNudyCore.h"
#include "TNudyEndfRecoPoint.h"
#include "TNudyCapture.h"
#include "TNudyEndfMat.h"
#include "TNudyElement.h"
#include <iostream>
#include <string>
#include <fstream>
#include <sstream>
#ifdef USE_ROOT
#include "TRandom3.h"
#endif
#include "TCanvas.h"
#include "TFile.h"
using namespace std;
#ifdef USE_ROOT
ClassImp(TNudyCapture)
#endif

const char TNudyCapture::fkElNam[119][4] = {
        "n",  "H",  "He", "Li", "Be", "B",  "C",  "N",  "O",  "F",  "Ne", "Na",  "Mg",  "Al",  "Si",  "P",   "S",
        "Cl", "Ar", "K",  "Ca", "Sc", "Ti", "V",  "Cr", "Mn", "Fe", "Co", "Ni",  "Cu",  "Zn",  "Ga",  "Ge",  "As",
        "Se", "Br", "Kr", "Rb", "Sr", "Y",  "Zr", "Nb", "Mo", "Tc", "Ru", "Rh",  "Pd",  "Ag",  "Cd",  "In",  "Sn",
        "Sb", "Te", "I",  "Xe", "Cs", "Ba", "La", "Ce", "Pr", "Nd", "Pm", "Sm",  "Eu",  "Gd",  "Tb",  "Dy",  "Ho",
        "Er", "Tm", "Yb", "Lu", "Hf", "Ta", "W",  "Re", "Os", "Ir", "Pt", "Au",  "Hg",  "Tl",  "Pb",  "Bi",  "Po",
        "At", "Rn", "Fr", "Ra", "Ac", "Th", "Pa", "U",  "Np", "Pu", "Am", "Cm",  "Bk",  "Cf",  "Es",  "Fm",  "Md",
        "No", "Lr", "Rf", "Db", "Sg", "Bh", "Hs", "Mt", "Ds", "Rg", "Cn", "Uut", "Uuq", "Uup", "Uuh", "Uus", "Uuo"};


TNudyCapture::TNudyCapture(){}

TNudyCapture::TNudyCapture(int ElemId, const char* rootfileName)
{
    ielemId =  ElemId;// for element number like say for H=1, U=2, the way we number the elements for materials having composite elements
    
    //elemId=ElemId;
    fRnd     = new TRandom3(0);
    const char* irENDF=rootfileName;//"/home/shiba/geant/rootendf/O16.root";
    rp = new TNudyEndfRecoPoint(ielemId, irENDF);
    mtValues = rp->MtValues[ielemId].size();
    sigmaPartial=0;
    sigmaTotal=0;
    sigmaCapture = 0;
}

void TNudyCapture::nCaptureXsec(double nuEn, TNudyElement *targetElement)
{
     kineticE = nuEn;
     charge = targetElement->GetatomicNumber();
     mass   = targetElement->GetatomicMass();
     for(unsigned int crsp = 0; crsp < mtValues; crsp++) {
        MT  = rp->MtValues[ielemId][crsp];
        switch (MT) {
               case 102:{
        
           MF4 = rp->GetMt4(ielemId, MT);
           MF5 = rp->GetMt5(ielemId, MT);
           MF6 = rp->GetMt6(ielemId, MT);
           LCT = rp->GetCos4Lct(ielemId, MT);
           
           std::cout<<"MF4: "<<MF4<<" MF5: "<<MF5<<"  MF6:"<<MF6<<"  LCT: "<<LCT<<std::endl;
           
          residueACa = mass + 1;
          residueZCa = charge;
          sigmaCapture = rp->GetSigmaPartial(ielemId, crsp, kineticE);
          sigmaTotal   = rp->GetSigmaTotal(ielemId, kineticE);
          procNameCap  = "nCapture"; 
          name        =  "gamma";
          cosLab = 0;
          secEnergyCM = 0;
          secEnergyLab = 0;
          std::cout<<"capture x-section: "<<sigmaCapture<<std::endl;
          std::cout<<"***********************************************************"<<std::endl;
          std::cout<<"*"<<"                     IMPORTANT NOTE                      *"<<std::endl;
          std::cout<<"*"<<" In nCapture process the angle of out going gamma is not *"<<std::endl;
          std::cout<<"*"<<" yet set: assumed fixed value is considered              *"<<std::endl;
          std::cout<<"***********************************************************"<<std::endl;
          productsName.push_back("gamma");
          residueType = residueName(residueZCa, residueACa);
          productsName.push_back(residueType);
          residueKineticEnergy = 0.0;
          residuecosAng = 0.0;
          GetSecParameters(targetElement, rp);
          nCaptureProductsEnergy.push_back(secEnergyLab);
          nCaptureProductsEnergy.push_back(residueKineticEnergy);
          nCaptureProductscosAng.push_back(cosLab);
          nCaptureProductscosAng.push_back(residuecosAng);
          }break;
        }
      }
 }       
//______________________________________________________________________________
void TNudyCapture::GetSecParameters(TNudyElement *targetElement, TNudyEndfRecoPoint *recoPoint)
{
// Here fe3  : photon energy, 
// Incident particle energy and Q-value are given in eV

    double vc = 3*pow(10.0,8.0);
    double vc2 = vc*vc;
     double u = 931.494;//MeV/C^2;
    double fm1 = 1.0 * u; // mass of projectile as neutron
    double fm2 = mass*u; // mass of target
    double fm4 = (mass + 1.0)*u; // mass of residue
    double fm3 =0.0;// (fm1 + fm2 - fm4); // mass of secondary particle
    double fqval =  1.0E-6*recoPoint->GetQValue(ielemId, MT) ;
    double cosAng = 0.1 ; // cos(theta) of emitted gamma in Lab 
    //std::cout<<"fm1: "<<fm1<< "   fm2 "<<fm2<<" fm3 : "<<fm3<< " fm4: "<<fm4<<std::endl;
    double phiResidue; // Angle of residue
    double gammaEnergy;
    double recoilEnergy;
    double sinAng = sqrt(1- cosAng*cosAng);        
    double fact1,fact2;      
        
        kineticE = kineticE*1.0E-6;  
        gammaEnergy = (fqval*(fm1+fm2+fm4)*0.5 + fm2*kineticE)/(fm1+fm2 + kineticE - cosAng*sqrt(kineticE*(2*fm1+kineticE)));
        //cout<<"Energy of gamma : "<<gammaEnergy<<endl;
        //thermal neutron capture kinematics
        //Egamma = (fm4)*(sqrt(1+(2*fqval/fm4)) -1);
        recoilEnergy = gammaEnergy*gammaEnergy/(2*fm4);
        //cout<<"Energy of gamma in MeV : "<<gammaEnergy<<" Recoil energy in MeV: "<<recoilEnergy<< " " << gammaEnergy + recoilEnergy<<endl;
        fact1 = 1.0/gammaEnergy;
        fact2 = fact1* sqrt(2*fm1*kineticE);
        phiResidue = atan(sinAng/(fact2 - cosAng));
        cosLab = cosAng;
        secEnergyLab = gammaEnergy*1.0E6; // converting MeV to eV
        residueKineticEnergy = recoilEnergy*1.0E6; // converting MeV to eV
        residuecosAng = cos(phiResidue);
   
}
//______________________________________________________________________________     
std::string TNudyCapture::GetCaptureProcessName()
{ 
     return procNameCap;
}
//______________________________________________________________________________
double TNudyCapture::GetCaptureXsec()
{ 
       return sigmaCapture;
} 
//______________________________________________________________________________
std::vector<std::string> TNudyCapture::nCaptureGetsecParticles()
{ 
     return productsName;
}
//______________________________________________________________________________
std::vector<double> TNudyCapture::nCaptureGetEnergy(){ 
       return nCaptureProductsEnergy;
  }
//______________________________________________________________________________
std::vector<double> TNudyCapture::nCaptureGetcosAng()
{ 
      return nCaptureProductscosAng;
}
//_____________________________________________________________________________
std::string TNudyCapture::residueName(int Z, int A)
{
  
  stringstream ssZ;
  ssZ << Z;
  string strZ = ssZ.str();
  
  stringstream ssA;
  ssA << A;
  string strA = ssA.str();
  
  std::string nameRes;
  std::string finalName;
  nameRes = fkElNam[Z];
  finalName = strZ + nameRes + strA;
  return finalName;
}
//______________________________________________________________________________
void TNudyCapture::nCaptureProcessReset()
{
  nCaptureProductsEnergy.clear();
  nCaptureProductscosAng.clear();
  productsName.clear();
}
//______________________________________________________________________________
TNudyCapture::~TNudyCapture()
{
delete fRnd;
delete rp;
}
//______________________________________________________________________________

