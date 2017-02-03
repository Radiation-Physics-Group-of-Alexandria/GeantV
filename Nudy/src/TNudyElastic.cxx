//SPB
//Elastic process of neutron, provides cross-section at temprature of 293.606087K , energy and angle in Lab Frame
#include "TNudyCore.h"
#include "TNudyEndfRecoPoint.h"
#include "TNudyElastic.h"
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
ClassImp(TNudyElastic)
#endif

const char TNudyElastic::fkElNam[119][4] = {
        "n",  "H",  "He", "Li", "Be", "B",  "C",  "N",  "O",  "F",  "Ne", "Na",  "Mg",  "Al",  "Si",  "P",   "S",
        "Cl", "Ar", "K",  "Ca", "Sc", "Ti", "V",  "Cr", "Mn", "Fe", "Co", "Ni",  "Cu",  "Zn",  "Ga",  "Ge",  "As",
        "Se", "Br", "Kr", "Rb", "Sr", "Y",  "Zr", "Nb", "Mo", "Tc", "Ru", "Rh",  "Pd",  "Ag",  "Cd",  "In",  "Sn",
        "Sb", "Te", "I",  "Xe", "Cs", "Ba", "La", "Ce", "Pr", "Nd", "Pm", "Sm",  "Eu",  "Gd",  "Tb",  "Dy",  "Ho",
        "Er", "Tm", "Yb", "Lu", "Hf", "Ta", "W",  "Re", "Os", "Ir", "Pt", "Au",  "Hg",  "Tl",  "Pb",  "Bi",  "Po",
        "At", "Rn", "Fr", "Ra", "Ac", "Th", "Pa", "U",  "Np", "Pu", "Am", "Cm",  "Bk",  "Cf",  "Es",  "Fm",  "Md",
        "No", "Lr", "Rf", "Db", "Sg", "Bh", "Hs", "Mt", "Ds", "Rg", "Cn", "Uut", "Uuq", "Uup", "Uuh", "Uus", "Uuo"};

TNudyElastic::TNudyElastic(){}

TNudyElastic::TNudyElastic(int ElemId, const char* rootfileName)
{
    ielemId =  ElemId;// for element number like say for H=1, U=2, the way we number the elements for materials having composite elements
    const char* irENDF=rootfileName;//"/home/shiba/geant/rootendf/O16.root";
    rp = new TNudyEndfRecoPoint(ielemId, irENDF);
    mtValues = rp->MtValues[ielemId].size();
    sigmaPartial=0;
    sigmaTotal=0;
    sigmaElastic = 0;
}

void TNudyElastic::nElasticXsec(double nuEn, TNudyElement *targetElement)
{
     kineticE = nuEn;
     charge = targetElement->GetatomicNumber();
     mass   = targetElement->GetatomicMass();
     //cout<<"Target of mass---->"<<mass<<"---charge---"<<charge<<"\t"<<kineticE<<endl;
     double gamma = 0;
     double massRatio = 0;
     double projectileMass = 1.0; // for neutron 
     for(unsigned int crsp = 0; crsp < mtValues; crsp++) {
        MT  = rp->MtValues[ielemId][crsp];
        MF4 = rp->GetMt4(ielemId, MT);
        MF5 = rp->GetMt5(ielemId, MT);
        MF6 = rp->GetMt6(ielemId, MT);
        LCT = rp->GetCos4Lct(ielemId, MT);
        div_t divr;
        //if(MT == 2){
        switch (MT) {
               case 2:{
          sigmaElastic    = rp->GetSigmaPartial(ielemId, crsp, kineticE);
          sigmaTotal      = rp->GetSigmaTotal(ielemId, kineticE);
          processName     = "nElastic";
          secParticleName = "n";
          residueAEl      = mass;// target mass
          residueZEl      = charge;
          cosCME          = rp->GetCos4(ielemId, MT, kineticE);
          cosLabE         = TNudyCore::Instance()->cmToLabElasticCosT(cosCME, mass);
          secondarycosAng = TNudyCore::Instance()->cmToLabElasticCosT(cosCME, mass);
          secEnergyLabE   = TNudyCore::Instance()->cmToLabElasticE(kineticE, cosCME, mass);
          
          massRatio = mass/projectileMass;
          gamma     = 1.0;//massRatio;
          secondarycosAng2 = (1 - gamma*cosCME)/sqrt((gamma*gamma) + 1 - (2*gamma*cosCME)); 
          
          energyProductsLab.push_back(secEnergyLabE);
          energyProductsLab.push_back(kineticE - secEnergyLabE);
          cosAngProductsLab.push_back(cosLabE);
          cosAngProductsLab.push_back(secondarycosAng2);
          
          std::cout<<secEnergyLabE<<"   Residue energy in lab: "<<kineticE - secEnergyLabE<<"  "<< kineticE*(mass/pow((mass+1),2.0))*2*(1-cosCME)<<std::endl;
          
          residueType = residueNameElastic(residueZEl, residueAEl);
          productsName.push_back(secParticleName);
          productsName.push_back(residueType);
          
           }break;
         } 
     }
 }
//______________________________________________________________________________
std::string TNudyElastic::residueNameElastic(int Z, int A)
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
std::string TNudyElastic::GetElasticProcessName()
{ 
     return processName;
}
//______________________________________________________________________________
double TNudyElastic::GetElasticXsec()
{ 
       return sigmaElastic;
}    
//______________________________________________________________________________
std::vector<std::string> TNudyElastic::GetParticleElastic()
{ 
     return productsName;
}
//______________________________________________________________________________
std::vector<double> TNudyElastic::GetsecParticleKiEn()
{ 
       return energyProductsLab;
}
//______________________________________________________________________________
std::vector<double> TNudyElastic::GetsecParticlecosAng()
{ 
      return  cosAngProductsLab;
}
//______________________________________________________________________________
TNudyElastic::~TNudyElastic()
{
  delete rp;
}
//______________________________________________________________________________
void TNudyElastic::nElasticProcessReset()
{
  cosAngProductsLab.clear();
  energyProductsLab.clear();
  productsName.clear();
}

