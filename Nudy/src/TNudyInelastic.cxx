//SPB
// Various inelastic process of neutron: Selection of secondary particle type, its energy and angle
//The final state product residues mass and charges are assigned considering the incident particle
//as a neutron. The charge and mass of the residue will change depending upon incident particle as
//a charge particle or photon.
//NOTE: For some of the neutron inelastic reaction channel although the cross section available but secondary products energy and angle informatio are not available.
//For some of the target nuclei the endf data files are not available. 
//
//Secondary products are available for some of the reaction channel MT-16 = 2n, MT-17 = 3n, MT-22 = n+alpha, MT- 23 = n+3alpha, MT- 24 =  2n+alpha, MT-25 = 3n+alpha , MT-28 = n+p, MT-29 = n+2alpha,
// MT-50-91 = n', MT-103 = p,
//MT-104 = d, MT-105 = t, MT-107 = alpha, MT-108 = 2alpha, MT-11= 2p, Mt -112 = p+alpha, MT-115 = p+d, MT-117 = d+alpha,

// From this class, we are getting the total neutron inelastic cross-section at temprature of 293.606087K.

//In some of inelastic process, the energy and angle information for neutron, charge particles, residues as well as gamma are given in File-6 of ENDF data file. In this class, I have not consider the gamma information, but that has to be considered
//_________________________________________________________________________________________ 
#include "TNudyEndfEnergyAng.h"
#include "TNudyEndfFile.h"
#include "TNudyCore.h"
#include "TNudyEndfRecoPoint.h"
#include "TNudyInelastic.h"
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
ClassImp(TNudyInelastic)
#endif


const char TNudyInelastic::fkElNam[119][4] = {
        "n",  "H",  "He", "Li", "Be", "B",  "C",  "N",  "O",  "F",  "Ne", "Na",  "Mg",  "Al",  "Si",  "P",   "S",
        "Cl", "Ar", "K",  "Ca", "Sc", "Ti", "V",  "Cr", "Mn", "Fe", "Co", "Ni",  "Cu",  "Zn",  "Ga",  "Ge",  "As",
        "Se", "Br", "Kr", "Rb", "Sr", "Y",  "Zr", "Nb", "Mo", "Tc", "Ru", "Rh",  "Pd",  "Ag",  "Cd",  "In",  "Sn",
        "Sb", "Te", "I",  "Xe", "Cs", "Ba", "La", "Ce", "Pr", "Nd", "Pm", "Sm",  "Eu",  "Gd",  "Tb",  "Dy",  "Ho",
        "Er", "Tm", "Yb", "Lu", "Hf", "Ta", "W",  "Re", "Os", "Ir", "Pt", "Au",  "Hg",  "Tl",  "Pb",  "Bi",  "Po",
        "At", "Rn", "Fr", "Ra", "Ac", "Th", "Pa", "U",  "Np", "Pu", "Am", "Cm",  "Bk",  "Cf",  "Es",  "Fm",  "Md",
        "No", "Lr", "Rf", "Db", "Sg", "Bh", "Hs", "Mt", "Ds", "Rg", "Cn", "Uut", "Uuq", "Uup", "Uuh", "Uus", "Uuo"};
TNudyInelastic::TNudyInelastic(){
}
TNudyInelastic::TNudyInelastic(int ElemId, const char* rootfileName)
{
  ielemId =  ElemId;
  //for element number like say for H=1, U=2, the way we
  //number the elements for materials having composite elements
  //elemId=ElemId;
  fRnd     = new TRandom3(0);
  const char* irENDF=rootfileName;//"/home/shiba/geant/rootendf/O16.root";
  rp = new TNudyEndfRecoPoint(ielemId, irENDF);
  nMTs = rp->MtValues[ielemId].size();
  sigmaPartial=0;
  sigmaTotal=0;
  sigmaInelastic = 0;
}
void TNudyInelastic::nInelasticXsec(double nuEn, TNudyElement *targetElement)
{
  //std::cout<<"product elemets: "<<recoEnergyAng->GetZd6V.size()<<std::endl;
  kineticE = nuEn;
  int reactionMT;
  double sumXsec=0;
  charge = targetElement->GetatomicNumber(); // target charge
  mass   = targetElement->GetatomicMass(); //target mass
  projectileMass = 1.0; // Projectile is neutron
  projectileCharge = 0.0; // Projectile is neutron;
 
  int count1=0;
  for(unsigned int crsp = 0; crsp < nMTs; crsp++){
    reactionMT=rp->MtValues[ielemId][crsp];
    //cout<<reactionMT<<"   n-inelastic partial cross-section: "<<rp->GetSigmaPartial(ielemId, crsp, kineticE)<<endl;
    if(reactionMT != 2 && reactionMT !=18 && reactionMT !=102 ){
      
      cout<<reactionMT<<"   n-inelastic partial cross-section: "<<rp->GetSigmaPartial(ielemId, crsp, kineticE)<<endl;
      
      procNameInelastic   = "nInelastic";
      name        =  "2n";
      sumXsec=sumXsec+rp->GetSigmaPartial(ielemId, crsp, kineticE);
      sigmaInelastic=sumXsec; //total n-inelastic cross-section
      //cout<<"inelastic xec: "<<sumXsec<<endl;
      mtTemp.push_back(reactionMT);
      crssIel.push_back(rp->GetSigmaPartial(ielemId, crsp, kineticE));
      nMTinelastic[ielemId][count1] =reactionMT;
      count1=count1+1;
      
    }
  }
  for (unsigned int iel = 0; iel < mtTemp.size(); iel++){
    crss.push_back(crssIel[iel]/sigmaInelastic);
  }
  double sum0=0;
  double rnd0 = fRnd->Uniform(1);
  for (unsigned int jel = 0; jel < mtTemp.size(); jel++){
    sum0 += crss[jel];
    if(rnd0 <= sum0){
      isel = jel;
      break;
      }
  }
  MT  = nMTinelastic[ielemId][isel];//rp->MtValues[ielemId][isel];
  MF4 = rp->GetMt4(ielemId, MT);
  MF5 = rp->GetMt5(ielemId, MT);
  MF6 = rp->GetMt6(ielemId, MT);
  LCT = rp->GetCos4Lct(ielemId, MT);
  //if(MT==16)
  //std::cout<<" nInelasticXsec:    "<<MT<<"  "<<MF4<<" "<<MF5<<" "<<MF6<<" "<<LCT<<std::endl;
  int productZ;
  int productA;
  double cosangle1=0;
  //cosangle1= rp->GetCos4(ielemId, MT, kineticE);
   if (MF4 == 99 && MF5 == -1 && MF6 == 6){
   rp->ResetRecoEnergyAng();
   tempproductsName.clear();
  }
  if((MF4 == 4 && MF5 == 5) || (MF4 == 4 && MF5 == -1) || (MF4 == 99 && MF5 == 5)){
  tempproductsName.clear();
  }
//______________________________________________________________________________
//O-16, MT = 16,22,23,28,32,41,44,45,108,112
//17Cl35, MT = 111
//Ca-40, MT = 29,108,111,115,117
//Ca-42, Mt = 115,117
//Ca43 , MT=24, also 108, 111, 112
//24Cr52, MT = 103 , 107
//40Zr95, MT = 105
//25Mn55, MT = 37 // large Q-value
//46Pd110, MT = 37 // large Q-value
//13Al-27, MT = 32,33
//28Ni59, MT = 34
//30Zn65, MT =16, 22, 24,28,32, 41
//35Br81, MT = 17 
//Li6 ,  MT = 24
//Li7 ,  MT= 25
// In some of the reaction channels the threshhold energy is greater than 20.0 MeV
switch (MT){
    //case 2:{
     // }break;
      //productMass, mass of products except neutron 
      case 11: // 2nd
      residueA = mass + projectileMass - 4;
      residueZ = charge + projectileCharge - 1;
      productsName.push_back("n");
      productsName.push_back("n");
      productsName.push_back("1H2");
      residueType = residueName(residueZ, residueA);
      productsName.push_back(residueType);
      neutronsN = 2;
      productMass.push_back(2.0);
      GetSecParameter(targetElement, rp);
      break;
      case 16: // 2n
      residueA = mass + projectileMass - 2;
      residueZ = charge + projectileCharge;
      // For some of the target materials, either energy(ENDF file-5) or angle information(ENDF file-4) for neutron given.
      //In that case 1st if condition is used
      if((MF4 == 4 && MF5 == 5) || (MF4 == 4 && MF5 == -1) || (MF4 == 99 && MF5 == 5)){
      tempproductsName.push_back("n");
      tempproductsName.push_back("n");
      residueType = residueName(residueZ, residueA);
      tempproductsName.push_back(residueType);
      // For some of the target materials, either energy or angle information for neutron, charge particles,residues and gamma( not considered in this class) is in File-6 of ENDF data file given.
      //In that case 2nd if condition is used
      }else if(MF4 == 99 && MF5 == -1 && MF6 == 6){
      tempproductsName.push_back("n");
      residueType = residueName(residueZ, residueA);
      tempproductsName.push_back(residueType);
      tempproductsName.push_back("n");
      }
      productsName = tempproductsName; 
      neutronsN = 2;
      productMass.push_back(0.0);
      GetSecParameter(targetElement, rp);
      break;
      case 17: // 3n
      residueA = mass + projectileMass - 3;
      residueZ = charge + projectileCharge;
      if((MF4 == 4 && MF5 == 5) || (MF4 == 4 && MF5 == -1) || (MF4 == 99 && MF5 == 5)){
      tempproductsName.push_back("n");
      tempproductsName.push_back("n");
      tempproductsName.push_back("n");
      residueType = residueName(residueZ, residueA);
      tempproductsName.push_back(residueType);
      }else if(MF4 == 99 && MF5 == -1 && MF6 == 6){
      tempproductsName.push_back("n");
      residueType = residueName(residueZ, residueA);
      tempproductsName.push_back(residueType);
      tempproductsName.push_back("n");
      tempproductsName.push_back("n");
      }
      productsName = tempproductsName; 
      neutronsN = 3;
      productMass.push_back(0.0);
      GetSecParameter(targetElement, rp);
      break;
      //case 18:{}
      //break;
      case 22: // n+alpha
      residueA = mass + projectileMass - 5;
      residueZ = charge + projectileCharge - 2;
      tempproductsName.push_back("n");
      tempproductsName.push_back("2He4");
      residueType = residueName(residueZ, residueA);
      tempproductsName.push_back(residueType);
      productsName = tempproductsName; 
      neutronsN = 1;
      productMass.push_back(4.0);
      GetSecParameter(targetElement, rp);
      break;
      case 23: // n+3a
      residueA = mass + projectileMass - 13;
      residueZ = charge + projectileCharge - 6;
      productsName.push_back("n");
      productsName.push_back("2He4");
      residueType = residueName(residueZ, residueA);
      productsName.push_back(residueType);
      productsName.push_back("2He4");
      productsName.push_back("2He4");
      neutronsN = 1;
      productMass.push_back(4.0);
      productMass.push_back(4.0);
      productMass.push_back(4.0);
      GetSecParameter(targetElement, rp);
      break;
      case 24: // 2n+a
      residueA = mass + projectileMass - 6;
      residueZ = charge + projectileCharge - 2;
      if((MF4 == 4 && MF5 == 5) || (MF4 == 4 && MF5 == -1) || (MF4 == 99 && MF5 == 5)){
      tempproductsName.push_back("n");
      tempproductsName.push_back("n");
      tempproductsName.push_back("2He4");
      residueType = residueName(residueZ, residueA);
      tempproductsName.push_back(residueType);
      }else if(MF4 == 99 && MF5 == -1 && MF6 == 6){
      tempproductsName.push_back("n");
      tempproductsName.push_back("2He4");
      residueType = residueName(residueZ, residueA);
      tempproductsName.push_back(residueType);
      tempproductsName.push_back("n");
      }
      productsName = tempproductsName;
      neutronsN = 2;
      productMass.push_back(4.0);
      GetSecParameter(targetElement, rp);
      break;
      case 25: // 3n+a // Target - 7Li
      residueA = mass + projectileMass - 7;
      residueZ = charge + projectileCharge - 2;
      if((MF4 == 4 && MF5 == 5) || (MF4 == 4 && MF5 == -1) || (MF4 == 99 && MF5 == 5)){
      tempproductsName.push_back("n");
      tempproductsName.push_back("n");
      tempproductsName.push_back("n");
      tempproductsName.push_back("2He4");
      residueType = residueName(residueZ, residueA);
      tempproductsName.push_back(residueType);
      }else if(MF4 == 99 && MF5 == -1 && MF6 == 6){
      tempproductsName.push_back("n");
      tempproductsName.push_back("2He4");
      residueType = residueName(residueZ, residueA);
      tempproductsName.push_back(residueType);
      tempproductsName.push_back("n");
      tempproductsName.push_back("n");
      }
      productsName = tempproductsName;
      neutronsN = 3;
      productMass.push_back(4.0);
      GetSecParameter(targetElement, rp);
      break;
      case 28:{ // n+p
      residueA = mass + projectileMass - 2;
      residueZ = charge + projectileCharge - 1;
      productsName.push_back("n");
      productsName.push_back("p");
      residueType = residueName(residueZ, residueA);
      productsName.push_back(residueType);
      neutronsN = 1;
      productMass.push_back(1.0);
      GetSecParameter(targetElement, rp);
      for(int i=0; i<secEnergyinLab.size();i++){
      //std::cout<<" from raection channel: "<<secEnergyinLab[i]<<std::endl;
      }
      }break;
      case 29: // n +2a
      residueA = mass + projectileMass - 9;
      residueZ = charge + projectileCharge - 4;
      productsName.push_back("n");
      productsName.push_back("2He4");
      residueType = residueName(residueZ, residueA);
      productsName.push_back(residueType);
      productsName.push_back("2He4");
      neutronsN = 1;
      productMass.push_back(4.0);
      productMass.push_back(4.0);
      GetSecParameter(targetElement, rp);
      break;
      case 30: // 2n+2a
      residueA = mass + projectileMass - 10;
      residueZ = charge + projectileCharge - 4;
      productsName.push_back("n");
      productsName.push_back("2He4");
      residueType = residueName(residueZ, residueA);
      productsName.push_back(residueType);
      productsName.push_back("n");
      productsName.push_back("2He4");
      neutronsN = 2;
      productMass.push_back(4.0);
      productMass.push_back(4.0);
      GetSecParameter(targetElement, rp);
      break;
      case 32: // n+d
      residueA = mass + projectileMass - 3;
      residueZ = charge + projectileCharge - 1;
      productsName.push_back("n");
      productsName.push_back("1H2");
      residueType = residueName(residueZ, residueA);
      productsName.push_back(residueType);
      neutronsN = 1;
      productMass.push_back(2.0);
      GetSecParameter(targetElement, rp);
      break;
      case 33: // n+t
      residueA = mass + projectileMass - 4;
      residueZ = charge + projectileCharge - 1;
      productsName.push_back("n");
      productsName.push_back("1H3");
      residueType = residueName(residueZ, residueA);
      productsName.push_back(residueType);
      neutronsN = 1;
      productMass.push_back(3.0);
      GetSecParameter(targetElement, rp);
      break;
      case 34: // n + He3
      residueA = mass + projectileMass - 4;
      residueZ = charge + projectileCharge - 2;
      productsName.push_back("n");
      productsName.push_back("2He3");
      residueType = residueName(residueZ, residueA);
      productsName.push_back(residueType);
      neutronsN = 1;
      productMass.push_back(3.0);
      GetSecParameter(targetElement, rp);
      break;
      case 35: // n+d+2a
      residueA = mass + projectileMass - 11;
      residueZ = charge + projectileCharge - 5;
      productsName.push_back("n");
      productsName.push_back("1H2");
      productsName.push_back("2He4");
      residueType = residueName(residueZ, residueA);
      productsName.push_back(residueType);
      productsName.push_back("2He4");
      neutronsN = 1;
      productMass.push_back(2.0);
      productMass.push_back(4.0);
      productMass.push_back(4.0);
      GetSecParameter(targetElement, rp);
      break;
      case 36: // n+t+2a
      residueA = mass + projectileMass - 12;
      residueZ = charge + projectileCharge - 5;
      productsName.push_back("n");
      productsName.push_back("1H3");
      productsName.push_back("2He4");
      residueType = residueName(residueZ, residueA);
      productsName.push_back(residueType);
      productsName.push_back("2He4");
      neutronsN = 1;
      productMass.push_back(3.0);
      productMass.push_back(4.0);
      productMass.push_back(4.0);
      GetSecParameter(targetElement, rp);
      break;
      case 37: // 4n
      residueA = mass + projectileMass - 4;
      residueZ = charge + projectileCharge;
      if((MF4 == 4 && MF5 == 5) || (MF4 == 4 && MF5 == -1) || (MF4 == 99 && MF5 == 5)){
      tempproductsName.push_back("n");
      tempproductsName.push_back("n");
      tempproductsName.push_back("n");
      tempproductsName.push_back("n");
      residueType = residueName(residueZ, residueA);
      tempproductsName.push_back(residueType);
      }else if(MF4 == 99 && MF5 == -1 && MF6 == 6){
      tempproductsName.push_back("n");
      residueType = residueName(residueZ, residueA);
      tempproductsName.push_back(residueType);
      tempproductsName.push_back("n");
      tempproductsName.push_back("n");
      tempproductsName.push_back("n");
      }
      productsName = tempproductsName;
      neutronsN = 4;
      productMass.push_back(0.0);
      GetSecParameter(targetElement, rp);
      break;
      case 41: // 2n+p
      residueA = mass + projectileMass - 3;
      residueZ = charge + projectileCharge - 1;
      productsName.push_back("n");
      productsName.push_back("n");
      productsName.push_back("p");
      residueType = residueName(residueZ, residueA);
      productsName.push_back(residueType);
      neutronsN = 2;
      productMass.push_back(1.0);
      GetSecParameter(targetElement, rp);
      break;
      case 42: // 3n+p
      residueA = mass + projectileMass - 4;
      residueZ = charge + projectileCharge - 1;
      productsName.push_back("n");
      productsName.push_back("n");
      productsName.push_back("n");
      productsName.push_back("p");
      residueType = residueName(residueZ, residueA);
      productsName.push_back(residueType);
      neutronsN = 3;
      productMass.push_back(1.0);
      GetSecParameter(targetElement, rp);
      break;
      case 44: // n+2p
      residueA = mass + projectileMass - 3;
      residueZ = charge + projectileCharge - 1;
      productsName.push_back("n");
      productsName.push_back("p");
      productsName.push_back("p");
      residueType = residueName(residueZ, residueA);
      productsName.push_back(residueType);
      neutronsN = 1;
      productMass.push_back(1.0);
      productMass.push_back(1.0);
      GetSecParameter(targetElement, rp);
      break;
      case 45: // n+p+a
      residueA = mass + projectileMass - 6;
      residueZ = charge + projectileCharge - 3;
      productsName.push_back("n");
      productsName.push_back("p");
      productsName.push_back("2He4");
      residueType = residueName(residueZ, residueA);
      productsName.push_back(residueType);
      neutronsN = 1;
      productMass.push_back(1.0);
      productMass.push_back(4.0);
      GetSecParameter(targetElement, rp);
      break;
      case 50:{//production of neutron and residue at the ground state
      //Note: for this channel incident particle should be different from neutron
      residueA = mass + projectileMass - 1;
      residueZ = charge + projectileCharge;
      productsName.push_back("n and residue at ground state");
      residueType = residueName(residueZ, residueA);
      productsName.push_back(residueType);
      neutronsN = 1;
      productMass.push_back(0.0);
      GetSecParameter(targetElement, rp);
      }break;
      case 51:{//production of neutron and residue at the 1st excited state, this will continue upto 90
      residueA = mass + projectileMass - 1;
      residueZ = charge + projectileCharge;
      productsName.push_back("n and residue at 1st excited state");
      residueType = residueName(residueZ, residueA);
      productsName.push_back(residueType);
          //cout<<"51mass + projectileMass: "<<mass + projectileMass<<" "<<residueA<<endl;
      neutronsN = 1;
      productMass.push_back(0.0);
      GetSecParameter(targetElement, rp);
      }break;
      case 52:{
      residueA = mass + projectileMass - 1;
      residueZ = charge + projectileCharge;
      productsName.push_back("n and residue at 2nd excited state");
      residueType = residueName(residueZ, residueA);
      productsName.push_back(residueType);
      neutronsN = 1;
      productMass.push_back(0.0);
      GetSecParameter(targetElement, rp);
      }break;
      case 53:{
      residueA = mass + projectileMass - 1;
      residueZ = charge + projectileCharge;
      productsName.push_back("n and residue at 3rd excited state");
      residueType = residueName(residueZ, residueA);
      productsName.push_back(residueType);
      neutronsN = 1;
      productMass.push_back(0.0);
      GetSecParameter(targetElement, rp);
      }break;
      case 54:{
      residueA = mass + projectileMass - 1;
      residueZ = charge + projectileCharge;
      productsName.push_back("n and residue at 4th excited state");
      residueType = residueName(residueZ, residueA);
      productsName.push_back(residueType);
      neutronsN = 1;
      productMass.push_back(0.0);
      GetSecParameter(targetElement, rp);
      }break;
      case 55:{
      residueA = mass + projectileMass - 1;
      residueZ = charge + projectileCharge;
      productsName.push_back("n and residue at 5th excited state");
      residueType = residueName(residueZ, residueA);
      productsName.push_back(residueType);
      neutronsN = 1;
      productMass.push_back(0.0);
      GetSecParameter(targetElement, rp);
      }break;
      case 56:{
      residueA = mass + projectileMass - 1;
      residueZ = charge + projectileCharge;
      productsName.push_back("n and residue at 6th excited state");
      residueType = residueName(residueZ, residueA);
      productsName.push_back(residueType);
      neutronsN = 1;
      productMass.push_back(0.0);
      GetSecParameter(targetElement, rp);
      }break;
      case 57:{
      residueA = mass + projectileMass - 1;
      residueZ = charge + projectileCharge;
      productsName.push_back("n and residue at 7th excited state");
      residueType = residueName(residueZ, residueA);
      productsName.push_back(residueType);
      neutronsN = 1;
      productMass.push_back(0.0);
      GetSecParameter(targetElement, rp);
      break;}
      case 58:{
      residueA = mass + projectileMass - 1;
      residueZ = charge + projectileCharge;
      productsName.push_back("n and residue at 8th excited state");
      residueType = residueName(residueZ, residueA);
      productsName.push_back(residueType);
      neutronsN = 1;
      productMass.push_back(0.0);
      GetSecParameter(targetElement, rp);
      }break;
      case 59:{
      residueA = mass + projectileMass - 1;
      residueZ = charge + projectileCharge;
      productsName.push_back("n and residue at 9th excited state");
      residueType = residueName(residueZ, residueA);
      productsName.push_back(residueType);
      neutronsN = 1;
      productMass.push_back(0.0);
      GetSecParameter(targetElement, rp);
      }break;
      case 60:{
      residueA = mass + projectileMass - 1;
      residueZ = charge + projectileCharge;
      productsName.push_back("n and residue at 10th excited state");
      residueType = residueName(residueZ, residueA);
      productsName.push_back(residueType);
      neutronsN = 1;
      productMass.push_back(0.0);
      GetSecParameter(targetElement, rp);
      }break;
      case 61:{
      residueA = mass + projectileMass - 1;
      residueZ = charge + projectileCharge;
      productsName.push_back("n and residue at 11th excited state");
      residueType = residueName(residueZ, residueA);
      productsName.push_back(residueType);
      neutronsN = 1;
      productMass.push_back(0.0);
      GetSecParameter(targetElement, rp);
      }break;
      case 62:{
      residueA = mass + projectileMass - 1;
      residueZ = charge + projectileCharge;
      productsName.push_back("n and residue at 12th excited state");
      residueType = residueName(residueZ, residueA);
      productsName.push_back(residueType);
      neutronsN = 1;
      productMass.push_back(0.0);
      GetSecParameter(targetElement, rp);
      }break;
      case 63:{
      residueA = mass + projectileMass - 1;
      residueZ = charge + projectileCharge;
      productsName.push_back("n and residue at 13th excited state");
      residueType = residueName(residueZ, residueA);
      productsName.push_back(residueType);
      neutronsN = 1;
      productMass.push_back(0.0);
      GetSecParameter(targetElement, rp);
      }break;
      case 64:{
      residueA = mass + projectileMass - 1;
      residueZ = charge + projectileCharge;
      productsName.push_back("n and residue at 14th excited state");
      residueType = residueName(residueZ, residueA);
      productsName.push_back(residueType);
      neutronsN = 1;
      productMass.push_back(0.0);
      GetSecParameter(targetElement, rp);
      }break;
      case 65:{
      residueA = mass + projectileMass - 1;
      residueZ = charge + projectileCharge;
      productsName.push_back("n and residue at 15th excited state");
      residueType = residueName(residueZ, residueA);
      productsName.push_back(residueType);
      neutronsN = 1;
      productMass.push_back(0.0);
      GetSecParameter(targetElement, rp);
      }break;
      case 66:{
      residueA = mass + projectileMass - 1;
      residueZ = charge + projectileCharge;
      productsName.push_back("n and residue at 16th excited state");
      residueType = residueName(residueZ, residueA);
      productsName.push_back(residueType);
      neutronsN = 1;
      productMass.push_back(0.0);
      GetSecParameter(targetElement, rp);
      }break;
      case 67:{
      residueA = mass + projectileMass - 1;
      residueZ = charge + projectileCharge;
      productsName.push_back("n and residue at 17th excited state");
      residueType = residueName(residueZ, residueA);
      productsName.push_back(residueType);
      neutronsN = 1;
      productMass.push_back(0.0);
      GetSecParameter(targetElement, rp);
      }break;
      case 68:{
      residueA = mass + projectileMass - 1;
      residueZ = charge + projectileCharge;
      productsName.push_back("n and residue at 18th excited state");
      residueType = residueName(residueZ, residueA);
      productsName.push_back(residueType);
      neutronsN = 1;
      productMass.push_back(0.0);
      GetSecParameter(targetElement, rp);
      }break;
      case 69:{
      residueA = mass + projectileMass - 1;
      residueZ = charge + projectileCharge;
      productsName.push_back("n and residue at 19th excited state");
      residueType = residueName(residueZ, residueA);
      productsName.push_back(residueType);
      neutronsN = 1;
      productMass.push_back(0.0);
      GetSecParameter(targetElement, rp);
      }break;
      case 70:{
      residueA = mass + projectileMass - 1;
      residueZ = charge + projectileCharge;
      productsName.push_back("n and residue at 20th excited state");
      residueType = residueName(residueZ, residueA);
      productsName.push_back(residueType);
      neutronsN = 1;
      productMass.push_back(0.0);
      GetSecParameter(targetElement, rp);
      }break;
      case 71:{
      residueA = mass + projectileMass - 1;
      residueZ = charge + projectileCharge;
      productsName.push_back("n and residue at 21st excited state");
      residueType = residueName(residueZ, residueA);
      productsName.push_back(residueType);
      neutronsN = 1;
      productMass.push_back(0.0);
      GetSecParameter(targetElement, rp);
      }break;
      case 72:{
      residueA = mass + projectileMass - 1;
      residueZ = charge + projectileCharge;
      productsName.push_back("n and residue at 22nd excited state");
      residueType = residueName(residueZ, residueA);
      productsName.push_back(residueType);
      neutronsN = 1;
      productMass.push_back(0.0);
      GetSecParameter(targetElement, rp);
      }break;
      case 73:{
      residueA = mass + projectileMass - 1;
      residueZ = charge + projectileCharge;
      productsName.push_back("n and residue at 23rd excited state");
      residueType = residueName(residueZ, residueA);
      productsName.push_back(residueType);
      neutronsN = 1;
      productMass.push_back(0.0);
      GetSecParameter(targetElement, rp);
      }break;
      case 74:{
      residueA = mass + projectileMass - 1;
      residueZ = charge + projectileCharge;
      productsName.push_back("n and residue at 24th excited state");
      residueType = residueName(residueZ, residueA);
      productsName.push_back(residueType);
      neutronsN = 1;
      productMass.push_back(0.0);
      GetSecParameter(targetElement, rp);
      }break;
      case 75:{
      residueA = mass + projectileMass - 1;
      residueZ = charge + projectileCharge;
      productsName.push_back("n and residue at 25th excited state");
      residueType = residueName(residueZ, residueA);
      productsName.push_back(residueType);
      neutronsN = 1;
      productMass.push_back(0.0);
      GetSecParameter(targetElement, rp);
      }break;
      case 76:{
      residueA = mass + projectileMass - 1;
      residueZ = charge + projectileCharge;
      productsName.push_back("n and residue at 26th excited state");
      residueType = residueName(residueZ, residueA);
      productsName.push_back(residueType);
      neutronsN = 1;
      productMass.push_back(0.0);
      GetSecParameter(targetElement, rp);
      }break;
      case 77:{
      residueA = mass + projectileMass - 1;
      residueZ = charge + projectileCharge;
      productsName.push_back("n and residue at 27th excited state");
      residueType = residueName(residueZ, residueA);
      productsName.push_back(residueType);
      neutronsN = 1;
      productMass.push_back(0.0);
      GetSecParameter(targetElement, rp);
      }break;
      case 78:{
      residueA = mass + projectileMass - 1;
      residueZ = charge + projectileCharge;
      productsName.push_back("n and residue at 28th excited state");
      residueType = residueName(residueZ, residueA);
      productsName.push_back(residueType);
      neutronsN = 1;
      productMass.push_back(0.0);
      GetSecParameter(targetElement, rp);
      }break;
      case 79:{
      residueA = mass + projectileMass - 1;
      residueZ = charge + projectileCharge;
      productsName.push_back("n and residue at 29th excited state");
      residueType = residueName(residueZ, residueA);
      productsName.push_back(residueType);
      neutronsN = 1;
      productMass.push_back(0.0);
      GetSecParameter(targetElement, rp);
      }break;
      case 80:{
      residueA = mass + projectileMass - 1;
      residueZ = charge + projectileCharge;
      productsName.push_back("n and residue at 30th excited state");
      residueType = residueName(residueZ, residueA);
      productsName.push_back(residueType);
      neutronsN = 1;
      productMass.push_back(0.0);
      GetSecParameter(targetElement, rp);
      }break;
      case 81:{
      residueA = mass + projectileMass - 1;
      residueZ = charge + projectileCharge;
      productsName.push_back("n and residue at 31st excited state");
      residueType = residueName(residueZ, residueA);
      productsName.push_back(residueType);
      neutronsN = 1;
      productMass.push_back(0.0);
      GetSecParameter(targetElement, rp);
      }break;
      case 82:{
      residueA = mass + projectileMass - 1;
      residueZ = charge + projectileCharge;
      productsName.push_back("n and residue at 32nd excited state");
      residueType = residueName(residueZ, residueA);
      productsName.push_back(residueType);
      neutronsN = 1;
      productMass.push_back(0.0);
      GetSecParameter(targetElement, rp);
      }break;
      case 83:{
      residueA = mass + projectileMass - 1;
      residueZ = charge + projectileCharge;
      productsName.push_back("n and residue at 33rd excited state");
      residueType = residueName(residueZ, residueA);
      productsName.push_back(residueType);
      neutronsN = 1;
      productMass.push_back(0.0);
      GetSecParameter(targetElement, rp);
      }break;
      case 84:{
      residueA = mass + projectileMass - 1;
      residueZ = charge + projectileCharge;
      productsName.push_back("n and residue at 34th excited state");
      residueType = residueName(residueZ, residueA);
      productsName.push_back(residueType);
      neutronsN = 1;
      productMass.push_back(0.0);
      GetSecParameter(targetElement, rp);
      }break;
      case 85:{
      residueA = mass + projectileMass - 1;
      residueZ = charge + projectileCharge;
      productsName.push_back("n and residue at 35th excited state");
      residueType = residueName(residueZ, residueA);
      productsName.push_back(residueType);
      neutronsN = 1;
      productMass.push_back(0.0);
      GetSecParameter(targetElement, rp);
      }break;
      case 86:{
      residueA = mass + projectileMass - 1;
      residueZ = charge + projectileCharge;
      productsName.push_back("n and residue at 36th excited state");
      residueType = residueName(residueZ, residueA);
      productsName.push_back(residueType);
      neutronsN = 1;
      productMass.push_back(0.0);
      GetSecParameter(targetElement, rp);
      }break;
      case 87:{
      residueA = mass + projectileMass - 1;
      residueZ = charge + projectileCharge;
      productsName.push_back("n and residue at 37th excited state");
      residueType = residueName(residueZ, residueA);
      productsName.push_back(residueType);
      neutronsN = 1;
      productMass.push_back(0.0);
      GetSecParameter(targetElement, rp);
      }break;
      case 88:{
      residueA = mass + projectileMass - 1;
      residueZ = charge + projectileCharge;
      productsName.push_back("n and residue at 38th excited state");
      residueType = residueName(residueZ, residueA);
      productsName.push_back(residueType);
      neutronsN = 1;
      productMass.push_back(0.0);
      GetSecParameter(targetElement, rp);
      }break;
      case 89:{
      residueA = mass + projectileMass - 1;
      residueZ = charge + projectileCharge;
      productsName.push_back("n and residue at 39th excited state");
      residueType = residueName(residueZ, residueA);
      productsName.push_back(residueType);
      neutronsN = 1;
      productMass.push_back(0.0);
      GetSecParameter(targetElement, rp);
      }break;
      case 90:{
      residueA = mass + projectileMass - 1;
      residueZ = charge + projectileCharge;
      productsName.push_back("n and residue at 40th excited state");
      residueType = residueName(residueZ, residueA);
      productsName.push_back(residueType);
      neutronsN = 1;
      productMass.push_back(0.0);
      GetSecParameter(targetElement, rp);
      }break;
      case 91:{
      residueA = mass + projectileMass - 1;
      residueZ = charge + projectileCharge;
      productsName.push_back("n_continum"); //production of neutron from the continuum of the excited nucleus
      residueType = residueName(residueZ, residueA);
      productsName.push_back(residueType);
      neutronsN = 1;
      productMass.push_back(0.0);
      GetSecParameter(targetElement, rp);
      }break;
      //case 102:{
      //}break;
      case 103: // p
      residueA = mass + projectileMass - 1;
      residueZ = charge + projectileCharge - 1;
      productsName.push_back("p");
      residueType = residueName(residueZ, residueA);
      productsName.push_back(residueType);
      GetSecParameter(targetElement, rp);
      neutronsN = 0;
      productMass.push_back(1.0);
      break;
      case 104: // d
      residueA = mass + projectileMass - 2;
      residueZ = charge + projectileCharge - 1;
      productsName.push_back("1H2");
      residueType = residueName(residueZ, residueA);
      productsName.push_back(residueType);
      neutronsN = 0;
      productMass.push_back(2.0);
      GetSecParameter(targetElement, rp);
      break;
      case 105: // t
      residueA = mass + projectileMass - 3;
      residueZ = charge + projectileCharge - 1;
      productsName.push_back("1H3");
      residueType = residueName(residueZ, residueA);
      productsName.push_back(residueType);
      neutronsN = 0;
      productMass.push_back(3.0);
      GetSecParameter(targetElement, rp);
      //cout<<"Residue type: "<<residueType<<endl;
      break;
      case 106: // He3
      residueA = mass + projectileMass - 3;
      residueZ = charge + projectileCharge - 2;
      productsName.push_back("2He3");
      residueType = residueName(residueZ, residueA);
      productsName.push_back(residueType);
      neutronsN = 0;
      productMass.push_back(3.0);
      GetSecParameter(targetElement, rp);
      break;
      case 107: // alpha
      residueA = mass + projectileMass - 4;
      residueZ = charge + projectileCharge - 2;
      productsName.push_back("2He4");
      //cout<<"MT---------107:\t"<<cosLab<<"\t"<<secEnergyLab<<endl;
      residueType = residueName(residueZ, residueA);
      productsName.push_back(residueType);
      neutronsN = 0;
      productMass.push_back(4.0);
      GetSecParameter(targetElement, rp);
      break;
      case 108: // 2a
      residueA = mass + projectileMass - 8;
      residueZ = charge + projectileCharge - 4;
      productsName.push_back("2He4");
      productsName.push_back("2He4");
      residueType = residueName(residueZ, residueA);
      productsName.push_back(residueType);
      neutronsN = 0;
      productMass.push_back(4.0);
      productMass.push_back(4.0);
      GetSecParameter(targetElement, rp);
      break;
      case 109: // 3a
      residueA = mass + projectileMass - 12;
      residueZ = charge + projectileCharge - 6;
      productsName.push_back("2He4");
      productsName.push_back("2He4");
      productsName.push_back("2He4");
      residueType = residueName(residueZ, residueA);
      productsName.push_back(residueType);
      neutronsN = 0;
      productMass.push_back(4.0);
      productMass.push_back(4.0);
      productMass.push_back(4.0);
      GetSecParameter(targetElement, rp);
      break;
      case 111: // 2p
      residueA = mass + projectileMass - 2;
      residueZ = charge + projectileCharge - 2;
      if((MF4 == 4 && MF5 == 5) || (MF4 == 4 && MF5 == -1) || (MF4 == 99 && MF5 == 5)){
      tempproductsName.push_back("p");
      tempproductsName.push_back("p");
      residueType = residueName(residueZ, residueA);
      tempproductsName.push_back(residueType);
      }else if(MF4 == 99 && MF5 == -1 && MF6 == 6){
      tempproductsName.push_back("p");
      residueType = residueName(residueZ, residueA);
      tempproductsName.push_back(residueType);
      tempproductsName.push_back("p");
      }
      productsName = tempproductsName;
      neutronsN = 0;
      productMass.push_back(1.0);
      productMass.push_back(1.0);
      GetSecParameter(targetElement, rp);
      break;
      case 112: // p+a
      residueA = mass + projectileMass - 5;
      residueZ = charge + projectileCharge - 3;
      productsName.push_back("p");
      productsName.push_back("2He4");
      residueType = residueName(residueZ, residueA);
      productsName.push_back(residueType);
      neutronsN = 0;
      productMass.push_back(1.0);
      productMass.push_back(4.0);
      GetSecParameter(targetElement, rp);
      break;
      case 113: // t+2a
      residueA = mass + projectileMass - 11;
      residueZ = charge + projectileCharge - 5;
      productsName.push_back("3H");
      productsName.push_back("2He4");
      productsName.push_back("2He4");
      residueType = residueName(residueZ, residueA);
      productsName.push_back(residueType);
      neutronsN = 0;
      productMass.push_back(3.0);
      productMass.push_back(4.0);
      productMass.push_back(4.0);
      GetSecParameter(targetElement, rp);
      break;
      case 114: // d+2a
      residueA = mass + projectileMass - 10;
      residueZ = charge + projectileCharge - 5;
      productsName.push_back("1H2");
      productsName.push_back("2He4");
      productsName.push_back("2He4");
      residueType = residueName(residueZ, residueA);
      productsName.push_back(residueType);
      neutronsN = 0;
      productMass.push_back(2.0);
      productMass.push_back(4.0);
      productMass.push_back(4.0);
      GetSecParameter(targetElement, rp);
      break;
      case 115: // p+d
      residueA = mass + projectileMass - 3;
      residueZ = charge + projectileCharge - 2;
      productsName.push_back("p");
      productsName.push_back("1H2");
      residueType = residueName(residueZ, residueA);
      productsName.push_back(residueType);
      neutronsN = 0;
      productMass.push_back(1.0);
      productMass.push_back(2.0);
      GetSecParameter(targetElement, rp);
      break;
      case 116: // p+t
      residueA = mass + projectileMass - 4;
      residueZ = charge + projectileCharge - 2;
      productsName.push_back("p");
      productsName.push_back("1H3");
      residueType = residueName(residueZ, residueA);
      productsName.push_back(residueType);
      neutronsN = 0;
      productMass.push_back(1.0);
      productMass.push_back(3.0);
      GetSecParameter(targetElement, rp);
      break;
      case 117: // d+a
      residueA = mass + projectileMass - 6;
      residueZ = charge + projectileCharge - 3;
      productsName.push_back("1H2");
      productsName.push_back("2He4");
      residueType = residueName(residueZ, residueA);
      productsName.push_back(residueType);
      neutronsN = 0;
      productMass.push_back(2.0);
      productMass.push_back(4.0);
      GetSecParameter(targetElement, rp);
      break;
      case 152: // 5n
      residueA = mass + projectileMass - 5;
      residueZ = charge + projectileCharge;
      productsName.push_back("n");
      productsName.push_back("n");
      productsName.push_back("n");
      productsName.push_back("n");
      productsName.push_back("n");
      residueType = residueName(residueZ, residueA);
      productsName.push_back(residueType);
      neutronsN = 5;
      productMass.push_back(0.0);
      GetSecParameter(targetElement, rp);
      break;
      case 153: // 6n
      residueA = mass + projectileMass - 6;
      residueZ = charge + projectileCharge;
      productsName.push_back("n");
      productsName.push_back("n");
      productsName.push_back("n");
      productsName.push_back("n");
      productsName.push_back("n");
      productsName.push_back("n");
      residueType = residueName(residueZ, residueA);
      productsName.push_back(residueType);
      neutronsN = 6;
      productMass.push_back(0.0);
      GetSecParameter(targetElement, rp);
      break;
      case 154: // 2n+t
      residueA = mass + projectileMass - 5;
      residueZ = charge + projectileCharge - 1;
      productsName.push_back("n");
      productsName.push_back("n");
      productsName.push_back("1H3");
      residueType = residueName(residueZ, residueA);
      productsName.push_back(residueType);
      neutronsN = 2;
      productMass.push_back(3.0);
      GetSecParameter(targetElement, rp);
      break;
      case 155: // t+a
      residueA = mass + projectileMass - 7;
      residueZ = charge + projectileCharge - 3;
      productsName.push_back("1H3");
      productsName.push_back("2He4");
      residueType = residueName(residueZ, residueA);
      productsName.push_back(residueType);
      neutronsN = 0;
      productMass.push_back(3.0);
      productMass.push_back(4.0);
      GetSecParameter(targetElement, rp);
      break;
      case 156: // 4n+p
      residueA = mass + projectileMass - 5;
      residueZ = charge + projectileCharge - 1;
      productsName.push_back("n");
      productsName.push_back("n");
      productsName.push_back("n");
      productsName.push_back("n");
      productsName.push_back("p");
      residueType = residueName(residueZ, residueA);
      productsName.push_back(residueType);
      neutronsN = 4;
      productMass.push_back(1.0);
      GetSecParameter(targetElement, rp);
      break;
      case 157: // 3n+d
      residueA = mass + projectileMass - 5;
      residueZ = charge + projectileCharge - 1;
      productsName.push_back("n");
      productsName.push_back("n");
      productsName.push_back("n");
      productsName.push_back("1H2");
      residueType = residueName(residueZ, residueA);
      productsName.push_back(residueType);
      neutronsN = 3;
      productMass.push_back(2.0);
      GetSecParameter(targetElement, rp);
      break;
      case 158: // n+d+a
      residueA = mass + projectileMass - 7;
      residueZ = charge + projectileCharge - 3;
      productsName.push_back("n");
      productsName.push_back("1H2");
      productsName.push_back("2He4");
      residueType = residueName(residueZ, residueA);
      productsName.push_back(residueType);
      neutronsN = 1;
      productMass.push_back(2.0);
      productMass.push_back(4.0);
      GetSecParameter(targetElement, rp);
      break;
      case 159: // 2n+p+a
      residueA = mass + projectileMass - 7;
      residueZ = charge + projectileCharge - 3;
      productsName.push_back("n");
      productsName.push_back("n");
      productsName.push_back("p");
      productsName.push_back("2He4");
      residueType = residueName(residueZ, residueA);
      productsName.push_back(residueType);
      neutronsN = 2;
      productMass.push_back(1.0);
      productMass.push_back(4.0);
      GetSecParameter(targetElement, rp);
      break;
      case 160: // 7n
      residueA = mass + projectileMass - 7;
      residueZ = charge + projectileCharge;
      productsName.push_back("n");
      productsName.push_back("n");
      productsName.push_back("n");
      productsName.push_back("n");
      productsName.push_back("n");
      productsName.push_back("n");
      productsName.push_back("n");
      residueType = residueName(residueZ, residueA);
      productsName.push_back(residueType);
      neutronsN = 7;
      productMass.push_back(0.0);
      GetSecParameter(targetElement, rp);
      break;
      case 161: // 8n
      residueA = mass + projectileMass - 8;
      residueZ = charge + projectileCharge;
      productsName.push_back("n");
      productsName.push_back("n");
      productsName.push_back("n");
      productsName.push_back("n");
      productsName.push_back("n");
      productsName.push_back("n");
      productsName.push_back("n");
      productsName.push_back("n");
      residueType = residueName(residueZ, residueA);
      productsName.push_back(residueType);
      neutronsN = 8;
      productMass.push_back(0.0);
      GetSecParameter(targetElement, rp);
      break;
      case 162: // 5np
      residueA = mass + projectileMass - 6;
      residueZ = charge + projectileCharge - 1;
      productsName.push_back("n");
      productsName.push_back("n");
      productsName.push_back("n");
      productsName.push_back("n");
      productsName.push_back("n");
      productsName.push_back("p");
      residueType = residueName(residueZ, residueA);
      productsName.push_back(residueType);
      neutronsN = 5;
      productMass.push_back(1.0);
      GetSecParameter(targetElement, rp);
      break;
      case 163: // 6np
      residueA = mass + projectileMass - 7;
      residueZ = charge + projectileCharge - 1;
      productsName.push_back("n");
      productsName.push_back("n");
      productsName.push_back("n");
      productsName.push_back("n");
      productsName.push_back("n");
      productsName.push_back("n");
      productsName.push_back("p");
      residueType = residueName(residueZ, residueA);
      productsName.push_back(residueType);
      neutronsN = 6;
      productMass.push_back(1.0);
      GetSecParameter(targetElement, rp);
      break;
      case 164: // 7np
      residueA = mass + projectileMass - 8;
      residueZ = charge + projectileCharge - 1;
      productsName.push_back("n");
      productsName.push_back("n");
      productsName.push_back("n");
      productsName.push_back("n");
      productsName.push_back("n");
      productsName.push_back("n");
      productsName.push_back("n");
      productsName.push_back("p");
      residueType = residueName(residueZ, residueA);
      productsName.push_back(residueType);
      neutronsN = 7;
      productMass.push_back(1.0);
      GetSecParameter(targetElement, rp);
      break;
      case 165: // 4n+a
      residueA = mass + projectileMass - 8;
      residueZ = charge + projectileCharge - 2;
      productsName.push_back("n");
      productsName.push_back("n");
      productsName.push_back("n");
      productsName.push_back("n");
      productsName.push_back("2He4");
      residueType = residueName(residueZ, residueA);
      productsName.push_back(residueType);
      neutronsN = 4;
      productMass.push_back(4.0);
      GetSecParameter(targetElement, rp);
      break;
      case 166: // 5na
      residueA = mass + projectileMass - 9;
      residueZ = charge + projectileCharge - 2;
      productsName.push_back("n");
      productsName.push_back("n");
      productsName.push_back("n");
      productsName.push_back("n");
      productsName.push_back("n");
      productsName.push_back("2He4");
      residueType = residueName(residueZ, residueA);
      productsName.push_back(residueType);
      neutronsN = 5;
      productMass.push_back(4.0);
      GetSecParameter(targetElement, rp);
      break;
      case 167: // 6na
      residueA = mass + projectileMass - 10;
      residueZ = charge + projectileCharge - 2;
      productsName.push_back("n");
      productsName.push_back("n");
      productsName.push_back("n");
      productsName.push_back("n");
      productsName.push_back("n");
      productsName.push_back("n");
      productsName.push_back("2He4");
      residueType = residueName(residueZ, residueA);
      productsName.push_back(residueType);
      neutronsN = 6;
      productMass.push_back(4.0);
      GetSecParameter(targetElement, rp);
      break;
      case 168: // 7na
      residueA = mass + projectileMass - 11;
      residueZ = charge + projectileCharge - 2;
      productsName.push_back("n");
      productsName.push_back("n");
      productsName.push_back("n");
      productsName.push_back("n");
      productsName.push_back("n");
      productsName.push_back("n");
      productsName.push_back("n");
      productsName.push_back("2He4");
      residueType = residueName(residueZ, residueA);
      productsName.push_back(residueType);
      neutronsN = 7;
      productMass.push_back(4.0);
      GetSecParameter(targetElement, rp);
      break;
      case 169: // 4nd
      residueA = mass + projectileMass - 6;
      residueZ = charge + projectileCharge - 1;
      productsName.push_back("n");
      productsName.push_back("n");
      productsName.push_back("n");
      productsName.push_back("n");
      productsName.push_back("1H2");
      residueType = residueName(residueZ, residueA);
      productsName.push_back(residueType);
      neutronsN = 4;
      productMass.push_back(2.0);
      GetSecParameter(targetElement, rp);
      break;
      case 170: // 5nd
      residueA = mass + projectileMass - 7;
      residueZ = charge + projectileCharge - 1;
      productsName.push_back("n");
      productsName.push_back("n");
      productsName.push_back("n");
      productsName.push_back("n");
      productsName.push_back("n");
      productsName.push_back("1H2");
      residueType = residueName(residueZ, residueA);
      productsName.push_back(residueType);
      neutronsN = 5;
      productMass.push_back(2.0);
      GetSecParameter(targetElement, rp);
      break;
      case 171: // 6nd
      residueA = mass + projectileMass - 8;
      residueZ = charge + projectileCharge - 1;
      productsName.push_back("n");
      productsName.push_back("n");
      productsName.push_back("n");
      productsName.push_back("n");
      productsName.push_back("n");
      productsName.push_back("n");
      productsName.push_back("1H2");
      residueType = residueName(residueZ, residueA);
      productsName.push_back(residueType);
      neutronsN = 6;
      productMass.push_back(2.0);
      GetSecParameter(targetElement, rp);
      break;
      case 172: // 3nt
      residueA = mass + projectileMass - 6;
      residueZ = charge + projectileCharge - 1;
      productsName.push_back("n");
      productsName.push_back("n");
      productsName.push_back("n");
      productsName.push_back("1H3");
      residueType = residueName(residueZ, residueA);
      productsName.push_back(residueType);
      neutronsN = 3;
      productMass.push_back(3.0);
      GetSecParameter(targetElement, rp);
      break;
      case 173: // 4nt
      residueA = mass + projectileMass - 7;
      residueZ = charge + projectileCharge - 1;
      productsName.push_back("n");
      productsName.push_back("n");
      productsName.push_back("n");
      productsName.push_back("n");
      productsName.push_back("1H3");
      residueType = residueName(residueZ, residueA);
      productsName.push_back(residueType);
      neutronsN = 4;
      productMass.push_back(3.0);
      GetSecParameter(targetElement, rp);
      break;
      case 174: // 5nt
      residueA = mass + projectileMass - 8;
      residueZ = charge + projectileCharge - 1;
      productsName.push_back("n");
      productsName.push_back("n");
      productsName.push_back("n");
      productsName.push_back("n");
      productsName.push_back("n");
      productsName.push_back("1H3");
      residueType = residueName(residueZ, residueA);
      productsName.push_back(residueType);
      neutronsN = 5;
      productMass.push_back(3.0);
      GetSecParameter(targetElement, rp);
      break;
      case 175: // 6nt
      residueA = mass + projectileMass - 9;
      residueZ = charge + projectileCharge - 1;
      productsName.push_back("n");
      productsName.push_back("n");
      productsName.push_back("n");
      productsName.push_back("n");
      productsName.push_back("n");
      productsName.push_back("n");
      productsName.push_back("1H3");
      residueType = residueName(residueZ, residueA);
      productsName.push_back(residueType);
      neutronsN = 6;
      productMass.push_back(3.0);
      GetSecParameter(targetElement, rp);
      break;
      case 176: // 2n+He3
      residueA = mass + projectileMass - 5;
      residueZ = charge + projectileCharge - 2;
      productsName.push_back("n");
      productsName.push_back("n");
      productsName.push_back("2He3");
      residueType = residueName(residueZ, residueA);
      productsName.push_back(residueType);
      neutronsN = 2;
      productMass.push_back(3.0);
      GetSecParameter(targetElement, rp);
      break;
      case 177: // 3n + He3
      residueA = mass + projectileMass - 6;
      residueZ = charge + projectileCharge - 2;
      productsName.push_back("n");
      productsName.push_back("n");
      productsName.push_back("n");
      productsName.push_back("2He3");
      residueType = residueName(residueZ, residueA);
      productsName.push_back(residueType);
      neutronsN = 3;
      productMass.push_back(3.0);
      GetSecParameter(targetElement, rp);
      break;
      case 178: // 4n +He3
      residueA = mass + projectileMass - 7;
      residueZ = charge + projectileCharge - 2;
      productsName.push_back("n");
      productsName.push_back("n");
      productsName.push_back("n");
      productsName.push_back("n");
      productsName.push_back("2He3");
      residueType = residueName(residueZ, residueA);
      productsName.push_back(residueType);
      neutronsN = 4;
      productMass.push_back(3.0);
      GetSecParameter(targetElement, rp);
      break;
      case 179: // 3n2p
      residueA = mass + projectileMass - 5;
      residueZ = charge + projectileCharge - 2;
      productsName.push_back("n");
      productsName.push_back("n");
      productsName.push_back("n");
      productsName.push_back("p");
      productsName.push_back("p");
      residueType = residueName(residueZ, residueA);
      productsName.push_back(residueType);
      neutronsN = 3;
      productMass.push_back(1.0);
      productMass.push_back(1.0);
      GetSecParameter(targetElement, rp);
      break;
      case 180: // 3n2a
      residueA = mass + projectileMass - 11;
      residueZ = charge + projectileCharge - 4;
      productsName.push_back("n");
      productsName.push_back("n");
      productsName.push_back("n");
      productsName.push_back("2He4");
      productsName.push_back("2He4");
      residueType = residueName(residueZ, residueA);
      productsName.push_back(residueType);
      neutronsN = 3;
      productMass.push_back(4.0);
      productMass.push_back(4.0);
      GetSecParameter(targetElement, rp);
      break;
      case 181: // 3npa
      residueA = mass + projectileMass - 8;
      residueZ = charge + projectileCharge - 3;
      productsName.push_back("n");
      productsName.push_back("n");
      productsName.push_back("n");
      productsName.push_back("p");
      productsName.push_back("2He4");
      residueType = residueName(residueZ, residueA);
      productsName.push_back(residueType);
      neutronsN = 3;
      productMass.push_back(1.0);
      productMass.push_back(4.0);
      GetSecParameter(targetElement, rp);
      break;
      case 182: // dt
      residueA = mass + projectileMass - 5;
      residueZ = charge + projectileCharge - 2;
      productsName.push_back("1H2");
      productsName.push_back("1H3");
      residueType = residueName(residueZ, residueA);
      productsName.push_back(residueType);
      neutronsN = 0;
      productMass.push_back(2.0);
      productMass.push_back(3.0);
      GetSecParameter(targetElement, rp);
      break;
      case 183: // npd
      residueA = mass + projectileMass - 4;
      residueZ = charge + projectileCharge - 2;
      productsName.push_back("n");
      productsName.push_back("p");
      productsName.push_back("1H2");
      residueType = residueName(residueZ, residueA);
      productsName.push_back(residueType);
      neutronsN = 1;
      productMass.push_back(1.0);
      productMass.push_back(2.0);
      GetSecParameter(targetElement, rp);
      break;
      case 184: // npt
      residueA = mass + projectileMass - 5;
      residueZ = charge + projectileCharge - 2;
      productsName.push_back("n");
      productsName.push_back("p");
      productsName.push_back("1H3");
      residueType = residueName(residueZ, residueA);
      productsName.push_back(residueType);
      neutronsN = 1;
      productMass.push_back(1.0);
      productMass.push_back(3.0);
      GetSecParameter(targetElement, rp);
      break;
      case 185: // ndt
      residueA = mass + projectileMass - 6;
      residueZ = charge + projectileCharge - 2;
      productsName.push_back("n");
      productsName.push_back("1H2");
      productsName.push_back("1H3");
      residueType = residueName(residueZ, residueA);
      productsName.push_back(residueType);
      neutronsN = 1;
      productMass.push_back(2.0);
      productMass.push_back(3.0);
      GetSecParameter(targetElement, rp);
      break;
      case 186: // npHe3
      residueA = mass + projectileMass - 5;
      residueZ = charge + projectileCharge - 3;
      productsName.push_back("n");
      productsName.push_back("p");
      productsName.push_back("2He3");
      residueType = residueName(residueZ, residueA);
      productsName.push_back(residueType);
      neutronsN = 1;
      productMass.push_back(1.0);
      productMass.push_back(3.0);
      GetSecParameter(targetElement, rp);
      break;
      case 187: // ndHe3
      residueA = mass + projectileMass - 6;
      residueZ = charge + projectileCharge - 3;
      productsName.push_back("n");
      productsName.push_back("1H2");
      productsName.push_back("2He3");
      residueType = residueName(residueZ, residueA);
      productsName.push_back(residueType);
      neutronsN = 1;
      productMass.push_back(2.0);
      productMass.push_back(3.0);
      GetSecParameter(targetElement, rp);
      break;
      case 188: // ntHe3
      residueA = mass + projectileMass - 7;
      residueZ = charge + projectileCharge - 3;
      productsName.push_back("n");
      productsName.push_back("1H3");
      productsName.push_back("2He3");
      residueType = residueName(residueZ, residueA);
      productsName.push_back(residueType);
      neutronsN = 1;
      productMass.push_back(3.0);
      productMass.push_back(3.0);
      GetSecParameter(targetElement, rp);
      break;
      case 189: // nta
      residueA = mass + projectileMass - 8;
      residueZ = charge + projectileCharge - 3;
      productsName.push_back("n");
      productsName.push_back("1H3");
      productsName.push_back("2He4");
      residueType = residueName(residueZ, residueA);
      productsName.push_back(residueType);
      neutronsN = 1;
      productMass.push_back(3.0);
      productMass.push_back(4.0);
      GetSecParameter(targetElement, rp);
      break;
      case 190: // 2n2p
      residueA = mass + projectileMass - 4;
      residueZ = charge + projectileCharge - 2;
      productsName.push_back("n");
      productsName.push_back("n");
      productsName.push_back("p");
      productsName.push_back("p");
      residueType = residueName(residueZ, residueA);
      productsName.push_back(residueType);
      neutronsN = 2;
      productMass.push_back(1.0);
      productMass.push_back(1.0);
      GetSecParameter(targetElement, rp);
      break;
      case 191: // pHe3
      residueA = mass + projectileMass - 4;
      residueZ = charge + projectileCharge - 3;
      productsName.push_back("p");
      productsName.push_back("2He3");
      residueType = residueName(residueZ, residueA);
      productsName.push_back(residueType);
      neutronsN = 0;
      productMass.push_back(1.0);
      productMass.push_back(3.0);
      GetSecParameter(targetElement, rp);
      break;
      case 192: // dHe3
      residueA = mass + projectileMass - 5;
      residueZ = charge + projectileCharge - 3;
      productsName.push_back("1H2");
      productsName.push_back("2He3");
      residueType = residueName(residueZ, residueA);
      productsName.push_back(residueType);
      neutronsN = 0;
      productMass.push_back(2.0);
      productMass.push_back(3.0);
      GetSecParameter(targetElement, rp);
      break;
      case 193: // aHe3
      residueA = mass + projectileMass - 7;
      residueZ = charge + projectileCharge - 4;
      productsName.push_back("2He3");
      productsName.push_back("2He4");
      residueType = residueName(residueZ, residueA);
      productsName.push_back(residueType);
      neutronsN = 0;
      productMass.push_back(3.0);
      productMass.push_back(4.0);
      GetSecParameter(targetElement, rp);
      break;
      case 194: // 4n2p
      residueA = mass + projectileMass - 6;
      residueZ = charge + projectileCharge - 2;
      productsName.push_back("n");
      productsName.push_back("n");
      productsName.push_back("n");
      productsName.push_back("n");
      productsName.push_back("p");
      productsName.push_back("p");
      residueType = residueName(residueZ, residueA);
      productsName.push_back(residueType);
      neutronsN = 4;
      productMass.push_back(1.0);
      productMass.push_back(1.0);
      GetSecParameter(targetElement, rp);
      break;
      case 195: // 4n2a
      residueA = mass + projectileMass - 12;
      residueZ = charge + projectileCharge - 4;
      productsName.push_back("n");
      productsName.push_back("n");
      productsName.push_back("n");
      productsName.push_back("n");
      productsName.push_back("2He4");
      productsName.push_back("2He4");
      residueType = residueName(residueZ, residueA);
      productsName.push_back(residueType);
      neutronsN = 4;
      productMass.push_back(4.0);
      productMass.push_back(4.0);
      GetSecParameter(targetElement, rp);
      break;
      case 196: // 4npa
      residueA = mass + projectileMass - 9;
      residueZ = charge + projectileCharge - 3;
      productsName.push_back("n");
      productsName.push_back("n");
      productsName.push_back("n");
      productsName.push_back("n");
      productsName.push_back("p");
      productsName.push_back("2He4");
      residueType = residueName(residueZ, residueA);
      productsName.push_back(residueType);
      neutronsN = 4;
      productMass.push_back(1.0);
      productMass.push_back(4.0);
      GetSecParameter(targetElement, rp);
      break;
      case 197: // 3p
      residueA = mass + projectileMass - 3;
      residueZ = charge + projectileCharge - 3;
      productsName.push_back("p");
      productsName.push_back("p");
      productsName.push_back("p");
      residueType = residueName(residueZ, residueA);
      productsName.push_back(residueType);
      neutronsN = 0;
      productMass.push_back(1.0);
      productMass.push_back(1.0);
      productMass.push_back(1.0);
      GetSecParameter(targetElement, rp);
      break;
      case 198: // n3p
      residueA = mass + projectileMass - 4;
      residueZ = charge + projectileCharge - 3;
      productsName.push_back("n");
      productsName.push_back("p");
      productsName.push_back("p");
      productsName.push_back("p");
      residueType = residueName(residueZ, residueA);
      productsName.push_back(residueType);
      neutronsN = 1;
      productMass.push_back(1.0);
      productMass.push_back(1.0);
      productMass.push_back(1.0);
      GetSecParameter(targetElement, rp);
      break;
      case 199: // 3n2pa
      residueA = mass + projectileMass - 9;
      residueZ = charge + projectileCharge - 4;
      productsName.push_back("n");
      productsName.push_back("n");
      productsName.push_back("n");
      productsName.push_back("p");
      productsName.push_back("p");
      productsName.push_back("2He4");
      residueType = residueName(residueZ, residueA);
      productsName.push_back(residueType);
      neutronsN = 3;
      productMass.push_back(1.0);
      productMass.push_back(1.0);
      productMass.push_back(4.0);
      GetSecParameter(targetElement, rp);
      break;
      case 200: // 5n2p
      residueA = mass + projectileMass - 7;
      residueZ = charge + projectileCharge - 2;
      productsName.push_back("n");
      productsName.push_back("n");
      productsName.push_back("n");
      productsName.push_back("n");
      productsName.push_back("n");
      productsName.push_back("p");
      productsName.push_back("p");
      residueType = residueName(residueZ, residueA);
      productsName.push_back(residueType);
      neutronsN = 5;
      productMass.push_back(1.0);
      productMass.push_back(1.0);
      GetSecParameter(targetElement, rp);
      break;     
     }
      mtTemp.clear();
      crssIel.clear();
      crss.clear();
}
//______________________________________________________________________________
std::string TNudyInelastic::GetInelasticProcessName()
{
     return procNameInelastic;
}
//______________________________________________________________________________
double TNudyInelastic::GetInelasticXsec()
{
       return sigmaInelastic;
}
//______________________________________________________________________________
std::vector<std::string> TNudyInelastic::GetParticleInelastic()
{
       return productsName;
}
//______________________________________________________________________________
std::vector<double> TNudyInelastic::GetKiEnInelastic()
{
       return secEnergyinLab;
}
//______________________________________________________________________________
std::vector<double> TNudyInelastic::GetcosAngInelastic()
{
       return seccosAnginLab;
}
//______________________________________________________________________________
int TNudyInelastic::reactionChannelNumber()
{
       return MT;
}
//______________________________________________________________________________
void TNudyInelastic::GetInelasticParameters(double &xsec1, double &E)
{
  xsec1 = sigmaInelastic;
      E = En;
}
//______________________________________________________________________________
//LCT Flag to specify the frame of reference used
//LCT=1, the data are given in the LAB system
//LCT=2, the data are given in the CM system
//MF=File number
//MT=Reaction number due to various reaction type
void TNudyInelastic::GetSecParameter(TNudyElement *targetElement, TNudyEndfRecoPoint *recoPoint)
{
  double massT= targetElement->GetatomicMass();
  double projectileKineticEnergy;
  double nE[10];
  double nE2;
  double tempQValue = recoPoint->GetQValue(ielemId, MT);
  double netEnergy = 0.0 , sumEnergy;
  double EH,EHe;
  double product1Mass,lightMass;
  double product2Mass, heavyMass;
  double netMass;
  double lightParticleEnergy;
  double heavyParticleEnergy;
  double pi = acos(-1.0);
  double tempproductEnergy = 0;
  double tempresidueEnergy = 0;
  double tempcosAng1 = 0, tempcosAng2 = 0, tempcosAng3 = 0;
  nsecParticles = productsName.size(); 

  if (MF4 == 4 && MF5 == 5){
    if(productMass[0] > 0.0){
       check1:
       for(int imt = 0; imt < neutronsN; imt++){
          cosCM        = recoPoint->GetCos4(ielemId, MT, kineticE);
          secEnergyCM  = recoPoint->GetEnergy5(ielemId, MT, kineticE);
          secEnergyLab = TNudyCore::Instance()->cmToLabInelasticE(secEnergyCM, kineticE, cosCM, massT);
          nE[imt] = secEnergyLab;
          if (LCT == 2){
            cosLab = TNudyCore::Instance()->cmToLabInelasticCosT(secEnergyLab, secEnergyCM, kineticE, cosCM,
                                                           massT);
            }else if (LCT == 1) {
               cosLab       = cosCM;
               secEnergyLab = secEnergyCM;
            }
          seccosAnginLab.push_back(cosLab);
          secEnergyinLab.push_back(secEnergyLab);
         }
         sumEnergy=0;
         for(int jmt = 0; jmt < neutronsN; jmt++){
          sumEnergy = sumEnergy + nE[jmt];
         }
         if((kineticE + tempQValue - sumEnergy)<0.0) goto check1; 
         netEnergy = kineticE + tempQValue - sumEnergy;
         for(int jmt = 0; jmt < productMass.size(); jmt++){
           product1Mass = productMass[jmt];
           product2Mass = residueA;
           netMass = product1Mass + product2Mass;
          }
         lightParticleEnergy  = (product2Mass/netMass)*netEnergy; 
         heavyParticleEnergy = (product1Mass/netMass)*netEnergy;
         if( product1Mass < product2Mass){
            tempproductEnergy    = lightParticleEnergy;
            tempresidueEnergy = heavyParticleEnergy;
         }else if( product1Mass > product2Mass){
                tempproductEnergy    = lightParticleEnergy;
                tempresidueEnergy = heavyParticleEnergy;
           }else{
              tempproductEnergy = lightParticleEnergy;
              tempresidueEnergy = lightParticleEnergy;
            }
          secEnergyinLab.push_back(tempproductEnergy);
          secEnergyinLab.push_back(tempresidueEnergy);
          tempcosAng1 = 1 - 2*fRnd->Uniform(1);
          tempcosAng2 = cos(pi - acos(tempcosAng1));
          seccosAnginLab.push_back(tempcosAng1);
          seccosAnginLab.push_back(tempcosAng2); 
        //if(MT == 25 && tempresidueEnergy < 0.0 )cout<<kineticE + tempQValue<< "  sum energy:  "<<sumEnergy<<"    Remaining energy -neutron energy: "
          //             <<netEnergy<< " "<< "  Eh:  "<<tempproductEnergy<<" EHe: "<<tempresidueEnergy<<endl;
        }else if(productMass[0] == 0.0){
          neutronsN = neutronsN -1;
         for(int imt = 0; imt < neutronsN; imt++){
          cosCM        = recoPoint->GetCos4(ielemId, MT, kineticE);
          secEnergyCM  = recoPoint->GetEnergy5(ielemId, MT, kineticE);
          secEnergyLab = TNudyCore::Instance()->cmToLabInelasticE(secEnergyCM, kineticE, cosCM, massT);
          nE[imt] = secEnergyLab;
          if (LCT == 2){
            cosLab = TNudyCore::Instance()->cmToLabInelasticCosT(secEnergyLab, secEnergyCM, kineticE, cosCM,
                                                           massT);
            }else if (LCT == 1) {
               cosLab       = cosCM;
               secEnergyLab = secEnergyCM;
            }
          seccosAnginLab.push_back(cosLab);
          secEnergyinLab.push_back(secEnergyLab);
         }
         sumEnergy = 0;
         for(int jmt = 0; jmt < neutronsN; jmt++){
          sumEnergy = sumEnergy + nE[jmt];
         }
         netEnergy = kineticE + tempQValue - sumEnergy;
         for(int jmt = 0; jmt < productMass.size(); jmt++){
           product1Mass = 1;
           product2Mass = residueA;
           netMass = product1Mass + product2Mass;
         }
         lightParticleEnergy  = (product2Mass/netMass)*netEnergy; 
         heavyParticleEnergy = (product1Mass/netMass)*netEnergy;
         if( product1Mass < product2Mass){
            tempproductEnergy    = lightParticleEnergy;
            tempresidueEnergy = heavyParticleEnergy;
         }else if( product1Mass > product2Mass){
                tempproductEnergy    = heavyParticleEnergy;
                tempresidueEnergy = lightParticleEnergy;
           }else{
              tempproductEnergy = lightParticleEnergy;
              tempresidueEnergy = lightParticleEnergy;
            }
          secEnergyinLab.push_back(tempproductEnergy);
          secEnergyinLab.push_back(tempresidueEnergy);
        //cout<<kineticE + tempQValue<<"    Remaining energy -neutron energy: "
         //               <<netEnergy<<"  Eh:  "<<lightParticleEnergy<<" EHe: "<<tempresidueEnergy
         //               <<"  "<<lightParticleEnergy+heavyParticleEnergy<<"   "<<
         //               sumEnergy+lightParticleEnergy+heavyParticleEnergy<<endl;
          //residueEnergy.push_back(tempresidueEnergy); 
          //residueAngle.push_back(residueCosang); 
          tempcosAng1 = 1 - 2*fRnd->Uniform(1);
          tempcosAng2 = cos(pi - acos(tempcosAng1));
          seccosAnginLab.push_back(tempcosAng1);
          seccosAnginLab.push_back(tempcosAng2);
          secEnergyinLab.push_back(tempproductEnergy);
          secEnergyinLab.push_back(tempresidueEnergy);
        }  
    }else if (MF4 == 4 && MF5 == -1) {
    cosCM      = recoPoint->GetCos4(ielemId, MT, kineticE);
    double vCM,v3,v4,E3cm,E4cm, cosAngCM;// velocity and energy in CM
    double V3,V4, E3Lab,E4Lab,cosAngLab;// in Lab frame
    double Eci,cosAng3CM,sinAng3CM,ang3Lab,ang4Lab;
    double gamma, productCosangCM,productCosangLab;
    double fm1;  double fm2;  double fm3;  double fm4;  double fqval;double fEt;  
    double fm1m2;double fm3m4;double fA1;    
    double fB1;  double fC1;  double fD1;  double sinLab;
    double pi = acos(-1.0);
    double phi;
    double fact1,fact2,fact3;      
          fqval  = recoPoint->GetQValue(ielemId, MT); // in MeV // q-value
          fm1    = projectileMass;  // mass of projectile // neutron
          fm2    = massT; //mass of target// 6Li
          fm3    = massT + projectileMass - residueA; // mass of outgoing particle// 3H
          fm4    = residueA; // mass of daughter// 4He
          fm1m2  = fm1 + fm2;
          fm3m4  = fm3 + fm4;
          //cout<<fm1<<" "<<fm2<<" "<<fm3<<" "<<fm4<<std::endl;
          Eci    = (fm2/fm1m2)*kineticE;
          fEt    = Eci + fqval;  // final energy in CM
          vCM    = sqrt(2*fm1*kineticE)/fm1m2;//velocity of CM
          v3   = sqrt((2/fm3m4)*(fm4/fm3)*fEt); // velocity in CM-3
          v4   = - sqrt((2/fm3m4)*(fm3/fm4)*fEt); // velocity in CM-4
          E3cm = (fm4/fm3m4)*fEt; // energy in CM of Particle 3
          E4cm = (fm3/fm3m4)*fEt; // energy in CM of Particle 4
         if (LCT == 2) {
          productCosangCM =  cosCM;
          cosAng3CM  = productCosangCM;
          sinAng3CM = sqrt(1- productCosangCM* productCosangCM);  //cosine angle of particle 3 in CM
          gamma = vCM/v3; 
          ang3Lab = atan(sinAng3CM/(productCosangCM+ gamma)); 
          phi =  (pi - acos(cosAng3CM));
          V3   = vCM*cos(ang3Lab) + sqrt((v3*v3) - (vCM*vCM *pow(sin(ang3Lab),2.0))); //velocity of particle 3 in lab
          V4   = sqrt((v4*v4) + (vCM*vCM) - (2*v4*vCM*cos(phi))); //velocity of particle 4 in lab
          E3Lab = 0.5*fm3*V3*V3; // Energy in lab-3
          E4Lab = 0.5*fm4*V4*V4; // Energy in lab-4
         secEnergyLab = E3Lab;
         residueKineticEnergy = E4Lab;
         cosLab = cos(ang3Lab);
         ang4Lab = atan(sin(phi)/((vCM/v4) - cos(phi)));
         residueCosang = cos(ang4Lab); 
         } 
         if (LCT == 1) {
          productCosangLab =  cosCM;
          gamma = vCM/v3; 
          ang3Lab = acos(productCosangLab);
          fact1 = sin(ang3Lab)/(productCosangLab - gamma);
          fact2 = atan(fact1);
          cosAng3CM = cos(fact2);
          phi =  (pi - acos(cosAng3CM));
          V3   = vCM*cos(ang3Lab) + sqrt((v3*v3) - (vCM*vCM *pow(sin(ang3Lab),2.0))); //velocity of particle 3 in lab
          V4   = sqrt((v4*v4) + (vCM*vCM) - (2*v4*vCM*cos(phi))); //velocity of particle 4 in lab
          E3Lab = 0.5*fm3*V3*V3; // Energy in lab-3
          E4Lab = 0.5*fm4*V4*V4; // Energy in lab-4
          secEnergyLab = E3Lab;
          residueKineticEnergy = E4Lab;
          cosLab = cos(ang3Lab);
          ang4Lab = atan(sin(phi)/((vCM/v4) - cos(phi)));
          residueCosang = cos(ang4Lab); 
         } 
     seccosAnginLab.push_back(cosLab);
     seccosAnginLab.push_back(residueCosang);
     secEnergyinLab.push_back(secEnergyLab);
     secEnergyinLab.push_back(residueKineticEnergy);
     }else if (MF4 == 99 && MF5 == 5) {
    //This two-body kinematics is applicable for the case when incident or product particle is a gamma
     double fm1 = projectileMass;
     double fm2 = massT;
     double fm4 = residueA;
     double fm3 = massT - residueA;
     //double fqval =  recoPoint->GetQValue(ielemId, MT) ;
     double fsi  = fm1 + fm2;
     double fdi  = fm1 - fm2;
     double fsf  = fm3 + fm4;
     double fdf  = fm3 - fm4;
     double fs   = fsi * fsi + 2 * fm2 * kineticE;
     double fki2 = (fs - fsi * fsi) * (fs - fdi * fdi) / (4 * fs);
     double fkf2 = (fs - fsf * fsf) * (fs - fdf * fdf) / (4 * fs);
     double fz2  = sqrt(fki2 + fm2 * fm2);
     double fz3  = sqrt(fkf2 + fm3 * fm3);
     //double fz4 = sqrt(fkf2 + fm4 * fm4);
     double fe3 = (fz2 * fz3 - sqrt(fki2 * fkf2) / fm2) - fm3;
     //double fe4 = (fz2 * fz4 - sqrt(fki2 * fkf2)/fm2) - fm4;
     double fp12 = kineticE * kineticE + 2 * fm1 * kineticE;
     double fp32 = fe3 * fe3 + 2 * fm3 * fe3;
     //double fp42 = fe4 * fe4 + 2 * fm4 * fe4;
     cosLab = ((kineticE + fsi) * (fe3 + fm3) - fz3 * sqrt(fs)) / sqrt(fp12 * fp32);
     //double fcos4 = ((kineticE + fsi)*(fe4 + fm4)- fz4 * sqrt(fs))/sqrt(fp12*fp42);
     secEnergyLab = fe3;
    //
     seccosAnginLab.push_back(cosLab);
     secEnergyinLab.push_back(secEnergyLab);
     //cout<<"secondary energy in Lab: "<<secEnergyLab<<endl;
    }else if (MF4 == 99 && MF5 == -1 && MF6 == 6) {
      std::vector<int> law6= recoPoint->GetLawProducts(ielemId, MT);
      std::vector<int> particleZ = recoPoint->GetProductsZ(ielemId, MT);
      std::vector<int> particleA = recoPoint->GetProductsA(ielemId, MT);
      //for(int il = 0; il<law6.size();il++){
      //std::cout<<"Law::::::: "<<law6[il]<<std::endl;
      //law = law6[0];//
      //}
      int law = recoPoint->GetLaw6(ielemId, MT);
      //std::cout<<MT<<"   Law:  "<<law<<std::endl;
      //std::vector<double> particleE = rp->GetProductsEnergy(ielemId, MT,kineticE);
      int nProducts = particleZ.size();
      double fqval;
      fqval  = recoPoint->GetQValue(ielemId, MT);
      if(law == -1){
         std::cout<<"***************************************************************"<<std::endl;
         std::cout<<"*"<<" Reaction Channel No.:  "<<MT<<"                                   *"<<std::endl;
         std::cout<<"*"<<" Based on LAW =3, LAW= 4,  or LAW= 0 File-6 is not processed *"<<std::endl;
         std::cout<<"*"<<" which leads to segmentation violation                       *"<<std::endl;
         std::cout<<"***************************************************************"<<std::endl;
      }            
      switch (law) {
        case 2: {
        double fm1;  
        double fm2;
        double fm3;
        double fm4;
        double fEt;  double fm1m2;double fm3m4;double fA1;  double fB1;  double fC1;  double fD1;
        double sinLab;
        std::vector<double> productCosangCM = recoPoint->GetProductscosAng64(ielemId, MT, kineticE);
        for(int i=0; i< nProducts ; i++){
          fm1    = projectileMass;  // mass of projectile
          fm2    = massT; //mass of target
          fm3    = particleA[i]; // mass of outgoing particle
          fm4    = massT+ projectileMass - particleA[i]; // mass of daughter
          fm1m2  = fm1 + fm2;
          fm3m4  = fm3 + fm4;
          fEt    = kineticE + fqval; 
          fA1    = fm1 * fm4 * (kineticE / fEt) / (fm1m2 * fm3m4);
          fB1    = fm1 * fm3 * (kineticE / fEt) / (fm1m2 * fm3m4);
          fC1    = fm2 * fm3 * (1 + ((fm1 * fqval) / (fEt * fm2))) / (fm1m2 * fm3m4);
          fD1    = fm2 * fm4 * (1 + ((fm1 * fqval) / (fEt * fm2))) / (fm1m2 * fm3m4);
          //cout<<"sum should be 1: "<<fA1+fB1+fC1+fD1<<" "<<fA1*fC1<<" these two terms should be equal:"<<fB1*fD1<<endl;
         //Energy and Angle of light product
          secEnergyLab  = fEt * (fB1 + fD1 + 2 * sqrt(fA1 * fC1) * productCosangCM[i]);//E3
          sinLab = fD1 * productCosangCM[i] * fEt / secEnergyLab;
          cosLab        = sqrt(1 - sinLab * sinLab);
          //cout<<"kinetic Energy: "<<kineticE<<"  q-value: "<< fqval<<"  "<<secEnergyLab<< " "<<fEt<<endl;
          //std::cout<<"light product energy and angle: "<<secEnergyLab<<"  "<<cosLab<<std::endl;
          //Energy and Angle of heavy product
          double pi=acos(-1.0);
          double phi= pi - acos(productCosangCM[i]);// in CM 
          //std::cout<<acos(productCosangCM[i])<<" cosine angle: "<<phi<<" "<<cos(phi)<<std::endl;
          residueKineticEnergy =  fEt * (fA1 + fC1 + 2 * sqrt(fA1 * fC1) * cos(phi)); // E4: in Lab       
          residueCosang        =  sqrt(fm3*secEnergyLab/(fm4*residueKineticEnergy))*sinLab;
          //if(MT==103 || MT==107)std::cout<<"Heavy product energy and angle: "<<residueKineticEnergy<<"  "<<secEnergyLab<<"  "<<residueKineticEnergy+secEnergyLab<<std::endl;
          seccosAnginLab.push_back(cosLab);
          seccosAnginLab.push_back(residueCosang);
          secEnergyinLab.push_back(secEnergyLab);
          secEnergyinLab.push_back(residueKineticEnergy);
        }
        productCosangCM.clear();
        }break;
        case 3:
        break;
       case 1:{
       case 6:
            //check2:
          std::vector<double> tempCos6CM= recoPoint->GetProductscosAng6(ielemId, MT, kineticE);
          std::vector<double> tempKineticE6CM= recoPoint->GetProductsEnergy(ielemId, MT, kineticE);
          double tempEnergy = 0.0;
          for(int np=0; np< productsName.size() ; np++){
          secEnergyLab = TNudyCore::Instance()->cmToLabInelasticE(tempKineticE6CM[np], kineticE, tempCos6CM[np], massT);
          cosLab       = TNudyCore::Instance()->cmToLabInelasticCosT(secEnergyLab, tempKineticE6CM[np], kineticE,tempCos6CM[np], massT);
          seccosAnginLab.push_back(cosLab);
          secEnergyinLab.push_back(secEnergyLab);
          tempEnergy = tempEnergy + secEnergyLab;
         }
         tempCos6CM.clear();
         tempKineticE6CM.clear();
       }break;
       case 7:
       cosLab = recoPoint->GetCos6(ielemId, MT, kineticE);
       secEnergyLab = recoPoint->GetEnergy6(ielemId, MT, kineticE);
       std::vector<double> tempCos6Lab  = recoPoint->GetProductscosAng6(ielemId, MT, kineticE);
       std::vector<double> tempKineticE6Lab= recoPoint->GetProductsEnergy(ielemId, MT, kineticE);
       for(int np=0; np< productsName.size() ; np++){
         seccosAnginLab.push_back(tempCos6Lab[np]);
         secEnergyinLab.push_back(tempKineticE6Lab[np]);
       }
       tempCos6Lab.clear();
       tempKineticE6Lab.clear();
       break;
      }
      law6.clear();
      particleZ.clear();
      particleA.clear();
    }else if (MF4 == 99 && MF5 == -1 && MF6 == 99) {
    std::cout<<"***************************************************************"<<std::endl;
         std::cout<<"*"<<" Reaction Channel No.:  "<<MT<<"                                  *"<<std::endl;
         std::cout<<"*"<<" energy and angle information are not available in the ENDF  *"<<std::endl;
         std::cout<<"*"<<" data File-6 ---> which leads to segmentation violation      *"<<std::endl;
         std::cout<<"***************************************************************"<<std::endl;
         //Sometimes it is given in other reaction channel ( e.g. Reaction No. MT = 600-649). See the ENDF manual
    }
}
//______________________________________________________________________________
std::string TNudyInelastic::residueName(int Z, int A)
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
TNudyInelastic::~TNudyInelastic()
{
delete fRnd;
delete rp;
}
//______________________________________________________________________________
void TNudyInelastic::Reset()
{
  productMass.clear();
  seccosAnginLab.clear();
  secEnergyinLab.clear();
  productsName.clear();
  residueEnergy.clear(); 
  residueAngle.clear(); 
  
}
