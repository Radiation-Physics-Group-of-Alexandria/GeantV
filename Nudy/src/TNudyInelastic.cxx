// Various interaction process of neutron: Selection of secondary particle type, its energy and angle
//The final state product residues mass and charges are assigned considering the incident particle
//as a neutron. The charge and mass of the residue will change depending upon incident particle as
//a charge particle or photon.
#include <iostream>
#include "TNudyEndfEnergyAng.h"
#include "ElementProp.h"
#include "TNudyCore.h"
#include "TNudyEndfRecoPoint.h"
#include "TNudyInelastic.h"
#include "TNudyEndfMat.h"
#include "TNudyElement.h"

#ifdef USE_ROOT
#include "TRandom3.h"
#endif
#include "TCanvas.h"
#include "TFile.h"
using namespace std;
#ifdef USE_ROOT
ClassImp(TNudyInelastic)
#endif

TNudyInelastic::TNudyInelastic(){}

TNudyInelastic::TNudyInelastic(int ElemId, const char* rootfileName)
{
  ielemId =  ElemId;// for element number like say for H=1, U=2, the way we number the elements for material having composite elements
  elemId=ElemId;
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
  kineticE = nuEn;
  int reactionMT;
  double sumXsec=0;
  charge = targetElement->GetatomicNumber();
  mass   = targetElement->GetatomicMass();
  int count1=0;
  for(unsigned int crsp = 0; crsp < nMTs; crsp++){
    reactionMT=rp->MtValues[ielemId][crsp];
    if(reactionMT != 2 && reactionMT !=18 && reactionMT !=102 ){
      procNameIne   = "nInelastic";
      name        =  "2n";
      sumXsec=sumXsec+rp->GetSigmaPartial(ielemId, crsp, kineticE);
      sigmaInelastic=sumXsec;
      //cout<<"inelastic xec: "<<sumXsec<<endl;
      mtTemp.push_back(reactionMT);
      crssIel.push_back(rp->GetSigmaPartial(elemId, crsp, kineticE));
      nMTinelastic[elemId][count1] =reactionMT;
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
  MT  = nMTinelastic[elemId][isel];//rp->MtValues[ielemId][isel];
  MF4 = rp->GetMt4(ielemId, MT);
  MF5 = rp->GetMt5(ielemId, MT);
  MF6 = rp->GetMt6(ielemId, MT);
  LCT = rp->GetCos4Lct(ielemId, MT);
  cout<<"MT: "<<MT<<" MF4: "<<MF4<<" MF5: "<<MF5<<" MF6: "<<MF6<<" LCT: "<<LCT<<endl;
  //if(MT==16)cout<<"MT:16"<<endl;
  div_t divr;
  switch (MT){
    //case 2:{
     // }break;
      case 11: // 2nd
      residueA = mass - 3;
      residueZ = charge - 1;
      name="2nd";
      En = secEnergyLab;
      costhlab=cosLab;
      break;
      case 16: // 2n
      residueA = mass - 1;
      residueZ = charge;
      GetSecParameter(targetElement, rp);
      name="2n";
      En = secEnergyLab;
      costhlab=cosLab;
      //cout<<"case:16_2n\t"<<cosLab <<"\t"<< secEnergyLab<<endl;
      break;
      case 17: // 3n
      residueA = mass - 2;
      residueZ = charge;
      GetSecParameter(targetElement, rp);
      name="3n";
      En = secEnergyLab;
      costhlab=cosLab;
      //cout<<"case:17_3n\t"<<cosLab <<"\t"<< secEnergyLab<<endl;
      break;
      //case 18:{}
      //break;
      case 22: // n+alpha
      residueA = mass - 4;
      residueZ = charge - 2;
      GetSecParameter(targetElement, rp);
      name="n+alpha";
      En = secEnergyLab;
      costhlab=cosLab;
      cout<<"case:22_n_alpha\t"<<cosLab <<"\t"<< secEnergyLab<<endl;
      break;
      case 23: // n+3a
      residueA = mass - 12;
      residueZ = charge - 6;
      GetSecParameter(targetElement, rp);
      name="n+3alpha";
      En = secEnergyLab;
      costhlab=cosLab;
      break;
      case 24: // 2n+a
      residueA = mass - 5;
      residueZ = charge - 2;
      GetSecParameter(targetElement, rp);
      name="2n+alpha";
      En = secEnergyLab;
      costhlab=cosLab;
      break;
      case 25: // 3n+a
      residueA = mass - 6;
      residueZ = charge - 2;
      GetSecParameter(targetElement, rp);
      name="3n+alpha";
      En = secEnergyLab;
      costhlab=cosLab;
      break;
      case 28: // n+p
      residueA = mass - 1;
      residueZ = charge - 1;
      GetSecParameter(targetElement, rp);
      name="n+p";
      En = secEnergyLab;
      costhlab=cosLab;
      break;
      case 29: // n +2a
      residueA = mass - 8;
      residueZ = charge - 4;
      GetSecParameter(targetElement, rp);
      name="n+2alpha";
      En = secEnergyLab;
      costhlab=cosLab;
      break;
      case 30: // 2n+2a
      residueA = mass - 9;
      residueZ = charge - 4;
      GetSecParameter(targetElement, rp);
      name="2n+2alpha";
      En = secEnergyLab;
      costhlab=cosLab;
      break;
      case 32: // n+d
      residueA = mass - 2;
      residueZ = charge - 1;
      GetSecParameter(targetElement, rp);
      name="n+d";
      En = secEnergyLab;
      costhlab=cosLab;
      break;
      case 33: // n+t
      residueA = mass - 3;
      residueZ = charge - 1;
      GetSecParameter(targetElement, rp);
      name="n+t";
      En = secEnergyLab;
      costhlab=cosLab;
      break;
      case 34: // n + He3
      residueA = mass - 3;
      residueZ = charge - 2;
      GetSecParameter(targetElement, rp);
      name="n+He3";
      En = secEnergyLab;
      costhlab=cosLab;
      break;
      case 35: // n+d+2a
      residueA = mass - 10;
      residueZ = charge - 5;
      GetSecParameter(targetElement, rp);
      name="n+d+2alpha";
      En = secEnergyLab;
      costhlab=cosLab;
      break;
      case 36: // n+t+2a
      residueA = mass - 11;
      residueZ = charge - 5;
      GetSecParameter(targetElement, rp);
      name="n+t+2alpha";
      En = secEnergyLab;
      costhlab=cosLab;
      break;
      case 37: // 4n
      residueA = mass - 3;
      residueZ = charge;
      GetSecParameter(targetElement, rp);
      name="4";
      En = secEnergyLab;
      costhlab=cosLab;
      break;
      case 41: // 2n+p
      residueA = mass - 2;
      residueZ = charge - 1;
      GetSecParameter(targetElement, rp);
      name="2n+p";
      En = secEnergyLab;
      costhlab=cosLab;
      break;
      case 42: // 3n+p
      residueA = mass - 3;
      residueZ = charge - 1;
      GetSecParameter(targetElement, rp);
      name="3n+p";
      En = secEnergyLab;
      costhlab=cosLab;
      break;
      case 44: // n+2p
      residueA = mass - 2;
      residueZ = charge - 1;
      GetSecParameter(targetElement, rp);
      name="n+2p";
      En = secEnergyLab;
      costhlab=cosLab;
      break;
      case 45: // n+p+a
      residueA = mass - 5;
      residueZ = charge - 3;
      GetSecParameter(targetElement, rp);
      name="n+p+a";
      En = secEnergyLab;
      costhlab=cosLab;
      //
      break;
      case 50://production of neutron and residue at the ground state
      name="n+residue1stexci";
      case 51://production of neutron and residue at the 1st excited state, this will continue upto 90
      
      cout<<"En " <<En<<" angle_lab "<<costhlab<<endl;
      case 52:
      case 53:
      case 54:
      case 55:
      case 56:
      case 57:
      case 58:
      case 59:
      case 60:
      case 61:
      case 62:
      case 63:
      case 64:
      case 65:
      case 66:
      case 67:
      case 68:
      case 69:
      case 70:
      case 71:
      case 72:
      case 73:
      case 74:
      case 75:
      case 76:
      case 77:
      case 78:
      case 79:
      case 80:
      case 81:
      case 82:
      case 83:
      case 84:
      case 85:
      case 86:
      case 87:
      case 88:
      case 89:
      case 90:
      case 91:
      residueA = mass;
      residueZ = charge;
      GetSecParameter(targetElement, rp);
      name="n_continum";
      En = secEnergyLab;
      costhlab=cosLab;
      //std::cout <<"91:::::::::::\t"<< cosLab <<"  "<<secEnergyLab << std::endl;
      break;
      //case 102:{
      //}break;
      case 103: // p
      residueA = mass;
      residueZ = charge - 1;
      GetSecParameter(targetElement, rp);
      En = secEnergyLab;
      costhlab=cosLab;
      name="p";
      break;
      case 104: // d
      residueA = mass - 1;
      residueZ = charge - 1;
      GetSecParameter(targetElement, rp);
      name="d";
      En = secEnergyLab;
      costhlab=cosLab;
      break;
      case 105: // t
      residueA = mass - 2;
      residueZ = charge - 1;
      GetSecParameter(targetElement, rp);
      name="t";
      En = secEnergyLab;
      costhlab=cosLab;
      break;
      case 106: // He3
      residueA = mass - 2;
      residueZ = charge - 2;
      GetSecParameter(targetElement, rp);
      name="He3";
      En = secEnergyLab;
      costhlab=cosLab;
      break;
      case 107: // alpha
      residueA = mass - 3;
      residueZ = charge - 2;
      GetSecParameter(targetElement, rp);
      name="alpha";
      En = secEnergyLab;
      costhlab=cosLab;
      //cout<<"MT---------107:\t"<<cosLab<<"\t"<<secEnergyLab<<endl;
      break;
      case 108: // 2a
      residueA = mass - 7;
      residueZ = charge - 4;
      GetSecParameter(targetElement, rp);
      name="2alpha";
      En = secEnergyLab;
      costhlab=cosLab;
      break;
      case 109: // 3a
      residueA = mass - 11;
      residueZ = charge - 6;
      GetSecParameter(targetElement, rp);
      name="3alpha";
      En = secEnergyLab;
      costhlab=cosLab;
      break;
      case 111: // 2p
      residueA = mass - 1;
      residueZ = charge - 2;
      GetSecParameter(targetElement, rp);
      name="2p";
      En = secEnergyLab;
      costhlab=cosLab;
      break;
      case 112: // p+a
      residueA = mass - 4;
      residueZ = charge - 3;
      GetSecParameter(targetElement, rp);
      name="p+alpha";
      En = secEnergyLab;
      costhlab=cosLab;
      break;
      case 113: // t+2a
      residueA = mass - 10;
      residueZ = charge - 5;
      GetSecParameter(targetElement, rp);
      name="t+2alpha";
      En = secEnergyLab;
      costhlab=cosLab;
      break;
      case 114: // d+2a
      residueA = mass - 9;
      residueZ = charge - 5;
      GetSecParameter(targetElement, rp);
      name="d+2alpha";
      En = secEnergyLab;
      costhlab=cosLab;
      break;
      case 115: // p+d
      residueA = mass - 2;
      residueZ = charge - 2;
      GetSecParameter(targetElement, rp);
      name="p+d";
      En = secEnergyLab;
      costhlab=cosLab;
      break;
      case 116: // p+t
      residueA = mass - 3;
      residueZ = charge - 2;
      GetSecParameter(targetElement, rp);
      name="p+t";
      En = secEnergyLab;
      costhlab=cosLab;
      break;
      case 117: // d+a
      residueA = mass - 5;
      residueZ = charge - 3;
      GetSecParameter(targetElement, rp);
      name="d+alpha";
      En = secEnergyLab;
      costhlab=cosLab;
      break;
      case 152: // 5n
      residueA = mass - 4;
      residueZ = charge;
      GetSecParameter(targetElement, rp);
      name="5n";
      En = secEnergyLab;
      costhlab=cosLab;
      break;
      case 153: // 6n
      residueA = mass - 5;
      residueZ = charge;
      GetSecParameter(targetElement, rp);
      name="6n";
      En = secEnergyLab;
      costhlab=cosLab;
      break;
      case 154: // 2n+t
      residueA = mass - 4;
      residueZ = charge - 1;
      GetSecParameter(targetElement, rp);
      name="2n+t";
      En = secEnergyLab;
      costhlab=cosLab;
      break;
      case 155: // t+a
      residueA = mass - 6;
      residueZ = charge - 3;
      GetSecParameter(targetElement, rp);
      name="t+a";
      En = secEnergyLab;
      costhlab=cosLab;
      break;
      case 156: // 4n+p
      residueA = mass - 4;
      residueZ = charge - 1;
      GetSecParameter(targetElement, rp);
      name="4n+p";
      En = secEnergyLab;
      costhlab=cosLab;
      break;
      case 157: // 3n+d
      residueA = mass - 4;
      residueZ = charge - 1;
      GetSecParameter(targetElement, rp);
      name="3n+d";
      En = secEnergyLab;
      costhlab=cosLab;
      break;
      case 158: // n+d+a
      residueA = mass - 6;
      residueZ = charge - 3;
      GetSecParameter(targetElement, rp);
      name="n+d+a";
      En = secEnergyLab;
      costhlab=cosLab;
      break;
      case 159: // 2n+p+a
      residueA = mass - 6;
      residueZ = charge - 3;
      GetSecParameter(targetElement, rp);
      name="2n+p+a";
      En = secEnergyLab;
      costhlab=cosLab;
      break;
      case 160: // 7n
      residueA = mass - 6;
      residueZ = charge;
      GetSecParameter(targetElement, rp);
      name="7n";
      En = secEnergyLab;
      costhlab=cosLab;
      break;
      case 161: // 8n
      residueA = mass - 7;
      residueZ = charge;
      GetSecParameter(targetElement, rp);
      name="8n";
      En = secEnergyLab;
      costhlab=cosLab;
      break;
      case 162: // 5np
      residueA = mass - 5;
      residueZ = charge - 1;
      GetSecParameter(targetElement, rp);
      name="5n+p";
      En = secEnergyLab;
      costhlab=cosLab;
      break;
      case 163: // 6np
      residueA = mass - 6;
      residueZ = charge - 1;
      GetSecParameter(targetElement, rp);
      name="6n+p";
      En = secEnergyLab;
      costhlab=cosLab;
      break;
      case 164: // 7np
      residueA = mass - 7;
      residueZ = charge - 1;
      GetSecParameter(targetElement, rp);
      name="7n+p";
      En = secEnergyLab;
      costhlab=cosLab;
      break;
      case 165: // 4n+a
      residueA = mass - 7;
      residueZ = charge - 2;
      GetSecParameter(targetElement, rp);
      name="4n+a";
      En = secEnergyLab;
      costhlab=cosLab;
      break;
      case 166: // 5na
      residueA = mass - 8;
      residueZ = charge - 2;
      GetSecParameter(targetElement, rp);
      name="5n+a";
      En = secEnergyLab;
      costhlab=cosLab;
      break;
      case 167: // 6na
      residueA = mass - 9;
      residueZ = charge - 2;
      GetSecParameter(targetElement, rp);
      name="6n+a";
      En = secEnergyLab;
      costhlab=cosLab;
      break;
      case 168: // 7na
      residueA = mass - 10;
      residueZ = charge - 2;
      GetSecParameter(targetElement, rp);
      name="7n+a";
      En = secEnergyLab;
      costhlab=cosLab;
      break;
      case 169: // 4nd
      residueA = mass - 5;
      residueZ = charge - 1;
      GetSecParameter(targetElement, rp);
      name="4n+d";
      En = secEnergyLab;
      costhlab=cosLab;
      break;
      case 170: // 5nd
      residueA = mass - 6;
      residueZ = charge - 1;
      GetSecParameter(targetElement, rp);
      name="5n+d";
      En = secEnergyLab;
      costhlab=cosLab;
      break;
      case 171: // 6nd
      residueA = mass - 7;
      residueZ = charge - 1;
      GetSecParameter(targetElement, rp);
      name="6n+d";
      En = secEnergyLab;
      costhlab=cosLab;
      break;
      case 172: // 3nt
      residueA = mass - 5;
      residueZ = charge - 1;
      GetSecParameter(targetElement, rp);
      name="3n+t";
      En = secEnergyLab;
      costhlab=cosLab;
      break;
      case 173: // 4nt
      residueA = mass - 6;
      residueZ = charge - 1;
      GetSecParameter(targetElement, rp);
      name="4n+t";
      En = secEnergyLab;
      costhlab=cosLab;
      break;
      case 174: // 5nt
      residueA = mass - 7;
      residueZ = charge - 1;
      GetSecParameter(targetElement, rp);
      name="5n+t";
      En = secEnergyLab;
      costhlab=cosLab;
      break;
      case 175: // 6nt
      residueA = mass - 8;
      residueZ = charge - 1;
      GetSecParameter(targetElement, rp);
      name="6n+t";
      En = secEnergyLab;
      costhlab=cosLab;
      break;
      case 176: // 2n+He3
      residueA = mass - 4;
      residueZ = charge - 2;
      GetSecParameter(targetElement, rp);
      name="2n+He3";
      En = secEnergyLab;
      costhlab=cosLab;
      break;
      case 177: // 3n + He3
      residueA = mass - 5;
      residueZ = charge - 2;
      GetSecParameter(targetElement, rp);
      name="3n+He3";
      En = secEnergyLab;
      costhlab=cosLab;
      break;
      case 178: // 4n +He3
      residueA = mass - 6;
      residueZ = charge - 2;
      GetSecParameter(targetElement, rp);
      name="4n+He3";
      En = secEnergyLab;
      costhlab=cosLab;
      break;
      case 179: // 3n2p
      residueA = mass - 4;
      residueZ = charge - 2;
      GetSecParameter(targetElement, rp);
      name="3n+2p";
      En = secEnergyLab;
      costhlab=cosLab;
      break;
      case 180: // 3n2a
      residueA = mass - 10;
      residueZ = charge - 4;
      GetSecParameter(targetElement, rp);
      name="3n+2a";
      En = secEnergyLab;
      costhlab=cosLab;
      break;
      case 181: // 3npa
      residueA = mass - 7;
      residueZ = charge - 3;
      GetSecParameter(targetElement, rp);
      name="3n+p+a";
      En = secEnergyLab;
      costhlab=cosLab;
      break;
      case 182: // dt
      residueA = mass - 4;
      residueZ = charge - 2;
      GetSecParameter(targetElement, rp);
      name="d+t";
      En = secEnergyLab;
      costhlab=cosLab;
      break;
      case 183: // npd
      residueA = mass - 3;
      residueZ = charge - 2;
      GetSecParameter(targetElement, rp);
      name="n+p+d";
      En = secEnergyLab;
      costhlab=cosLab;
      break;
      case 184: // npt
      residueA = mass - 4;
      residueZ = charge - 2;
      GetSecParameter(targetElement, rp);
      name="n+p+t";
      En = secEnergyLab;
      costhlab=cosLab;
      break;
      case 185: // ndt
      residueA = mass - 5;
      residueZ = charge - 2;
      GetSecParameter(targetElement, rp);
      name="n+d+t";
      En = secEnergyLab;
      costhlab=cosLab;
      break;
      case 186: // npHe3
      residueA = mass - 4;
      residueZ = charge - 3;
      GetSecParameter(targetElement, rp);
      name="n+p+He3";
      En = secEnergyLab;
      costhlab=cosLab;
      break;
      case 187: // ndHe3
      residueA = mass - 5;
      residueZ = charge - 3;
      GetSecParameter(targetElement, rp);
      name="n+d+He3";
      En = secEnergyLab;
      costhlab=cosLab;
      break;
      case 188: // ntHe3
      residueA = mass - 6;
      residueZ = charge - 3;
      GetSecParameter(targetElement, rp);
      name="n+t+He3";
      En = secEnergyLab;
      costhlab=cosLab;
      break;
      case 189: // nta
      residueA = mass - 7;
      residueZ = charge - 3;
      GetSecParameter(targetElement, rp);
      name="n+t+a";
      En = secEnergyLab;
      costhlab=cosLab;
      break;
      case 190: // 2n2p
      residueA = mass - 3;
      residueZ = charge - 2;
      GetSecParameter(targetElement, rp);
      name="2n+2p";
      En = secEnergyLab;
      costhlab=cosLab;
      break;
      case 191: // pHe3
      residueA = mass - 3;
      residueZ = charge - 3;
      GetSecParameter(targetElement, rp);
      name="p+He3";
      En = secEnergyLab;
      costhlab=cosLab;
      break;
      case 192: // dHe3
      residueA = mass - 4;
      residueZ = charge - 3;
      GetSecParameter(targetElement, rp);
      name="d+He3";
      En = secEnergyLab;
      costhlab=cosLab;
      break;
      case 193: // aHe3
      residueA = mass - 6;
      residueZ = charge - 4;
      GetSecParameter(targetElement, rp);
      name="a+He3";
      En = secEnergyLab;
      costhlab=cosLab;
      break;
      case 194: // 4n2p
      residueA = mass - 5;
      residueZ = charge - 2;
      GetSecParameter(targetElement, rp);
      name="4n+4p";
      En = secEnergyLab;
      costhlab=cosLab;
      break;
      case 195: // 4n2a
      residueA = mass - 11;
      residueZ = charge - 4;
      GetSecParameter(targetElement, rp);
      name="4n+2a";
      En = secEnergyLab;
      costhlab=cosLab;
      break;
      case 196: // 4npa
      residueA = mass - 8;
      residueZ = charge - 3;
      GetSecParameter(targetElement, rp);
      name="4n+p+a";
      En = secEnergyLab;
      costhlab=cosLab;
      break;
      case 197: // 3p
      residueA = mass - 2;
      residueZ = charge - 3;
      GetSecParameter(targetElement, rp);
      name="3p";
      En = secEnergyLab;
      costhlab=cosLab;
      break;
      case 198: // n3p
      residueA = mass - 3;
      residueZ = charge - 3;
      GetSecParameter(targetElement, rp);
      name="n+3p";
      En = secEnergyLab;
      costhlab=cosLab;
      break;
      case 199: // 3n2pa
      residueA = mass - 8;
      residueZ = charge - 4;
      GetSecParameter(targetElement, rp);
      name="3n+2p+a";
      En = secEnergyLab;
      costhlab=cosLab;
      break;
      case 200: // 5n2p
      residueA = mass - 6;
      residueZ = charge - 2;
      GetSecParameter(targetElement, rp);
      name="5n+2p";
      En = secEnergyLab;
      costhlab=cosLab;
      break;
      }
      mtTemp.clear();
      crssIel.clear();
      crss.clear();
   }
//______________________________________________________________________________
std::string TNudyInelastic::GetInelasticProcessName()
{
     return procNameIne;
}
double TNudyInelastic::GetInelasticXsec()
{
       return sigmaInelastic;
}
std::string TNudyInelastic::GetParticleInelastic()
{
     return name;
}
//______________________________________________________________________________
double TNudyInelastic::GetKiEnInelastic()
{
       return En;
}
//______________________________________________________________________________
double TNudyInelastic::GetcosAngInelastic()
{
      return costhlab;
}


void TNudyInelastic::GetInelasticParameters(double &xsec1, double &E)
{

   xsec1 = sigmaInelastic;
      E = En;
       //return sigmaInelastic;
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
    double fm4 = residueA;
    double fm3 = massT - residueA;
    double fqval =  recoPoint->GetQValue(ielemId, MT)/1E6 ;
   // double fm3 = massT + fm1  - residueA - fqval;
    //if(fm3==0)fm3= fqval;// - (fm1+fm2) + fm4;
    double u =1;// 931;
    double fsi  = fm1 + fm2;
    double fdi  = fm1 - fm2;
    double fsf  = fm3 + fm4;
    double fdf  = fm3 - fm4;
    //cout<<"Mass of particles:"<<fm1<<"\t"<<fm2<<"\t"<<fm3<<"\t"<<fm4<<"\t"<<MT<<"  "<<recoPoint->GetQValue(ielemId, MT)<< " frame of reference: "<<LCT<<endl;
    //cout<<fsi<<" "<<fdi<<" "<<fsf<<" "<<fdf<<endl;
    double fs   = fsi * fsi + 2 * fm2 * kineticE;
    double fki2 = (fs - fsi * fsi) * (fs - fdi * fdi) / (4 * fs);
    double fkf2 = (fs - fsf * fsf) * (fs - fdf * fdf) / (4 * fs);
    double fz2  = sqrt(fki2 + fm2 * fm2);
    double fz3  = sqrt(fkf2 + fm3 * fm3);
    double fz4 = sqrt(fkf2 + fm4 * fm4);
    double fe3 = (fz2 * fz3 - sqrt(fki2 * fkf2) / fm2) - fm3;
    double fe4 = (fz2 * fz4 - sqrt(fki2 * fkf2)/fm2) - fm4;
    // double fp12 = kineticE * kineticE + 2 * fm1 * kineticE;
    // double fp32 = fe3 * fe3 + 2 * fm3 * fe3;
    // double fp42 = fe4 * fe4 + 2 * fm4 * fe4;
    // double fcos3 = ((kineticE + fsi)*(fe3 + fm3)- fz3 * sqrt(fs))/sqrt(fp12*fp32);
    // double fcos4 = ((kineticE + fsi)*(fe4 + fm4)- fz4 * sqrt(fs))/sqrt(fp12*fp42);
   //if(cosCM>1)std::cout<<"cos CM = "<< cosCM << std::endl;
    // secEnergyCM = recoPoint->GetEnergy5(ielemId, MT, kineticE);
    secEnergyLab = fe3;
    //cout<<"energy: "<< secEnergyLab<< "  "<<fe4<<endl;
    
    //In case of MT 50-90, the residue is in excited state, if the incident particle 
    //is same as the out going particle, but the target is in excited state, 
    //then the Q-value is same as the excitation energy, however if we express Q=m1+2-(m3+m4) = 0,
    // in this case, but the actual Q value is the Ex
            fm3   = 1;
     double Q     = fm1 + fm2 - fm3 - fm4;
     double Ex    = recoPoint->GetQValue(ielemId, MT);
     double Ecin  = kineticE*fm2/(fm1+fm2);
     double Ecout = Ecin + Q + Ex;
     double Ec3   = Ecout*fm4/(fm3+fm4);
     double Ec4   = Ecout*fm3/(fm3+fm4);
     secEnergyLab = TNudyCore::Instance()->cmToLabInelasticE(Ec3, kineticE, cosCM, massT);
    //cout<< Q <<"  "<< Ex << " " << Ecin<< "  "<< Ecout<< " "<< Ec3<< " "<<Ec4<< " "<<
    //                  secEnergyLab<<"\t"<<cosCM<<endl;
    cout<< "B4 secEnergyCM:  "<<secEnergyCM<<endl;
    secEnergyCM = Ec3;
    cout<< "After secEnergyCM:  "<<secEnergyCM<<endl;
    if (LCT == 2) {
      cosLab = TNudyCore::Instance()->cmToLabInelasticCosT(secEnergyLab, secEnergyCM, kineticE, cosCM,
                                                           massT);
     } else if (LCT == 1) {
       cosLab = cosCM;
     }
   }else if (MF4 == 99 && MF5 == 5) {
    double fm1 = 1;
    double fm2 = massT;
    double fm4 = residueA;
    double fm3 = massT - residueA;
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
    int  zd = recoPoint->GetZd6(ielemId, MT); 
    int  ad = recoPoint->GetAd6(ielemId, MT);
    
     //std::cout<<"law "<< law <<" zd: "<<zd<<" ad: "<<ad<<std::endl;
    switch (law) {
    case 2: {
      cosCM         = recoPoint->GetCos64(ielemId, MT, kineticE);
      double fm1    = 1;
      double fm2    = massT;
      double fm3    = massT - residueA;
      double fm4    = residueA;
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
      cout<<"Law2:Mass of particles:"<<massT<<"\t"<<fm1<<"\t"<<fm2<<"\t"<<fm3<<"\t"<<fm4<<"\t"<<MT<<endl;
       //std::cout<< cosLab <<"  "<<secEnergyLab << std::endl;
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
       
       cout<<"Law6------->"<<massT<<"\t"<<cosCM<<"\t"<<cosLab<<"\t"<<secEnergyCM<<"\t"<<secEnergyLab<<endl;
      break;
    case 7:
      cosLab = recoPoint->GetCos6(ielemId, MT, kineticE);
      // std::cout<< cosLab << std::endl;
      secEnergyLab = recoPoint->GetEnergy6(ielemId, MT, kineticE);
      // std::cout<< secEnergyLab << std::endl;
      cout<<"Law7------->"<<massT<<"\t"<<cosCM<<"\t"<<cosLab<<"\t"<<secEnergyCM<<"\t"<<secEnergyLab<<endl;
      break;
    }
  }
// }


}
void TNudyInelastic::FillHisto(double icosLab, double isecEnergyLab)
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
TNudyInelastic::~TNudyInelastic()
{
delete fRnd;
delete rp;
}
