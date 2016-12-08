#include <iostream>
#include <string>
#include <fstream>
#include <sstream>
#include <cmath>
#include <vector>
#include <iterator>

#include <Riostream.h>
#include <TTree.h>
#include <TCollection.h>
#include <TIterator.h>
#include <TFile.h>
#include <TMath.h>
#include <TString.h>
#include <TClass.h>
#include <TObject.h>
#include <TROOT.h>
#include <TSystem.h>
#include <TCollection.h>
#include <Rtypes.h>
#include <TCanvas.h>
#include <TNtuple.h>
#include <TTreeReader.h>
#include <TKey.h>
#include <TMacro.h>
#include <TGeoManager.h>
#include <TGeoNode.h>
#include <TGeoVolume.h>

#include "TNudyDB.h"
#include "TNudyElementTable.h"
#include "TNudyENDF.h"
#include "TNudyEndfCont.h"
#include "TNudyEndfFile.h"
#include "TNudyEndfINTG.h"
#include "TNudyEndfList.h"
#include "TNudyEndfMat.h"
#include "TNudyEndfTape.h"
#include "TNudyEndfRecord.h"
#include "TNudyEndfSec.h"
#include "TNudyEndfTab1.h"
#include "TNudyEndfTab2.h"
#include "TNudyLibrary.h"
#include "TNudyManager.h"
#include "TNudySubLibrary.h"
#include "TVNudyModel.h"
#include "TNudyEndfEnergy.h"
#include "TNudyEndfSigma.h"
//#include "TNudySampling.h"
#include "TNudyEndfRecoPoint.h"
#include "TNudyEndfTape.h"
#include "TNudyEndfAng.h"
#include "Particle.h"
#include "ElementProp.h"
#include "G4ParticleTable.hh"
#include "G4ParticleTable.hh"
#include "G4ParticleDefinition.hh"
#include "TNudyEndfThermal.h"
#include "TNudyInelastic.h"
#include "TNudyElastic.h"
#include "TNudyFission.h"
#include "TNudyCapture.h"
#include   "TNudyElastic.h"
#include   "TNudyElement.h"

#include   "Element.h"
#include   "Material.h"
#include   "SystemOfUnits.h"
#include   "PhysicalConstants.h"
#include   "Isotope.h"

class TFile;
class TNudyEndfTape;
class TNudyEndfMat;
class TNudyEndfFile;
class TNudyEndfSec;
class TNudyEndfCont;
class TNudyEndfList;
class TNudyEndfTab1;
class TNudyEndfTab2;
class TNudyEndfINTG;
class TNudyENDF;
class TNudyEndfAng;
class TNudyEndfEnergy;
class TNudyEndfSigma;
//class TNudySampling;
class TNudyEndfRecoPoint;
class TNudyEndfTape;
class TNudyEndfThermal;
class TNudyInelastic;
class TNudyElastic;
class TNudyFission;
class TNudyCapture;
class TNudyElastic;
class TNudyElement;
//class Element;
using namespace std;


 void setEnv(void) {
  gSystem->Load("libRIO");
  gSystem->Load("libHist");
  gSystem->Load("libGeom");
  gSystem->Load("libGraf");
  gSystem->Load("libGpad");
  gSystem->Load("libEG");
  gSystem->Load("libMathMore");
  gSystem->Load("libNudy");
  gSystem->Load("/usr/lib64/libgfortran.so");
}
//---------------------------------------------//
std::string atomicNumbertostring(int Z){
  int targetZ = Z;
  stringstream ss;
  ss << targetZ;
  string strZ = ss.str();
  int len=strZ.length(); 
  std::string str4;
  if(len==1){
  str4="00" + strZ;
  }else{str4="0"+strZ;}
  return str4;
 }

std::string atomicMasstostring(int A){

  int targetA = A;
  stringstream ssa;
  ssa << targetA;
  string strA = ssa.str();
  int lenA=strA.length(); 
  std::string strAf;
  if(lenA==1){
    strAf="00" + strA;
  }else if(lenA==2){strAf="0"+strA;
  }else{strAf=strA;}
  return strAf;

}

std::string myreplace(std::string &s,
                       std::string &toReplace,
                       std::string &replaceWith)
{
    return(s.replace(s.find(toReplace), toReplace.length(), replaceWith));
}

const char* fileName(std::string S,int tZ, int A1){
  std::string base="/home/shiba/endffile/n-ddd_XX_fff.endf";
  std::string str1;
  std::string str2;
  std::string str3;
  std::string name;
  int zTarget=tZ;
  std::string Zf;
  Zf=atomicNumbertostring(zTarget);
  //cout<<"atomic number to string:"<<Zf<<endl;
  int mTarget=A1;
  std::string Af;
  Af=atomicMasstostring(mTarget);
  //cout<<"atomic mass to string:"<<Af<<endl;
  std::string s = "XX";
  std::string symbol = S; 
  std::string Z1 = "ddd";
  std::string m1 = "fff";
  std::string m = Af;
  str1=myreplace(base, s, symbol);
  str2=myreplace(str1, Z1, Zf);
  str3=myreplace(str2, m1,m);
  const char* file=str3.c_str();
  return file;
}


//-------------------Main-------------------//
int main(){
    // G4ParticleDefinition* particle1
    //      = G4ParticleTable::GetParticleTable()->FindParticle("neutron");
   std::cout<<"Entered Main : "<< std::endl;
   TNudyENDF();
   TNudyENDF obj1;
   TNudyEndfRecoPoint *obj2;
   bool menu = true;
   int choice;
   Int_t matNum;
   const char* particle = "neutron";
   Int_t zVal, aVal, za;   // zVal = Z, aVal = A, za = 100 * Z + A
   Int_t op = 4;           // Option for Root file creation verbosity
   std::vector<std::vector<double> > CSTable;
// *************** Declare filenames ************************
   //const char* fENDF = "/home/shiba/endffile/n-026_Fe_057.endf";
   //const char* fENDF = "/home/shiba/thermal_scatt/tsl-graphite.endf";
   
   const char* fENDF = "/home/shiba/thermal_scatt/tsl-HinZrH.endf";
   const char* rENDF = "test.root";
   const char* mLib = "mem.dat";
   //-----------checking---------------
   const char* irENDF="test.root";
   setEnv();
  // makeRootFile(fENDF, rENDF,op);   // <------------WORKING
   int ieleId=0;
 //----------Lines added 22/09/2016--------------- 
     TNudyENDF *proc = new TNudyENDF(fENDF, rENDF, "recreate");
     // proc->SetPreProcess (0) ;
      proc->Process();
      std::string fENDFSUB = "/home/shiba/fission/nfy-094_Pu_241.endf";
      proc->SetEndfSub(fENDFSUB);
      proc->Process();
      Particle *obj3=new Particle[1];
      
      double en=obj3[0].energy;
      //cout<<"Energy of incident neutron :::::::::\t"<<en<<endl;
      TNudyEndfSigma();
      double isigDiff=0.001;
      TNudyEndfSigma xsec;
      xsec.GetData(irENDF, isigDiff);
      //--------added on 8th-oct-2016-----------
     // int mat=9543;
      //int mt=4;
     //TNudyEndfFile *file = new TNudyEndfFile(9543,4);
      //TNudyEndfAng *objf = new TNudyEndfAng(file);
      //--------added on 8th-oct-2016-----------
      TNudyEndfRecoPoint *rp = new TNudyEndfRecoPoint(ieleId, irENDF);
      rp->GetData(ieleId, irENDF);
      int nevent=10; 
      int mt;
      double sigma;
      double En = 1.79177;
      
      double costhlab;
       TNudyEndfThermal xsecth;
       xsecth.GetData(irENDF, isigDiff);
       cout<<"Thermal neutron x-section: "<<xsecth.nThermalElasticXsecion(En)<<endl;
return 0;
  }// for main
