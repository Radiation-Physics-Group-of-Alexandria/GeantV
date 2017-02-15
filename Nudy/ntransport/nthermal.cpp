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
   //const char* fENDF = "../../endffile/n-026_Fe_057.endf";
   //const char* fENDF = "../../thermal_scatt/tsl-graphite.endf";
   
   const char* fENDF = "../../thermal_scatt/tsl-HinZrH.endf";
   const char* rENDF = "test.root";
   const char* mLib = "mem.dat";
   //-----------checking---------------
   const char* irENDF="test.root";
   setEnv();
  // makeRootFile(fENDF, rENDF,op);   // <------------WORKING
   int ieleId=0;
     TNudyENDF *proc = new TNudyENDF(fENDF, rENDF, "recreate");
     // proc->SetPreProcess (0) ;
      proc->Process();
      std::string fENDFSUB = "../../fission/nfy-094_Pu_241.endf";
      proc->SetEndfSub(fENDFSUB);
      proc->Process();
      
      TNudyEndfSigma();
      double isigDiff=0.001;
      TNudyEndfSigma xsec;
      xsec.GetData(irENDF, isigDiff);
      TNudyEndfRecoPoint *rp = new TNudyEndfRecoPoint(ieleId, irENDF);
      rp->GetData(ieleId, irENDF);
      int nevent=10; 
      int mt;
      double sigma;
      double neutronEnergy =0.015;// Energy in eV
      double materialTemprature = 296.0;//  // Temprature in K
      double costhlab;
      double thermalneutronXsec;
       TNudyEndfThermal *xsecth = new TNudyEndfThermal() ;
       xsecth->GetData(irENDF, isigDiff);
       xsecth->storeElasticXsecion();
       for(int ne = 0; ne<3; ne++){
         neutronEnergy = 0.015*(ne+1)*0.1;
         thermalneutronXsec = xsecth->nThermalElasticXsecion(neutronEnergy, materialTemprature);
         //cout<<"Energy: "<< neutronEnergy<<"  x-section: "<<xsecth->nThermalElasticXsecion(neutronEnergy, materialTemprature )<<endl;
         cout<<"Main: neutron energy  "<<neutronEnergy<< " Xsec: "<< thermalneutronXsec<<endl;
         
       }
       delete xsecth;
return 0;
  }// for main
