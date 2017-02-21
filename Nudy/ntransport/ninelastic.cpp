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

  std::cout<<"Entered Main : "<< std::endl;
  setEnv();
  int ielemId;
  double nuEn=1.0E6;
  double sigmaIne;
  double isigDiff=0.001;
  const char* irENDF;
  const char* rENDF;
  const char* neutronENDFFilename;
  std::string symbolE;
  std::vector<std::string> nameX;
  std::string secParticle;
  std::string nameProc;
  std::stringstream stream;
  std::string fileName2;
  std::vector<std::string> process_names;
  int         targetZ1 = 48;
  std::string symbolT  ="Cd";
  int         targetM1 = 114;
  TNudyElement *elementProp;
  elementProp = new TNudyElement(symbolT,targetZ1,targetM1);
  neutronENDFFilename = elementProp->endfFileName();
  cout<<"file name: "<<neutronENDFFilename<<endl;
  stream << "../../endfrootfile/n" << targetZ1 << symbolT <<targetM1<<".root";
  fileName2= stream.str();
  rENDF = fileName2.c_str();
  irENDF = fileName2.c_str();
  TNudyENDF *proc = new TNudyENDF(neutronENDFFilename, rENDF, "recreate");
//proc->SetPreProcess (0) ;
  proc->Process();
  std::string fENDFSUB = "../../fission/nfy-094_Pu_241.endf";
  proc->SetEndfSub(fENDFSUB);
  proc->Process();
  TNudyEndfSigma();
  TNudyEndfSigma xsec;
  xsec.GetData(irENDF, isigDiff);
  ielemId = 0 ;
  std::stringstream str;
  std::string rootData;
  str << "../../endfrootfile/n" << targetZ1 << symbolT <<targetM1<<".root";
  rootData= str.str();
  const char* rootENDF = rootData.c_str();
  cout<<"Reading x-section file:\t"<<rootENDF<<endl;
  double xsec1;
  double energy1;
  std::vector<std::string> productsName;
  std::vector<double> secKineticEnergyLab;
  std::vector<double> seccosAngLab;
  TNudyInelastic *nProcIne= new TNudyInelastic(ielemId,rootENDF);
  int reactionMT;
  for(int i = 0; i <1; i++){//energy loop
    nuEn = 10.0E6;// + 0.25E6*i;// 0.002E6*i + 0.1E6;
//  n-Inelastic process
    for(int ie = 0 ; ie<100; ie++){ // event loop
      nProcIne->nInelasticXsec(nuEn, elementProp); //provide neutron energy
    	nProcIne->GetInelasticParameters(xsec1,energy1);
    	sigmaIne = nProcIne->GetInelasticXsec();
    	productsName = nProcIne->GetParticleInelastic();
    	secKineticEnergyLab = nProcIne->GetKiEnInelastic();
    	seccosAngLab        = nProcIne->GetcosAngInelastic();
    	reactionMT = nProcIne->reactionChannelNumber();
      std::cout<<"Incident neutron kinetic energy: "<< nuEn<<"  total n-Inelastic cross-section: "<<sigmaIne<<std::endl;
      std::cout << "<----------------------------------------------------->" << std::endl;
    	for(unsigned int j=0; j<productsName.size(); j++){
//      std::cout<<"no. of secondaries: "<<productsName.size()<<"   "<<secKineticEnergyLab.size()<<std::endl;
     	  std::cout<<"event No.: "<<ie <<"  Reaction Channel No. MT: "<<reactionMT<<"  Products Name: "<<productsName[j]
     	            <<""<<" Kinetic energy_Lab: "<<secKineticEnergyLab[j]<<""
     	           <<" cosine angle_lab: "<<seccosAngLab[j]<<std::endl;
     	}
     	nProcIne->Reset();
    }//event loop
    productsName.clear();
    secKineticEnergyLab.clear();
  }//energy loop*/
return 0;
}// for main
