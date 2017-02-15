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
#include "TNudyElastic.h"
#include "TNudyElement.h"

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
//______________________________________________________________________________

//-------------------Main-----------------------------------------------------//
int main(){
   
   std::cout<<"Entered Main : "<< std::endl;
   setEnv();   
   
   int ielemId;
   int targetZ; 
   int targetM;  
   int eventcount=0;
   int num_process_names=0;
   double nuEn=1.0E6;
   double sigmaEla;
   double En;
   double costhlab;
   std::vector<double> secE;
   std::vector<double> secpcosAng;
   double fissFragmass;
   double fissFragEn;
   double isigDiff=0.001;
   std::vector<double> xsecs ={0};
   const char* neutronENDFFilename; 
   const char* irENDF;
   const char* rENDF;
   std::string symbolE;
   std::vector<std::string> nameX;
   std::vector<std::string> secParticle;
   std::string nameProc;
   std::stringstream stream;
   std::string fileName2;
   std::vector<std::string> process_names;
  
   int         targetZ1 = 26;
   std::string symbolT  ="Fe";
   int         targetM1 = 56;

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
   
   
       ofstream fout;
       fout.open("../../output/junk.txt",ios::out);
       ielemId = 0 ;
       std::stringstream str;
       std::string rootData;
       str << "../../endfrootfile/n" << targetZ1 << symbolT <<targetM1<<".root";
       rootData= str.str();
       const char* rootENDF = rootData.c_str();
       cout<<"Reading x-section file:\t"<<rootENDF<<endl;
       TNudyCapture *nProc= new TNudyCapture(ielemId,rootENDF);     //n-Elastic
    //n-Elastic process
   for(int i = 0; i <1; i++){//event loop
        nuEn= 1.0E-5;//2.5E6 + 0.02E6*i;
        for(int j = 0; j <1; j++){ // energy loop
           nProc->nCaptureXsec(nuEn, elementProp); //provide neutron energy
           secParticle = nProc->nCaptureGetsecParticles();
           secE	      = nProc->nCaptureGetEnergy();
           secpcosAng  = nProc->nCaptureGetcosAng();
           for(int k = 0; k<secParticle.size(); k++){
             cout<<"Event No. : "<<i <<"   Products name:  "<< secParticle[k]<<"  cosAng_Lab:   "<<secpcosAng[k]<<"  Energy_Lab(eV) "<<secE[k]<<endl;
           }
         nProc->nCaptureProcessReset();      
        }// energy loop
    }//event loop
return 0;
  }// for main
