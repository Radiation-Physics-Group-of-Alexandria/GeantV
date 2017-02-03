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
   int targetZ; 
   int targetM;  
   double nuEn=1.0E6;
   double sigmaIne;
   double costhlab;
   double secE;
   double isigDiff=0.001;
   const char* fENDFn; 
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
   //int         targetZ1 = 30;
   //std::string symbolT  ="Zn";
   //int         targetM1 = 64;
   
   //int         targetZ1 = 57;
   //std::string symbolT  ="La";
   //int         targetM1 = 138;
   
   //int         targetZ1 = 25;
   //std::string symbolT  ="Mn";
   //int         targetM1 = 55;
   
   int         targetZ1 = 26;
   std::string symbolT  ="Fe";
   int         targetM1 = 56;
   
   //int         targetZ1 = 26;
   //std::string symbolT  ="Fe";
   //int         targetM1 = 54;
   
   //int         targetZ1 = 34;
   //std::string symbolT  ="Se";
   //int         targetM1 = 74;
   
   //int         targetZ1 = 3;
   //std::string symbolT  ="Li";
   //int         targetM1 = 6;
      
   //For MT = 25, file-4 and file-5
   //int         targetZ1 = 3;
   //std::string symbolT  ="Li";
   //int         targetM1 = 7;
   
   //int         targetZ1 = 27;
   //std::string symbolT  ="Co";
   //int         targetM1 = 58;
  
   //int         targetZ1 = 64;
   //std::string symbolT  ="Gd";
   //int         targetM1 = 154;
   
    //int         targetZ1 = 14;
    //std::string symbolT  ="Si";
    //int         targetM1 = 29;
   
   //int         targetZ1 = 6;
   //std::string symbolT  ="C";
   //int         targetM1 = 12;
   
   //int         targetZ1 = 20;
   //std::string symbolT  ="Ca";
   //int         targetM1 = 40;
   
   //int         targetZ1 = 24;
   //std::string symbolT  ="Cr";
   //int         targetM1 = 52;
   
    //int         targetZ1 = 35;
    //std::string symbolT  ="Br";
    //int         targetM1 = 81;
    
    //int         targetZ1 = 47;
    //std::string symbolT  ="Ag";
    //int         targetM1 = 107;
   
    //int         targetZ1 = 79;
    //std::string symbolT  ="Au";
    //int         targetM1 = 197;
    
    //int         targetZ1 = 20;
    //std::string symbolT  ="Ca";
    //int         targetM1 = 48;
    
   
   TNudyElement *elementProp;
   elementProp = new TNudyElement(symbolT,targetZ1,targetM1);
   neutronENDFFilename = elementProp->endfFileName();
   cout<<"file name: "<<neutronENDFFilename<<endl;
   stream << "/home/shiba/geantOct16/endfrootfile/n" << targetZ1 << symbolT <<targetM1<<".root";
            fileName2= stream.str();
            rENDF = fileName2.c_str();
            irENDF = fileName2.c_str();
            TNudyENDF *proc = new TNudyENDF(elementProp->endfFileName(), rENDF, "recreate");
            //proc->SetPreProcess (0) ;
            proc->Process();
            std::string fENDFSUB = "/home/shiba/fission/nfy-094_Pu_241.endf";
            proc->SetEndfSub(fENDFSUB);
            proc->Process();
            TNudyEndfSigma();
            TNudyEndfSigma xsec;
            xsec.GetData(irENDF, isigDiff);
   //ofstream fout;
   //fout.open("/home/shiba/output/Fe54_EnAng_20MeV.txt",ios::out);
   ielemId = 0 ;
   std::stringstream str;
   std::string rootData;
   str << "/home/shiba/geantOct16/endfrootfile/n" << targetZ1 << symbolT <<targetM1<<".root";
   rootData= str.str();
   const char* rootENDF = rootData.c_str();
   cout<<"Reading x-section file:\t"<<rootENDF<<endl;
   double xsec1;
   double energy1;
   std::vector<std::string> productsName;
   std::vector<double> secKineticEnergyLab;
   std::vector<double> seccosAngLab;
   TNudyInelastic *nProcIne= new TNudyInelastic(ielemId,rootENDF); //n-Capture
   int reactionMT;
   for(int i = 0; i <1; i++){//energy loop
       nuEn=19.0E6;// + 0.25E6*i;// 0.002E6*i + 0.1E6;
       //n-Inelastic process
       for(int ie = 0 ; ie<1; ie++){ // event loop
    	  nProcIne->nInelasticXsec(nuEn, elementProp); //provide neutron energy
    	  nProcIne->GetInelasticParameters(xsec1,energy1);
    	  sigmaIne = nProcIne->GetInelasticXsec();
    	  productsName = nProcIne->GetParticleInelastic();
    	  secKineticEnergyLab = nProcIne->GetKiEnInelastic();
    	  seccosAngLab        = nProcIne->GetcosAngInelastic();
    	  reactionMT = nProcIne->reactionChannelNumber();
    	  sigmaIne = nProcIne->GetInelasticXsec();
     	  for(int j=0; j<productsName.size(); j++){
           /*if (reactionMT ==17){
           std::cout<<"event No. "<<ie <<"Reaction Channel No. MT: "<<reactionMT<<"  Products Name: "<<productsName[j]
     	            <<""<<" Kinetic energy_Lab: "<<secKineticEnergyLab[j]<<""
     	           <<" cosine angle_lab: "<<seccosAngLab[j]<<std::endl;//}*/
     	          // if(reactionMT ==16 && productsName[j] == "n")fout<<secKineticEnergyLab[j]<<"\t"<<seccosAngLab[j]<<endl;
     	 }
     	 nProcIne->Reset();
       }//event loop
    }//energy loop*/
return 0;
  }// for main
