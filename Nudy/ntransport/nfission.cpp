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


//-------------------Main-------------------//
int main(){
   
   std::cout<<"Entered Main : "<< std::endl;
   setEnv();   
   
   int ielemId;
   int targetZ; 
   int targetM;  
   int eventcount=0;
   int num_process_names=0;
   double nuEn=1.0E6;
   double sigmaEla,sigmaFis,sigmaCap,sigmaIne;
   double En;
   double costhlab;
   double secE;
   double secpcosAng;
   double fissFragmass;
   double fissFragEn;
   double isigDiff=0.001;
   std::vector<int> fissionfragmentsMass;
   std::vector<int> fissionfragmentsCharge;
   std::vector<double> neutronEnergy;
   std::vector<double> neutronAngle;  
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

   int         targetZ1 = 94;
   std::string symbolT  ="Pu";
   int         targetM1 = 241;
   
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
            
   ofstream fout1;
   //fout1.open("/home/shiba/FissionFragments/Pu241_nsEnergyAngle.txt",ios::out);
   
   ofstream fout2;
   //fout2.open("/home/shiba/FissionFragments/Pu241_FFsmasscharge.txt",ios::out);
   
   ielemId = 0 ;
   std::stringstream str;
   std::string rootData;
   str << "/home/shiba/geantOct16/endfrootfile/n" << targetZ1 << symbolT <<targetM1<<".root";
   rootData= str.str();
   const char* rootENDF = rootData.c_str();
   cout<<"Reading x-section file:\t"<<rootENDF<<endl;
   
   
   
   TNudyFission *nProcFis= new TNudyFission(ielemId,rootENDF); //n-Fission
   for(int i = 0; i <1; i++){//energy loop
        nuEn= 2.0E6 + 1.0E6*i;
    	//n-Fission process
    	for(int j=0; j<1 ; j++){//event loop
        nProcFis->nFissionXsec(nuEn, elementProp); //provide neutron energy
     	//cout<< nuEn<<"\t"<<nProcFis->GetFissionXsec()<<"\t"<<endl;
     	//cout<<"prompt neutron: "<<nProcFis->GetPromptneutron()<<" fission fragment mass "<< nProcFis->GetFissionFragmentmass()<<endl;//nProcFis->GetFissionProcessName()
     	
    	//nProcFis->nFissionXsec(nuEn, elementProp); 
    	neutronEnergy = nProcFis->GetKineticEnergyNeutron();
    	neutronAngle  = nProcFis->GetcosAngNeutron();
    	for(int ie=0; ie<neutronEnergy.size(); ie++){     
    	//cout<<ie<<"  neutron energy and angle produced in fission process: "<<neutronEnergy[ie]<<"  "<<neutronAngle[ie]<<endl;
    	
    	//fout1<<neutronEnergy[ie]<<"\t"<<neutronAngle[ie]<<endl;
    	
    	}   
    	fissionfragmentsMass   = nProcFis->GetFissionFragmentsmass();
    	fissionfragmentsCharge = nProcFis->GetFissionFragmentscharge();
    	for(int im=0; im<fissionfragmentsMass.size(); im++){
    	//cout<<"event no. " <<j<<"   fission fragments mass : "<<fissionfragmentsMass[im]<<" charge:  "<<fissionfragmentsCharge[im]<<endl;
    	//fout2<<fissionfragmentsMass[im]<<"\t"<<fissionfragmentsCharge[im]<<endl;
    	
    	}
    	nProcFis->processReset();
    	}//event loop
     }//energy loop*/
return 0;
  }// for main
