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
   std::vector<double> xsecs ={0};
   const char* fENDFn; 
   const char* irENDF;
   const char* rENDF;
   
   std::string symbolE;
   std::vector<std::string> nameX;
   std::string secParticle;
   std::string nameProc;
   std::stringstream stream;
   std::string fileName2;
   std::vector<std::string> process_names;
// Create some NIST materials 
  using geantphysics::Material;
  using geantphysics::Element;
  using geantphysics::Isotope;

     // NIST amaterial hydrogen
     //    Material *nist_mat_H     = Material::NISTMaterial("NIST_MAT_H");
  // NIST amaterial aluminum
     //Material *nist_mat_Al    = Material::NISTMaterial("NIST_MAT_Al");
  // NIST amaterial carbon
  // Material *nist_mat_C    = Material::NISTMaterial("NIST_MAT_C");
  // NIST amaterial oxygen
    // Material *nist_mat_O    = Material::NISTMaterial("NIST_MAT_O");  
  // NIST amaterial silicon
  //Material *nist_mat_Si    = Material::NISTMaterial("NIST_MAT_Si");
  
  // NIST amaterial iron
    //Material *nist_mat_Fe    = Material::NISTMaterial("NIST_MAT_Fe");
  
  // NIST amaterial water
  //Material *nist_mat_Water = Material::NISTMaterial("NIST_MAT_WATER");
  // NIST amaterial Concrete
  //Material *nist_mat_Concr = Material::NISTMaterial("NIST_MAT_CONCRETE");
  
  // NIST amaterial Gadolinium
    //Material *nist_mat_Gd     = Material::NISTMaterial("NIST_MAT_Gd");
  
  
  // NIST amaterial Plutonium
    //  Material *nist_mat_Pu     = Material::NISTMaterial("NIST_MAT_Pu");
  
  
  // printing out the material table (i.e. all materials that were created)
  //
  //std::cout<< Material::GetTheMaterialTable();
 /* const std::vector<Material*> theMaterialT = Material::GetTheMaterialTable();
  for (unsigned int imat=0; imat<theMaterialT.size(); ++imat) {
     Material *mat = theMaterialT[imat];
     const std::vector<Element*> theElemVect = mat->GetElementVector();
     //std::cout<< " === Material name = " << mat->GetName() 
          //    << " has "<< theElemVect.size() << " elements:" << std::endl;
     for (unsigned int ielem=0; ielem<theElemVect.size(); ++ielem) {
          Element *elem =theElemVect[ielem];
          const std::vector<Isotope*> theIsoVect = elem->GetIsotopeVector();
          const double* theIsoAb = elem->GetRelativeAbundanceVector();
          //std::cout<< "    ----> " << ielem <<"-th element name is " << elem->GetName()<<"\t"<<elem->GetZ() 
            //       << " has "<< theIsoVect.size() << " isotopes:"<<std::endl;
          symbolE = elem->GetSymbol();
          targetZ = elem->GetZ(); 
          for (unsigned int iiso=0; iiso<theIsoVect.size(); ++iiso) {
            Isotope *theIso = theIsoVect[iiso];
            targetM = theIso->GetN(); //provides number of nucleons 
            //std::cout<< "          ***** "<< theIso->GetName()<<"\t"<<theIso->GetZ() 
            //         << "          ---- "<<theIso->GetN()<<"\t"<<theIso->GetA()<<"\t"<<
            //         " with rel. ab. = " <<theIsoAb[iiso]/geant::perCent << " [%] "<<std::endl; 
            //cout<<"targetM----->"<< theIso->GetA()<<endl;
            fENDFn=fileName(symbolE,targetZ,targetM);
            cout<<"file name is:\t"<<fENDFn<<endl;
            std::stringstream stream;
            std::string fileName2;
            stream << "/home/shiba/geantOct16/endfrootfile/n" << targetZ << symbolE <<targetM<<".root";
            fileName2= stream.str();
            rENDF = fileName2.c_str();
            irENDF = fileName2.c_str();
            TNudyENDF *proc = new TNudyENDF(fileName(symbolE,targetZ,targetM), rENDF, "recreate");
            //proc->SetPreProcess (0) ;
            proc->Process();
            std::string fENDFSUB = "/home/shiba/fission/nfy-094_Pu_241.endf";
            proc->SetEndfSub(fENDFSUB);
            proc->Process();
            TNudyEndfSigma();
            TNudyEndfSigma xsec;
            xsec.GetData(irENDF, isigDiff);
          }          
     }
  } */
   int         targetZ1 = 94;
   std::string symbolT  ="Pu";
   int         targetM1 = 241;
   
   fENDFn=fileName(symbolT,targetZ1,targetM1);
            cout<<"file name is:\t"<<fENDFn<<endl;
           // std::stringstream stream;
           // std::string fileName2;
            stream << "/home/shiba/geantOct16/endfrootfile/n" << targetZ1 << symbolT <<targetM1<<".root";
            fileName2= stream.str();
            rENDF = fileName2.c_str();
            irENDF = fileName2.c_str();
            TNudyENDF *proc = new TNudyENDF(fileName(symbolT,targetZ1,targetM1), rENDF, "recreate");
            //proc->SetPreProcess (0) ;
            proc->Process();
            std::string fENDFSUB = "/home/shiba/fission/nfy-094_Pu_241.endf";
            proc->SetEndfSub(fENDFSUB);
            proc->Process();
            TNudyEndfSigma();
            TNudyEndfSigma xsec;
            xsec.GetData(irENDF, isigDiff);
   // ofstream fout;
   //fout.open("/home/shiba/output/Fe56n_nalpha.txt",ios::out);
   ielemId = 0 ;
   std::stringstream str;
   std::string rootData;
   str << "/home/shiba/geantOct16/endfrootfile/n" << targetZ1 << symbolT <<targetM1<<".root";
   rootData= str.str();
   const char* rootENDF = rootData.c_str();
   cout<<"Reading x-section file:\t"<<rootENDF<<endl;
   TNudyElement *elementProp;
   elementProp = new TNudyElement(symbolT,targetZ1,targetM1);
   TNudyFission *nProcFis= new TNudyFission(ielemId,rootENDF); //n-Fission
   for(int i = 0; i <1; i++){//energy loop
        nuEn= 2.0E6 + 1.0E6*i;
    	//n-Fission process
        nProcFis->nFissionXsec(nuEn, elementProp); //provide neutron energy
     	//cout<< nuEn<<"\t"<<nProcFis->GetFissionXsec()<<"\t"<<endl;
     	//cout<<"prompt neutron: "<<nProcFis->GetPromptneutron()<<" fission fragment mass "<< nProcFis->GetFissionFragmentmass()<<endl;//nProcFis->GetFissionProcessName()
     	
    	     // nProcFis->GetParticleFission()<<"\t"<<nProcFis->GetKiEnFission()<<"\t"<<
    	//for(int j=0; j<1 ; j++){
    	//nProcFis->nFissionXsec(nuEn, elementProp);     
    	//fout<<nProcFis->GetKiEnFission()<<"\t"<<nProcFis->GetcosAngFission()<<endl;
    	//}
     }//energy loop*/
return 0;
  }// for main
