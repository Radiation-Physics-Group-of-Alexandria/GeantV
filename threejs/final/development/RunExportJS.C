#include "TGeoManager.h"
#include "TGeoVolume.h"
#include "TGeoNode.h"
#ifndef VAR_STORE
#define  VAR_STORE
#include "varStore.h"
#endif


void RunExportJS(int MaxVisiLevel=4,TString ChosenTopVolume="Top",TString FileName=""){
   //MaxVisiLevel -> only volumes positioned in the tree level up to MaxVisiLevel will be considered
   //ChosenTopVolume -> Only the volumes that belong to the tree branch starting by ChosenTopVolume will be considered
   //OneCopy -> if true:  Only the first copy number of the volume and its branch will be considered
   //FileName --\> json outpout will be stored in the file: FileName

   // Load some root geometry

   TGeoNode *TopVol;
   TGeoVolume *vol;
   TGeoShape *shape;
   TString path;
   varStore prm;
      Int_t ns=0;

   //GetListOfVolumes()
   //GetListOfShapes()

   if (prm.nMaxVol==0) {
      TObjArray  *myvol=gGeoManager->GetListOfVolumes();
      TIter next(myvol);
      while ((vol=(TGeoVolume*)next())) {
         ns++;
         if (vol->GetNumber()>prm.nMaxVol)prm.nMaxVol=vol->GetNumber();
      }
   }
    
   for (int i=0;i<100000;i++) {
      prm.CdF[i]=false;
      prm.LastParent[i]=-1;
      prm.newCloneId[i]=-1;
      prm.MapVol[i]=true;
   }

   for (int i=0;i<100000;i++) prm.MapVol[i]=true;

   if (FileName=="")FileName="Cogevito.json";
   if (FileName.Index(".json")<=0)FileName=FileName+".json";

   //cout <<"Output File name: "<<FileName<<"\n";
   
   if (!gGeoManager ) {
      cout << "\nERROR    NO GEOMETRY DEFINED gGeoManager=null   STOP !!!!\n";
      return(0);
   }
   TGeoVolume *top = gGeoManager->GetTopVolume();
   //gGeoManager->DefaultColors();

   if (ChosenTopVolume!="Top"){
      top=gGeoManager->GetVolume(ChosenTopVolume);
   }
   if (!top) {
         cout << "\nBAD LUCK MY FRIEND, Top volume "<<ChosenTopVolume<<" does not exist !!!!\n";
         return;
   }
   Int_t nTotalNodes = gGeoManager->CountNodes(top,1000,0);

   cout << "\n\n-->Top volume: " << top->GetName() << " nNodes: " << nTotalNodes << " MaxVisiLevel: "<<MaxVisiLevel<<"\n\n";


   prm.kVolume=0;
   prm.top=top; 
   prm.MaxVisiLevel=MaxVisiLevel;
   prm.PrimFilename=FileName;
   prm.ClonesFilename=prm.PrimFilename(0,prm.PrimFilename.Index("."))+"-clones.json";
   prm.nTotalNodes=gGeoManager->CountNodes(top,1000,0);
   
   ExportPhyVols(&prm);

   return(0);
}

//============== JSON GEOMETRY ======== End

