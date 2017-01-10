#include "TGeoManager.h"
#include "Riostream.h"
#include "TString.h"
#include "TGeoVolume.h"

void RunExportJson(int MaxVisiLevel=4,TString ChosenTopVolume="Top",TString FileName="")
{
   // Load some root geometry
   TGeoVolume *top = gGeoManager->GetTopVolume();
   TGeoNode *current;
   TGeoNode *TopVol;
   TGeoVolume *vol;
   TGeoShape *shape;
   TString path;
   TGeoVolume *CurrentVolume;
   TString guid;
   ofstream JsonFile;
   TString JsonFileName = "geo-box.json";
   int MapVol1    [300000];
   TString       MapVol    [30000];
   TString       MapGuidVol[30000];
   TString       MapGuidMat[1000];


   int nAsblVol[100];
   Int_t kVolume = 0;

   // Materials (for semplicity) are simply defined by the LineColour of the volume/node
   // Color N has an GUID identifier stored in MapGuidMat[N]
   for (int i=0;i<1000;i++){
      MapGuidMat[i]="";
   }

   if (FileName!="") {
      JsonFileName=FileName;
   }

   if (ChosenTopVolume!="Top"){
      TGeoNode *TopVol=top->GetNode(ChosenTopVolume);
      if (!TopVol) {
         cout << "\nBAD LUCK MY FRIEND, Top volume "<<ChosenTopVolume<<" does not exist !!!!\n";
         return;
      }
      top=TopVol->GetVolume();
   } 
   Int_t nTotalNodes = gGeoManager->CountNodes(top,10000,0);

   //ExportJsonTree();

   cout << "\n\n-->Top volume: " << top->GetName() << " Nodes: " << nTotalNodes << "\n\n";

   JsonFile.open(JsonFileName);

   //Initialize Json file
   JsonFile << "{\n\"metadata\": { \n\"version\": 4.3,\"type\": \"Object\",\n\"generator\": \"ObjectExporter\"},\n\"geometries\": [\n";


   // Check volumes in the tree structure
   ofstream JsonTree;
   JsonTree.open("Volumes1.txt");
   //CheckLogicalVolumes(top,MaxVisiLevel,nAsblVol,0,0,0,nTotalNodes, MapGuidMat,MapVol1,MapGuidVol,JsonTree);
   cout <<"\n Count volumes ...\n";

   TIter next(gGeoManager->GetListOfVolumes());
   kass=0;
   kVolume = 0;
   while ((CurrentVolume = (TGeoVolume *)next())) {
      
      TString chvol = CheckVolume(CurrentVolume->GetShape());
      if (chvol != "") {
         if (!CurrentVolume->GetShape()->IsAssembly()) {
            kass++;
         }
         kVolume++;
         JsonTree << kVolume << "\t"<<CurrentVolume->GetNumber()<< "\t"<<CurrentVolume->GetShape()->GetPointerName()<<"\t"<< CurrentVolume->GetNtotal()<<"\t"<<CurrentVolume->GetShape()->ClassName()<<"\n";

      }
    }
   JsonTree.close();
   //   End definition of Geometries

   cout <<"\n"<<kVolume<<" "<<kass<<" Volumes have been defined Ntotal: "<<"\n";


   return(0);

   //   Start definition of geometries
   cout <<"\n ******   START EXPORTING Logical Volumes   *****\n";
   //cout << "  10%  20%  30%  40%  50%  60%  70%  80%  90% 100%\n";

   kVolume=ExportLogicalVolumes(top,MaxVisiLevel,nAsblVol,0,0,0,nTotalNodes, MapGuidMat,MapVol,MapGuidVol,JsonFile);

   JsonFile << "]," << "\n"; 
   cout <<"\n"<<kVolume<<" Geometries have been defined in the JSON file: "<< JsonFileName<<"\n";

   //   End definition of Geometries

    
   //   Start definition of Materials
   JsonFile << "\"materials\": [" << "\n"; 
   cout  << "\n Start definition of Materials .... ";

   TIter nextmat(gGeoManager->GetListOfMaterials());
   TGeoMaterial *lmaterial;

   int kMaterial = 0;
   for (int i=0;i<100;i++){
      if (MapGuidMat[i]!=""){
         if (kMaterial > 0) JsonFile << ",\n";
         guid = doGuidStuff(GuidGenerator());
         MapGuidMat[i] = guid;
         ExportMaterials(i,guid,JsonFile);
         kMaterial++;
      }
   } 
   cout  << kMaterial<< " Materials definedin the JSON file: "<< JsonFileName<<"\n";

 
   JsonFile << "]," << "\n"; 
   // End definition of Materials


   // Start definition of scene and Physical volumes/nodes

   JsonFile << "\"object\": {" << "\n"; 
   JsonFile << "\"uuid\": \"" << doGuidStuff(GuidGenerator()) << "\",\n";
   JsonFile << "\"type\": \"Scene\",\n";
   JsonFile << "\"name\": \"Scene\",\n";
   JsonFile << "\"matrix\": [1,0,0,0,0,1,0,0,0,0,1,0,0,0,0,1]";

   cout <<"\n ******   START EXPORTING Physical Volumes   *****\n";
   cout << "  10%  20%  30%  40%  50%  60%  70%  80%  90% 100%\n";

   int nPhysicalVol= ExportPhysicalVolumes(top,MaxVisiLevel,nAsblVol,0,0,0,nTotalNodes,MapGuidMat,MapVol,MapGuidVol,kVolume,JsonFile);

   cout <<"\n"<<nPhysicalVol<<" Physical Volumes have been defined in the JSON file: "<< JsonFileName<<"\n";

   //   End definition of Scene and Physical volumes/nodes

   JsonFile << "}\n}";

   JsonFile.close(); // close the file handle

   std::cout << "\n\n--> JSON FILE: " << JsonFileName << " CREATED " << std::endl << std::endl;

}

//============== JSON GEOMETRY ======== End

