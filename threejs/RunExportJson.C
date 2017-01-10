#include "TGeoManager.h"
#include "Riostream.h"
#include "TString.h"
#include "TGeoVolume.h"

void RunExportJson(int MaxVisiLevel=4,TString ChosenTopVolume="Top",bool OneCopy=true,TString FileName="")
{
   //MaxVisiLevel -> only volumes positioned in the tree level up to MaxVisiLevel will be considered
   //ChosenTopVolume -> Only the volumes that belong to the tree branch starting by ChosenTopVolume will be considered
   //OneCopy -> if true:  Only the first copy number of the volume and its branch will be considered
   //FileName --\> json outpout will be stored in the file: FileName

   // Load some root geometry
   TGeoNode *current;
   TGeoNode *TopVol;
   TGeoVolume *vol;
   TGeoShape *shape;
   TString path;
   TGeoVolume *CurrentVolume;
   TString guid,Vguid;
   ofstream JsonFile;
   TString JsonFileName = "geo-box.json";
   TGeoShape *MapVol1  [100000];
   TGeoShape *MapVol   [100000];
   TString   MapGuidVol[100000];
   TString   MapGuidMat[30000];
   
   if (!gGeoManager ) {
      cout << "\nERROR    NO GEOMETRY DEFINED gGeoManager=null   STOP !!!!\n";
      return(0);
   }
   TGeoVolume *top = gGeoManager->GetTopVolume();
   int nAsblVol[100];
   Int_t kVolume = 0;
   int nPhyVol=0;


   // Materials (for semplicity) are simply defined by the LineColour of the volume/node
   // Color N has an GUID identifier stored in MapGuidMat[N]
   for (int i=0;i<30000;i++){
      MapGuidMat[i]="";
   }

   if (FileName!="") {
      JsonFileName=FileName;
   }

   //cout << "\n\n-->Geometry for: " << top->GetName() << " Nodes: " << gGeoManager->GetNNodes() << "\n";

   if (ChosenTopVolume!="Top"){
      top=gGeoManager->GetVolume(ChosenTopVolume);
      //TGeoNode *TopVol=top->FindNode(ChosenTopVolume);
      //if (!TopVol) {
      //top=TopVol->GetVolume();
   }
   if (!top) {
         cout << "\nBAD LUCK MY FRIEND, Top volume "<<ChosenTopVolume<<" does not exist !!!!\n";
         return;
   }
   Int_t nTotalNodes = gGeoManager->CountNodes(top,1000,0);

   //ExportJsonTree();

   cout << "\n\n-->Top volume: " << top->GetName() << " Nodes: " << nTotalNodes << "Max Visi Lev:"<<MaxVisiLevel<<"\n\n";

   JsonFile.open(JsonFileName);


   //Initialize Json file
   JsonFile << "{\n\"metadata\": { \n\"version\": 4.3,\"type\": \"Object\",\n\"generator\": \"ObjectExporter\"},\n\"geometries\": [\n";


   // Check volumes in the tree structure
   ofstream JsonTree;
   JsonTree.open("VolumesInAlice-number3.txt");

   int pl=0;
   for (int k=0;k<100;k++){nAsblVol[k]=0;}
   int kvol=ExportJsonTree(top,MaxVisiLevel,nAsblVol,0,pl,0,nTotalNodes, MapGuidMat,MapVol,MapGuidVol,nPhyVol,OneCopy,JsonTree);
   JsonTree.close();


   nPhyVol=0;
   //   Start definition of geometries


   cout <<"\n ******   START EXPORTING Logical Volumes   *****\n";
   pl=0;
   for (int k=0;k<100;k++){nAsblVol[k]=0;}

   //====================================
   //Create an box geometry to be used for Assembly Volumes
   MapGuidVol[kVolume] = doGuidStuff(GuidGenerator());
   JsonFile << "{\n\t\"uuid\": \""<<MapGuidVol[kVolume]<<"\",\n\t\"type\": \"BoxGeometry\",\n\t\"width\": 0.001,\n\t\"height\": 0.001,\n\t\"depth\": 0.001,\n\t\"widthSegments\": 1,\n\t\"heightSegments\": 1,\n\t\"depthSegments\": 1\n}";
   //====================================
   kVolume++;
   kVolume=ExportLogicalVolumes(top,MaxVisiLevel,nAsblVol,0,pl,kVolume,nTotalNodes, MapGuidMat,MapVol,MapGuidVol,nPhyVol,OneCopy,JsonFile);

   JsonFile << "]," << "\n"; 
   cout <<"\n"<<kVolume<<" Geometries have been defined in the JSON file: "<< JsonFileName<<"\n";
   cout <<"\n"<<nPhyVol<<" Physical volumes have been counted.\n";

   //   End definition of Geometries

    
   //   Start definition of Materials
   JsonFile << "\"materials\": [" << "\n"; 
   cout  << "\n Start definition of Materials .... ";

   //  TIter nextmat(gGeoManager->GetListOfMaterials());
   //  TGeoMaterial *lmaterial;

   int kMaterial = 0;

   //=======================================
   // We define a transparent material for vacuum and assmbly volumes
   Vguid = doGuidStuff(GuidGenerator());
   JsonFile << "{\n\t\"uuid\": \"" << Vguid << "\",\n";
   JsonFile << "\t\"type\": \"MeshLambertMaterial\",\n";
   JsonFile << "\t\"name\":\"Transparent \",\n";

   JsonFile << "\t\"color\":0,\n";
   JsonFile << "\t\"emissive\":0,\n";

   JsonFile << "\t\"opacity\": 0.0,\n";
   JsonFile << "\t\"transparent\": true,\n";
   JsonFile << "\t\"wireframe\": true\n},\n";

   //=======================================

   for (int i=0;i<30000;i++){
      //cout <<"\t\t\t\n Loop: "<<i<<" "<<MapGuidMat[i]<<"\n";
      if (MapGuidMat[i]!=""){
         //cout <<"\n Color: "<<i<<" "<<MapGuidMat[i]<<"\n";

         if (kMaterial > 0) JsonFile << ",\n";
         guid = doGuidStuff(GuidGenerator());
         MapGuidMat[i] = guid;
         ExportMaterials(i,guid,JsonFile);
         kMaterial++;
      }
   } 
   cout  << kMaterial<< " Materials defined in the JSON file: "<< JsonFileName<<"\n";

   JsonFile << "]," << "\n"; 
   // End definition of Materials

   // Start definition of scene and Physical volumes/nodes

   JsonFile << "\"object\": {" << "\n"; 
   JsonFile << "\"uuid\": \"" << doGuidStuff(GuidGenerator()) << "\",\n";
   JsonFile << "\"type\": \"Scene\",\n";
   JsonFile << "\"name\": \"Scene\",\n";
   JsonFile << "\"matrix\": [1,0,0,0,0,1,0,0,0,0,1,0,0,0,0,1]";

   //cout <<"\n ******   START EXPORTING Physical Volumes   *****\n"<<top->GetName();
   //cout << "  10%  20%  30%  40%  50%  60%  70%  80%  90% 100%\n";

   pl=0;
   for (int k=0;k<100;k++){nAsblVol[k]=0;}
   int nPhysicalVol= ExportPhysicalVolumes(top,MaxVisiLevel,nAsblVol,0,pl,0,nPhyVol,MapGuidMat,MapVol,MapGuidVol,kVolume,OneCopy,Vguid,JsonFile);
   for (int i = 0; i < pl; i++) {
    TabLine(pl - i, JsonFile);
    JsonFile << "}]\n";
   }

   cout <<"\n"<<nPhysicalVol<<" Physical Volumes have been defined in the JSON file: "<< JsonFileName<<"\n";

   //   End definition of Scene and Physical volumes/nodes

   JsonFile << "}\n}";

   JsonFile.close(); // close the file handle

   std::cout << "\n\n--> JSON FILE: " << JsonFileName << " CREATED " << std::endl << std::endl;

   return(0);
}

//============== JSON GEOMETRY ======== End

