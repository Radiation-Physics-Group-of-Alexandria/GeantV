void ExportPhysicalVolumes(TGeoVolume *top,TString guidMat,TGeoShape *MapVol[],TString *MapGuidVol, int CurrentLevel, int StartingLevel, int nn,int nnd,bool SkipOneVol,ofstream &JsonFile){

	int kLevel = 0;
  int vstep=nnd/50;
  int iVol=0;
  int nSkipTab1=0;
  int CurrentLevel=0;
  int PreviousLevel=0;
  int RelLevel=0;
 
  TGeoIterator iter(top);
  if (SkipOneVol) iter->next(); // Skip the examination of an Assembly node (examined previously)

  while ((current = iter())) {
      RelLevel=StartingLevel+iter.GetLevel()
    	nAsblVol[RelLevel]=1;

    	//print out percentange of nodes that have been examined.-------------start
    	nn++;
      if (nn>vstep){
         cout <<"#";
         fflush(stdout); 
         nn=0;
      }
      //---------------------------------------------------------------------end

    	vol = current->GetVolume();
    	shape = vol->GetShape();

      nSkipTab1=0;
      for (int j = 0; j<RelLevel); j++) {
          if (nAsblVol[j]<0) {
             nSkipTab1=nSkipTab1+nAsblVol[j];
          }
      }
      CurrentLevel=iter.GetLevel()+nSkipTab1;

    	if (vol->IsAssembly()) {
    	  if (!iter.GetMotherVolume()->IsAssembly()){
					ExportPhysicalVolumes(vol,guidMat,MapVol[],MapGuidVol,CurrentLevel,StartingLevel,nn,nnd,SkipOneVol,JsonFile);
        	iter.Skip();         
   	  	}
    	} else {
       	if ((iter.GetLevel()+nSkipTab1<MaxVisiLevel)){
         	TString chvol = CheckVolume(shape);
         	if (chvol != "") {
            //   find the guid of the material of the current volume
            guidMat = "";
            for (int i = 0; i < kMaterial; i++) {
               if (vol->GetMaterial()->GetName() == MapMat[i]) {
                  //cout << vol->GetMaterial()->GetName();
                  guidMat = MapGuidMat[i];
                  break;
               }
            }
            if (guidMat == "") {
               guidMat = MapGuidMat[0];
               cout << "ERROR **** UNDEFINED material for:" << current->GetName() << "  ASSIGNED index 0\n";
            }

        		SetTabInJson(CurrentLevel,PreviousLevel,JsonFile);
            //cout << "\n---->"<<"Volume:"<<vol->GetName()<<"\t"<<iter.GetLevel()<<" - "<<PreviousLevel<<"\t"<<shape->GetName()<<"\n";

            if (!iter.GetMotherVolume()->IsAssembly()){
            	const Double_t *MatrixRot   = vol->GetTransform()->GetRotationMatrix();
              const Double_t *MatrixTrans = vol->GetTransform()->GetTranslation();
            }else{
            	const Double_t *MatrixRot   = iter->GetMatrix()->GetRotationMatrix();
              const Double_t *MatrixTrans = iter->GetMatrix()->GetTranslation();
            }
            int nSkipLev=0;

            iVol=ExportCurrentVolume(vol,guidMat, MapVol,MapGuidVol,kVolume, iVol, shape, MatrixRot, MatrixTrans, CurrentLevel, JsonFile);
            kLevel = iter.GetLevel();
            PreviousLevel=CurrentLevel;
          }
        }
      }
    }
  }

  JsonFile << "\n";
  for (int i = StartingLevel; i < PreviousLevel; i++) {
    TabLine(PreviousLevel - i, JsonFile);
    JsonFile << "}]\n";
  }
  JsonFile << "}\n}";
}