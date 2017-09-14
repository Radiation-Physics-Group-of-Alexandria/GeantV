  int numCells  = 25;
  int numRows   = 25;
  int numLayers = 25;

 //define world
  vecgeom::UnplacedBox *worldUnplaced      = new vecgeom::UnplacedBox(fWorldSizeX/2, fWorldSizeYZ/2, fWorldSizeYZ/2); 
  vecgeom::LogicalVolume *world = new vecgeom::LogicalVolume("world", worldUnplaced);
  world->SetRegion(worldRegion);
  world->SetMaterialPtr(fWorldMaterial); 

// absorbers
  vecgeom::UnplacedBox *fAbsBoxes[i];   // two different thicknesses: x17 (2.4 mm) x8 (4.8 mm)
  for (int i=1; i<=maxAbsorbers; i++) {  
      fAbsBoxes[i] = new vecgeom::UnplacedBox("Absorber", fAbsThickness[i]/2,fCaloSizeXY/2,fCaloSizeXY/2);
  }
  vecgeom::LogicalVolume *fAbsLogic[maxAbsorbers];

  char *volName = new char;
  char *cellInRowName = new char;
  char *rowInLayerName = new char;
  char *layerName = new char;

  double zfront=-0.5*fTotThickness;
  double xcenter_c=0.;
  double ycenter_c = 0;
  double zcenter = 0;
  double zcenter_layer = 0;

 // one cell
  vecgeom::UnplacedBox *fCellBox = new vecgeom::UnplacedBox("Cell", fCellThickness[i]/2,fCellSizeX/2,fCellSizeY/2);
  vecgeom::LogicalVolume *fCellLogic = new vecgeom::LogicalVolume("Cell", fCellBox);
  fCellLogic->SetMaterialPtr(fCellMaterial); 
  fCellLogic->SetRegion(calRegion);

//one cell row
  vecgeom::UnplacedBox *fRowBox = new vecgeom::UnplacedBox("CellRow", fCellThickness[i]/2,fCellSizeX/2,25*fCellSizeY/2);
  vecgeom::LogicalVolume *fRowLogic = new vecgeom::LogicalVolume("CellRow", fRowBox);

  for (int ix=0;ix<numCells;ix++) {
      sprintf(cellInRowName,"cellInRow%d",k);
      xcenter_c = fCellThickness*0.5 + ix*fCellXY;
      vecgeom::Transformation3D *cellPlaceInRow = new vecgeom::Transformation3D(xcenter_c,0,zfront, 0,0,0);
      fRowLogic->PlaceDaughter(cellInRowName,fCellBox,cellPlaceInRow);
  }

//one Layer
  vecgeom::UnplacedBox *fLayerBox = new vecgeom::UnplacedBox("CellLayer", fCellThickness[i]/2,25*fCellSizeX/2,25*fCellSizeY/2);
  vecgeom::LogicalVolume *fLayerLogic = new vecgeom::LogicalVolume("CellLayer", fLayerBox);

  for (int iy=0;iy<numRows;iy++) {
      sprintf(rowInLayerName,"abs%d",k);
      ycenter_c = fCellThickness*0.5 + iy*fCellXY;
      vecgeom::Transformation3D *rowPlace = new vecgeom::Transformation3D(0,ycenter_c,zfront, 0,0,0);
      fLayerLogic->PlaceDaughter(rowInLayerName,fRowLogic,rowPlace);
  }

//place layers in world
  for (int k=0; k<numLayers; k++) {
      sprintf(absName,"abs%d",k);
      sprintf(layerName,"abs%d",k);
      fAbsLogic[k] = new vecgeom::LogicalVolume(volName,fAbsBoxes[k]);
      fAbsLogic[k]->SetMaterialPtr(fAbsMaterial[k]);
      fAbsLogic[k]->SetRegion(calRegion);
      fAbsLogicVolumeID[k]=fAbsLogic[k]->id();
      
      zcenter = zfront+0.5*fAbsThickness[k] + fCellThickness;
      zcenter_layer = zfront+0.5*fCellThickness;
      zfront += ((k+1)*fCellThickness + fAbsThickness[k]);
      vecgeom::Transformation3D *absPlace = new vecgeom::Transformation3D(0,0,zcenter, 0,0,0);
      world->PlaceDaughter(volName,fAbsLogic[k],absPlace);
      vecgeom::Transformation3D *layerPlace = new vecgeom::Transformation3D(0,0,zcenter_layer, 0,0,0);
      world->PlaceDaughter(layerName,fLayerLogic,layerPlace);
  }
