///////////////////////////////////////////////////////////////////////////////////////////
//
//			ClicDetectorConstruction.cxx
//			Created: September 2017
//			Author: S. Vallecorsa from examples by M Novak and R Schmitz
//
// Description: A (linear) calorimeter implemented using VecGeom libraries, this detector construction
// 		is fully customizable and can easily pass data about its structure to other
// 		classes.
// 		
// 		Customizable features include:
// 		-Number of calorimeter layers
// 		-Number of absorbers per layer
//		-Production cuts (on energy, length, or particle-specific for gammas, electrons, or positrons)
//		-YZ calorimeter cross-section
//
//		Other features like incident particle energy and particle type may be set in the
//		macro which executes this application, caloAppRP.cc, as well as in its accompanying
//		macro.
//
///////////////////////////////////////////////////////////////////////////////////////////

#include "ClicDetectorConstruction.h"
#include "GeantVDetectorConstruction.h"
#include <iostream>
#include <vector>

//Material includes
#include "Isotope.h"
#include "Element.h"
#include "Material.h"
#include "MaterialProperties.h"
#include "NISTElementData.h"
#include "SystemOfUnits.h"
#include "PhysicalConstants.h"

//Region and productionCut includes
#include "Region.h"
#include "PhysicsParameters.h"
#include "MaterialCuts.h"
/////////////////////////////////
//VecGeom includes
#include "management/GeoManager.h"
#include "volumes/Box.h"
#include "volumes/LogicalVolume.h"
/////////////////////////////////

namespace userapplication {

ClicDetectorConstruction::~ClicDetectorConstruction()
{
}

void ClicDetectorConstruction::CreateMaterials() {

  //unit definitions
  using geant::g;
  using geant::mg;
  using geant::mole;
  using geant::cm3;
//  using geant::bar;
  using geant::pascal;
  using geant::kelvin;
  using geant::atmosphere;
  using geant::kUniverseMeanDensity;
//  using geant::perCent;
//  using geant::kSTPTemperature;

  //useful variable declarations
  std::string name   = "";
  std::string symbol = "";
  double a       = 1.*g/mole,   // molar mass in internal [weight/mole] unit
         z       = 1.,          // mean numnber of protons
         density = 1.*g/cm3;    // material density in internal [weight/length^3] unit
//         abundance;             // relative abundance of the i-th isotope
//  int isoZ, isoN;               // number of protons/nucleons in an isotope;
//  int numcomponents,            // number of components the naterial is built up
//      numatoms;                 // number of i-th atoms in the molecule (for Way 1) some useful variable declaration

  double pressure    = 1.*atmosphere;  // pressure
  double temperature = 273.15*kelvin; // temperature

  //Isotope from geantphysics namespace
  using geantphysics::Isotope;
  // Element from geantphysics namespace
  using geantphysics::Element;
  // Material from the geantphysics namespace
  using geantphysics::Material;
  // MaterialState from geantphysics namespace
  using geantphysics::MaterialState;
  using geantphysics::NISTElementData;
  using geantphysics::MaterialProperties;
  
    //Define vacuum
    density     = kUniverseMeanDensity;
    pressure    = 3.e-18*pascal;
    temperature = 2.73*kelvin;
    Material *Galactic = new Material(name="Galactic", z=1., a=1.01*g/mole, density,
                 MaterialState::kStateGas,temperature,pressure);
  
    //Explicit Pb definition
 //   density = 11.35*g/cm3;
//    a = 207.19*g/mole;
//    Material* Pb = new Material(name="Lead"     , z=82., a, density);

    //Declare calo materials
    fWorldMaterialName = "Galactic";
    fWorldMaterial     =  Galactic ;
}
/////////////////////////////////////////////////////////////////////////////////////////////////////////

void ClicDetectorConstruction::SetProductionCutsByEnergy(double energy){

  prodCutByLength=false; 
  fGammaCut=energy*geant::eV;
  fElectronCut=energy*geant::eV;
  fPositronCut=energy*geant::eV;
}

void ClicDetectorConstruction::SetProductionCutsByLength(double length){
  fGammaCut=length*geant::mm;
  fElectronCut=length*geant::mm;
  fPositronCut=length*geant::mm;
}

//all these settings are in length
void ClicDetectorConstruction::SetDetectorGammaProductionCut(double length){  fGammaCut=length*geant::mm;}
void ClicDetectorConstruction::SetDetectorElectronProductionCut(double length){  fElectronCut=length*geant::mm;}
void ClicDetectorConstruction::SetDetectorPositronProductionCut(double length){  fPositronCut=length*geant::mm;}
/////////////////////////////////////////////////////////////////////////////////////////////////////////

void ClicDetectorConstruction::SetAbsorberMaterialName(int absNum, std::string matName){	
	if (absNum<=numAbsorbers){
		userMaterial[absNum]=true;
		fAbsMaterialName[absNum]=matName;
	} else {std::cerr << "ERROR: There are too few absorbers! Increase the number of absorbers or set a different absorber's material.\n";}
}
std::string ClicDetectorConstruction::GetAbsorberMaterialName(int absNum){return fAbsMaterialName[absNum];}
geantphysics::Material* ClicDetectorConstruction::GetAbsorberMaterial(int absNum){
	if (absNum<=numAbsorbers){
		return fAbsMaterial[absNum];
	} else {std::cerr << "ERROR: There are too few absorbers! Increase the number of absorbers or set a different absorber's material.\n"; 
	return nullptr;
	}
}
void ClicDetectorConstruction::SetDetectorMaterials(){ //private, set each absorber to its material name
	for (int i=0; i<numAbsorbers; i++){
		fAbsMaterial[i]=geantphysics::Material::NISTMaterial(fAbsMaterialName[i]);
	}
	fWorldMaterial=geantphysics::Material::NISTMaterial(fWorldMaterialName);
	fCellMaterial=geantphysics::Material::NISTMaterial(fCellMaterialName);
}

/////////////////////////////////////////////////////////////////////////////////////////////////////////

void ClicDetectorConstruction::SetNumLayers(int nLayers){
	userLayerNum=true;
	numLayers=nLayers;
}
void ClicDetectorConstruction::SetNumAbsorbers(int nAbsorbers){
	userAbsorberNum=true;
	numAbsorbers=nAbsorbers;
}
int ClicDetectorConstruction::GetNumLayers() {return numLayers;}
int ClicDetectorConstruction::GetNumAbsorbers() {return numAbsorbers;}

/////////////////////////////////////////////////////////////////////////////////////////////////////////

void ClicDetectorConstruction::SetAbsorberThickness(int absNum, double thickness){
	if (absNum<=numAbsorbers){
		userThickness[absNum]=true;
		fAbsThickness[absNum]=thickness;
	} else {std::cerr << "ERROR: There are too few absorbers! Increase the number of absorbers or set a different absorber's material.\n"; }
}
double ClicDetectorConstruction::GetAbsorberThickness(int absNum){
	if (absNum<=numAbsorbers){
		return fAbsThickness[absNum];
	} else {std::cerr << "ERROR: There are too few absorbers! Increase the number of absorbers or set a different absorber's material.\n";
	return 0;
	}
}

/////////////////////////////////////////////////////////////////////////////////////////////////////////

void ClicDetectorConstruction::SetDetectorXY(double xy){
	userCaloXY=true;
	fCaloSizeXY=xy;
}

double ClicDetectorConstruction::GetDetectorZ() {return fCellThickness*numLayers;}
double ClicDetectorConstruction::GetDetectorXY() {return fCaloSizeXY;}
double ClicDetectorConstruction::GetWorldZ() {return fWorldSizeZ;}
double ClicDetectorConstruction::GetWorldXY() {return fWorldSizeXY;}

/////////////////////////////////////////////////////////////////////////////////////////////////////////

int ClicDetectorConstruction::GetAbsorberLogicalVolumeID(int absorber){ return fAbsLogicVolumeID[absorber]; }
int ClicDetectorConstruction::GetLayerLogicalVolumeID(int layer){ return layer; }
int ClicDetectorConstruction::GetDetectorRegionIndex() {return fDetectorRegionIndex; }

/////////////////////////////////////////////////////////////////////////////////////////////////////////

void ClicDetectorConstruction::CreateGeometry() {

  using geant::mm;
  using geant::cm;

  //create default detector parameters if the user hasn't defined them
  numCells = 25;
  numRows = 25;
  numAbsorbers=numLayers = 25;
  numAbsorbers1 = 17;
  numAbsorbers2 = 8;
  
  fCellThickness=0.5*mm;
  fCellX = 5.1*mm;
  fCellY = 5.1*mm;

  fCellMaterialName = "NIST_MAT_Si";

  for (int i = 0; i<numAbsorbers1;i++) fAbsThickness[i] = 2.4*mm;
  for (int i = numAbsorbers1; i<=(numAbsorbers1+numAbsorbers2);i++) fAbsThickness[i] = 4.8*mm;
  for (int i = 0; i<maxAbsorbers;i++) fAbsMaterialName[i]="NIST_MAT_W";

  for (int p=0;p<maxAbsorbers;p++)
	std::cout << "Absorber number " << p << " has material " << fAbsMaterialName[p] << " and thickness " << fAbsThickness[p] << "cm"<< std::endl;
  
  fCaloSizeXY = 25*fCellX;

  fTotThickness = 25*fCellThickness;
  for (int i=1;i<=maxAbsorbers;i++){fTotThickness+=fAbsThickness[i];}
  fWorldSizeZ=1.2*fTotThickness ;
  fWorldSizeXY=1.2*fCaloSizeXY;
  SetDetectorMaterials();  
  std::cout<<"SetDetectorMaterials OK"<<std::endl;
  
//define regions
//for now both regions have the same cuts; can change this if a user wants, but it's not very relevant for this example
  vecgeom::Region *worldRegion = new vecgeom::Region("WorldRegion",prodCutByLength,fGammaCut,fElectronCut,fPositronCut);
  vecgeom::Region *calRegion = new vecgeom::Region("CalRegion",prodCutByLength,fGammaCut,fElectronCut,fPositronCut);
  fDetectorRegionIndex = calRegion->GetIndex();



 //define world
  vecgeom::UnplacedBox *worldUnplaced      = new vecgeom::UnplacedBox(fWorldSizeXY/2, fWorldSizeXY/2, fWorldSizeZ/2); 
  vecgeom::LogicalVolume *world = new vecgeom::LogicalVolume("world", worldUnplaced);
  world->SetRegion(worldRegion);
  world->SetMaterialPtr(fWorldMaterial); 

// absorbers
  vecgeom::UnplacedBox *fAbsBoxes[maxAbsorbers];   // two different thicknesses: x17 (2.4 mm) x8 (4.8 mm)
  for (int i=1; i<=maxAbsorbers; i++) {   
      fAbsBoxes[i] = new vecgeom::UnplacedBox("Absorber", fCaloSizeXY/2,fCaloSizeXY/2, fAbsThickness[i]/2);
  }
  vecgeom::LogicalVolume *fAbsLogic[maxAbsorbers];

  char *absName = new char;
  char *cellInRowName = new char;
  char *rowInLayerName = new char;
  char *layerName = new char;

  double zfront=-0.5*fTotThickness;
  double xcenter_c=0.;
  double ycenter_c = 0;
  double zcenter = 0;
  double zcenter_layer = 0;
  std::cout<<"Logic absorbers OK"<<std::endl;
//one cell
  vecgeom::UnplacedBox *fCellBox = new vecgeom::UnplacedBox("Cell",fCellX/2,fCellY/2,fCellThickness/2);
  vecgeom::LogicalVolume *fCellLogic = new vecgeom::LogicalVolume("Cell", fCellBox);
  fCellLogic->SetMaterialPtr(fCellMaterial); 
  fCellLogic->SetRegion(calRegion);
  std::cout<<"Logic cells OK"<<std::endl;

//one cell row
  vecgeom::UnplacedBox *fRowBox = new vecgeom::UnplacedBox("CellRow",fCellX/2,25*fCellY/2, fCellThickness/2);
  vecgeom::LogicalVolume *fRowLogic = new vecgeom::LogicalVolume("CellRow", fRowBox);
  fRowLogic->SetMaterialPtr(fCellMaterial); 
  fRowLogic->SetRegion(calRegion);

  std::cout<<"Logic rows OK"<<std::endl;

  for (int ix=0;ix<numCells;ix++) {
      sprintf(cellInRowName,"cellInRow%d",ix);
      xcenter_c = fCellThickness*0.5 + ix*fCellX;
      vecgeom::Transformation3D *cellPlaceInRow = new vecgeom::Transformation3D(xcenter_c,0,zfront, 0,0,0);
      fRowLogic->PlaceDaughter(cellInRowName,fCellLogic,cellPlaceInRow);
  }
  std::cout<<"Cells in rows OK"<<std::endl;
//one Layer
  vecgeom::UnplacedBox *fLayerBox = new vecgeom::UnplacedBox("CellLayer", 25*fCellX/2,25*fCellY/2, fCellThickness/2);
  vecgeom::LogicalVolume *fLayerLogic = new vecgeom::LogicalVolume("CellLayer", fLayerBox);
  fLayerLogic->SetMaterialPtr(fCellMaterial); 
  fLayerLogic->SetRegion(calRegion);

  for (int iy=0;iy<numRows;iy++) {
      sprintf(rowInLayerName,"row%d",iy);
      ycenter_c = fCellThickness*0.5 + iy*fCellY;
      vecgeom::Transformation3D *rowPlace = new vecgeom::Transformation3D(0,ycenter_c,zfront, 0,0,0);
      fLayerLogic->PlaceDaughter(rowInLayerName,fRowLogic,rowPlace);
  }

//place layers in world
  for (int k=0; k<numLayers; k++) {
      sprintf(absName,"abs%d",k);
      sprintf(layerName,"layer%d",k);
      fAbsLogic[k] = new vecgeom::LogicalVolume(absName,fAbsBoxes[k]);
      fAbsLogic[k]->SetMaterialPtr(fAbsMaterial[k]);
      fAbsLogic[k]->SetRegion(calRegion);
      fAbsLogicVolumeID[k]=fAbsLogic[k]->id();

      zcenter = zfront+0.5*fAbsThickness[k] + fCellThickness;
      zcenter_layer = zfront+0.5*fCellThickness;
      zfront += ((k+1)*fCellThickness + fAbsThickness[k]);
      vecgeom::Transformation3D *absPlace = new vecgeom::Transformation3D(0,0,zcenter, 0,0,0);
      world->PlaceDaughter(absName,fAbsLogic[k],absPlace);
      vecgeom::Transformation3D *layerPlace = new vecgeom::Transformation3D(0,0,zcenter_layer, 0,0,0);
      world->PlaceDaughter(layerName,fLayerLogic,layerPlace);
  }


//place world volume, close geometry
  world->PrintContent();
  vecgeom::VPlacedVolume *w = world->Place();
  vecgeom::GeoManager::Instance().SetWorld(w);
  vecgeom::GeoManager::Instance().CloseGeometry();
}

}
