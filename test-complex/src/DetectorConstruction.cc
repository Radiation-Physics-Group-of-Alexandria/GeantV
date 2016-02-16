
#include "DetectorConstruction.hh"
#include "TGeoManager.h"
#include "TabulatedDataManager.hh"
#include "MaterialConverter.hh"

#include "G4UniformMagField.hh"
#include "G4FieldManager.hh"
#include "G4ChordFinder.hh"
#include "G4TransportationManager.hh"

#include "CMSMagneticFieldG4.hh"

TGeoManager * DetectorConstruction::fgGeomMgrRoot= 0 ; // Pointer to the geometry manager   
CMSMagneticFieldG4* DetectorConstruction::fpMasterCmsField= 0;

DetectorConstruction::DetectorConstruction(G4VPhysicalVolume *setWorld, bool uniformMagField ) :
   fWorld(setWorld), fUseUniformField( uniformMagField ), fpMagField(0),
   // fDistanceConst( 2.0 * CLHEP::millimeter ),
   fFieldFilename("cms2015.txt")
{   
   // fWorld = setWorld;
/*      char* gdmlFileName = getenv("VP_GEOM_GDML");
      if( gdmlFileName ){
        std::cout << " Creating empty TGeoManager by reading Root geometry from file " << gdmlFileName  << G4endl;
        fgGeomMgrRoot = TGeoManager::Import(gdmlFileName);
      } else {
        std::cout << " Creating empty TGeoManager " << std::endl;
        fgGeomMgrRoot = new TGeoManager();
      }
*/
//      std::cout << " Creating empty TGeoManager " << std::endl;
  fgGeomMgrRoot = new TGeoManager();

  TabulatedDataManager::SetTGeomManager( fgGeomMgrRoot );
  MaterialConverter::SetTGeomManager( fgGeomMgrRoot );
  MaterialConverter::Instance();// call just to initialize

  if( fUseUniformField ) { 
     // initialize magnetic field :: same value as in the prototype
     SetUniformBzMagField( fBzFieldValue );
  } else { 
     // To use the CMS field map - must create it from the field-dump file
     if( ! fpMasterCmsField ) 
        fpMasterCmsField= new CMSMagneticFieldG4( fFieldFilename );
     fpMagField= fpMasterCmsField;     
  }
}

DetectorConstruction::~DetectorConstruction()
{
  if( fpMagField != fpMasterCmsField ) delete fpMagField;
  delete fpMasterCmsField;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void DetectorConstruction::SetUniformBzMagField(G4double BzValue)
{
  //apply a global uniform magnetic field along Z axis
  G4FieldManager* fieldMgr
   = G4TransportationManager::GetTransportationManager()->GetFieldManager();

//  if(fpMagField) delete fpMagField;                //delete the existing magn field

  fBzFieldValue=    BzValue;
  fUseUniformField= true;
  
  delete fpMagField;
  delete fieldMgr->GetChordFinder();
  
  if(BzValue!=0.)                        // create a new one if non nul
  { fpMagField = new G4UniformMagField(G4ThreeVector(0.,0.,BzValue));
    fieldMgr->SetDetectorField(fpMagField);
    fieldMgr->CreateChordFinder(fpMagField);
  } else {
    fpMagField = 0;
    fieldMgr->SetDetectorField(fpMagField);
  }

  /***   The new way of creating it -- following examples/basic 2015
  G4ThreeVector fieldVector = G4ThreeVector(0.,0.,BzValue));
  fMagFieldMessenger = new G4GlobalMagFieldMessenger(BzValue);
  fMagFieldMessenger->SetVerboseLevel(1);  
  ***/
}

// Called in each thread to create a magnetic field
//
void DetectorConstruction::ConstructSDandField()
{
  // Sensitive detectors

   
  // Create global magnetic field messenger.
  // Uniform magnetic field is then created automatically if
  // the field value is not zero.
  if( fUseUniformField ){
     SetUniformBzMagField( fBzFieldValue );
  } else {

     if( fpMagField != fpMasterCmsField ) delete fpMagField;
     if( ! fpMasterCmsField ) 
        fpMasterCmsField= new CMSMagneticFieldG4( fFieldFilename );

     fpMagField= fpMasterCmsField; // ->CloneOrSafeSelf();
      // Re-use the same object in all threads -- methods must remain 'truly' const

     // fpMagField= G4CachedMagneticField( fpMasterCmsField, fDistanceConst ); 
     
     G4FieldManager* fieldMgr
        = G4TransportationManager::GetTransportationManager()->GetFieldManager();
     delete fieldMgr->GetChordFinder();

     fieldMgr->SetDetectorField(fpMagField);
     fieldMgr->CreateChordFinder(fpMagField);
  }

  // G4AutoDelete::Register(fMagFieldMessenger);
}
