#ifndef GEANT_CMS_Detector_Construction
#define GEANT_CMS_Detector_Construction

#include "UserDetectorConstruction.h"
#include "CMSDetectorConstruction.h"

#ifdef  USE_ROOT_TObject
#ifndef ROOT_TObject
#include "TObject.h"
#endif
#endif

#include <string>
#include "Geant/Error.h"
#include "CMSDetectorConstruction.h"
#include "CMSmagField.h"

class CMSmagField;

class CMSDetectorConstruction : public UserDetectorConstruction 
#ifdef  USE_ROOT_TObject
                               , public TObject
#endif
{
  public:
    /** @brief Destructor */
    CMSDetectorConstruction():  fFieldFilename(std::string("")), fCMSfield(0) {}
    // CMSDetectorConstruction(const char* fieldFilename);
    // CMSDetectorConstruction(std::string fieldFilename);
    ~CMSDetectorConstruction() { delete fCMSfield; } ;
    
    /** @brief Destructor */
    void SetFileForField(const char *filename){ fFieldFilename= filename; }
    void SetFileForField(std::string filename){ fFieldFilename= filename; }

    /** @brief Method to register a B-field, and create integrator for it. */
    // bool CreateFieldAndSolver(bool useRungeKutta= true);  // override final;

  private:
    std::string   fFieldFilename;
    CMSmagField*  fCMSfield;
    TUniformMagField*  fUniformField; // Alternative - for debugging only
    /** Field is created and owned by this class */

    // ClassDef(CMSDetectorConstruction, 1) // User application

// };

//  Implementations made inline, in order to cope with need to load dynamically,
//   using ROOT v6.
public:

// CMSDetectorConstruction::
CMSDetectorConstruction(const char* fieldFilename) :
  fFieldFilename(fieldFilename),
  fCMSfield(0)
{}

// CMSDetectorConstruction::
CMSDetectorConstruction(std::string fieldFilename) :
  fFieldFilename(fieldFilename),
  fCMSfield(0)
{}

// ClassImp(CMSDetectorConstruction);

bool
// CMSDetectorConstruction::
CreateFieldAndSolver(bool useRungeKutta)

{
  // static std::string OnOff[]= { "Off", "On" };
  Geant::Print("CMSDetectorConstruction::CreateFieldAndSolver", " Called with Arg: useRungeKutta=");
   //   (useRungeKutta ? "on" : "Off" ) ); 

  if(useRungeKutta )   { printf("on"); }  else { printf("Off"); }

#if 1
  using FieldType = CMSmagField;  
  std::cout << "    Calling CMSmagField constructor with filename= " << fFieldFilename << std::endl;
  fCMSfield= new CMSmagField(fFieldFilename);
  fUniformField= nullptr;
  
  auto fieldPtr = fCMSfield;  
#else
  using FieldType = TUniformMagField;
  double BzValue = 38.1 * fieldUnits::kilogauss;
  printf("    Called CMSmagField c-tor to create Uniform field with value %f \n", BzValue );
  fCMSfield= nullptr;

  vecgeom::Vector3D<double>  fieldValue( 0.0, 0.0, BzValue ); 
  fUniformField= new TUniformMagField( fieldValue);

  auto fieldPtr = fUniformField;
#endif
  
  Geant::Print("CMSDetectorConstruction::CreateFieldAndSolver", "CMSmagfield created.");
  
  if(useRungeKutta){
    CreateSolverForField<FieldType>(fieldPtr);
    printf("%s", "CMSdetectorConstruction - Configured field propagation with Runge Kutta.");    
  } else {
    printf("%s", "CMSdetectorConstruction - NOT configuring field propagation with Runge Kutta.");
  }  

  return true;
}

};  // Added for ROOT
#endif
