#include <string>
#include "Geant/Error.h"
#include "CMSDetectorConstruction.h"
#include "CMSmagField.h"

CMSDetectorConstruction::CMSDetectorConstruction() :
  fFieldFilename(std::string("")),
  fCMSfield(0)
{}

CMSDetectorConstruction::CMSDetectorConstruction(const char* fieldFilename) :
  fFieldFilename(fieldFilename),
  fCMSfield(0)
{}

CMSDetectorConstruction::CMSDetectorConstruction(std::string fieldFilename) :
  fFieldFilename(fieldFilename),
  fCMSfield(0)
{}

CMSDetectorConstruction::~CMSDetectorConstruction()
{
  delete fCMSfield;
}   

// ClassImp(CMSDetectorConstruction);

bool
CMSDetectorConstruction::
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
