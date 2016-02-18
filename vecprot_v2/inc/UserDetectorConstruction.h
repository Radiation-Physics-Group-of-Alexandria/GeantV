//===--- UserDetectorConstruction.h - Geant-V ----------------------------*- C++ -*-===//
//
//                     Geant-V Prototype
//
//===----------------------------------------------------------------------===//
/**
 * @file   UserDetectorConstruction.h
 * @brief  Base class for the user's mandatory initialization class
 *         for detector setup.
 * @author John Apostolakis
 */
//===----------------------------------------------------------------------===//

#ifndef USERDETECTORCONSTRUCTION_H
#define USERDETECTORCONSTRUCTION_H 1

#include "base/Vector3D.h"
// #include "Geant/Typedefs.h"

#include "FieldPropagatorFactory.h"
#include "Units.h"  //  For fieldUnits

#include "TUniformMagField.h"  //  For now - plan to strip it into separate 'simple' det constr.

class UserDetectorConstruction
{
  public:
    // UserDetectorConstruction();  // RootDL
    virtual ~UserDetectorConstruction() {};
    // virtual bool Construct(const char *geomfile);
    // virtual bool CreateFieldAndSolver(bool useRungeKutta= true); // RootDL

    /** Register a constanct B-field */ 
    // void UseConstantMagField( float value[3], const char* Unit= 0 ); // Default unit is kilogauss // RootDL

    /** Register a B-field, and create integrator for it. */ 
    // template <class Field_t>                   // RootDL
    // bool CreateSolverForField(Field_t* field); // RootDL

    void SetEpsilonRK(double val) { fEpsilonRK = val; }
    void SetMinimumStep( double v) { fMinimumStepInField = v; }

    // Inquiry methods
    bool   IsFieldUniform() { return fUseUniformField; }
    double GetEpsilonRK() { return fEpsilonRK; }
    double GetMinimumStep() { return fMinimumStepInField; }
    bool   IsFieldCreated() { return fCreatedField; }
  
  private:
    double           fEpsilonRK;
    double           fMinimumStepInField;
    vecgeom::Vector3D<float>  fUniformMagField;
    bool             fUseUniformField;
    bool             fZeroField;

  protected: 
    bool             fCreatedField;
    bool             fCalled;

  public:
    static constexpr double   fEpsilonDefault = 3.0e-5; 
    static constexpr double   fMinimumStepInFieldDef= 1.0e-4; // GV units = cm
    // vecgeom::Vector3D<float>  fUniformMagFieldValue;

// };   // RootComm

// --> Changed to accomodate Root needs for 
public: // RootAdded
   
// UserDetectorConstruction:: // RootComm
UserDetectorConstruction() : 
   fEpsilonRK(fEpsilonDefault), 
   fMinimumStepInField(fMinimumStepInFieldDef),
   fUseUniformField(false),
   fZeroField(true),
   fCreatedField(false),
   fCalled(false)
   {}

// bool UserDetectorConstruction::Construct(const char *geomfile)
// {
//    return LoadGeometry(geomfile));
// }

template <class Field_t>
bool
// UserDetectorConstruction:: // RootComm
CreateSolverForField(Field_t* ptrField)
{
  FieldPropagatorFactory::CreatePropagator<Field_t>( *ptrField,
                                                     fEpsilonRK,
                                                     fMinimumStepInField);
  fCreatedField= true;  
  return true;
}

void
// UserDetectorConstruction:: // RootComm
UseConstantMagField( float fieldVal[3],  const char* Units =0 )
{
  const char *methodName= "UserDetectorConstruction::UseConstantMagField";
  bool defaultUsed= false;
  double unit= 1;
  
  if( Units == 0  || strcmp(Units,"kilogauss") == 0 ) {
    unit= fieldUnits::kilogauss;
    defaultUsed = (Units == 0);
  } else if( ( strcmp(Units,"gauss") == 0 ) || ( strcmp(Units,"Gauss") == 0 ) ) {
    unit= fieldUnits::gauss;
  } else if( ( strcmp(Units,"tesla") == 0 ) || ( strcmp(Units,"Tesla") == 0 ) ) {
    unit= fieldUnits::gauss;
  } else {
    unit= fieldUnits::kilogauss;
    defaultUsed = (Units == 0);     
  }

  if( defaultUsed )
     printf("%s - WARNING: No units provided - using kilogauss as default unit", 
            "UserDetectorConstruction::UseConstantMagField");

  fUniformMagField= vecgeom::Vector3D<float>( fieldVal[0] * unit, fieldVal[1] * unit, fieldVal[2] * unit );
  
  printf("%s - Info: UseConstantMagField called. Field value = %9.3g , %9.3g  %9.3g  kiloGauss\n",
         methodName,
         fUniformMagField[0] / fieldUnits::kilogauss,       fUniformMagField[1] / fieldUnits::kilogauss,
         fUniformMagField[2] / fieldUnits::kilogauss );

  fUseUniformField= true;
  fZeroField = ( fUniformMagField.Mag2() == 0.0 );
}

/**  This method must exist for derived classes, 
  *    and must call CreateSolverForField() for concrete field class
  */
bool
// UserDetectorConstruction:: // RootComm
CreateFieldAndSolver(bool /*useRungeKutta*/ )
{
  static const char *method="UserDetectorConstruction::CreateFieldAndSolver";
  bool rtv= false;
   
  // Geant::Print(method, "%s - method called.  Uniform= %d  Zero-value= %d.", method, fUseUniformField, fZeroField );

  if( fUseUniformField )
  {
    auto gvField= new TUniformMagField( fUniformMagField );
    rtv= CreateSolverForField<TUniformMagField>(gvField);

    if (fZeroField) {
      Geant::Print(method," Zero Magnetic Field configured.");
    } else {
      printf("Creating uniform B-field: %6.2g  %6.2g  %6.2g kilo-gauss\n",
           fUniformMagField[0]/fieldUnits::kilogauss,
           fUniformMagField[1]/fieldUnits::kilogauss,
           fUniformMagField[2]/fieldUnits::kilogauss );
      printf("  Recall  kilogauss = %8.2g in field Units\n", fieldUnits::kilogauss );       
    }
    fCalled = true;
  } else {
    fCalled = true;
    Geant::Error(method,"No user Magnetic Field is registered.");
  }
  return rtv;
}

}; // RootAdded   
#endif
