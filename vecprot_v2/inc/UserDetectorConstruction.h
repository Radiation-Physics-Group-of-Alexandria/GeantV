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
#include "TUniformMagField.h"  //  For now - plan to strip it into separate 'simple' det constr.

class UserDetectorConstruction
{
  public:
    UserDetectorConstruction();
    virtual ~UserDetectorConstruction() {};
    // virtual bool Construct(const char *geomfile);
    virtual bool CreateFieldAndSolver(bool useRungeKutta= true);

    /** Register a constanct B-field */ 
    void UseConstantMagField( float value[3] );

    /** Register a B-field, and create integrator for it. */ 
    template <class Field_t>
    bool CreateSolverForField(Field_t* field);

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

  public:
    static constexpr double   fEpsilonDefault = 3.0e-5; 
    static constexpr double   fMinimumStepInFieldDef= 1.0e-4; // GV units = cm
    // vecgeom::Vector3D<float>  fMagField;
};

UserDetectorConstruction::UserDetectorConstruction() : 
   fEpsilonRK(fEpsilonDefault), 
   fMinimumStepInField(fMinimumStepInFieldDef),
   fUseUniformField(false),
   fZeroField(true),
   fCreatedField(false)
   {}

// bool UserDetectorConstruction::Construct(const char *geomfile)
// {
//    return LoadGeometry(geomfile));
// }

template <class Field_t>
bool UserDetectorConstruction::CreateSolverForField(Field_t* ptrField)
{
  FieldPropagatorFactory::CreatePropagator<Field_t>( *ptrField,
                                                     fEpsilonRK,
                                                     fMinimumStepInField);
  return true;
}

void UserDetectorConstruction::
UseConstantMagField( float fieldVal[3] ) // vecgeom::Vector3D<float> value )
{
  fUniformMagField= vecgeom::Vector3D<float>( fieldVal[0], fieldVal[1], fieldVal[2] );
  fUseUniformField= true;
  fZeroField = ( fUniformMagField.Mag2() == 0.0 );
}

/**  This method must exist for derived classes, 
  *    and must call CreateSolverForField() for concrete field class
  */
bool
UserDetectorConstruction::
CreateFieldAndSolver(bool useRungeKutta)
{
  static const char *method="UserDetectorConstruction::CreateFieldAndSolver";
  bool rtv= false;

  if( fUseUniformField && !fZeroField )
  {
    if( useRungeKutta )
    {
        auto gvField= new TUniformMagField( fieldUnits::kilogauss * fUniformMagField );
        rtv= CreateSolverForField<TUniformMagField>(gvField);
    }
    fCreatedField= true;
  }
  if (fZeroField) {
    Geant::Print(method," Zero Magnetic Field configured.");
    fCreatedField= true;
  }else{    
    Geant::Error(method,"No user Magnetic Field is registered.");
  }
  return rtv;
}
#endif