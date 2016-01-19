//===----------------------------------------------------------------------===//
/**
 * @file GUVVectorField.h
 * @brief  Abstract field class for Geant-V prototype
 * @author John Apostolakis
 */
//===----------------------------------------------------------------------===//

//
//
// class GUVVectorField
//
// Class description:
//
// Abstract class for any kind of Field.
// It allows any kind of field (vector, scalar, tensor and any set of them)
// to be defined by implementing the inquiry function interface.
//
// The key method is  GetFieldValue( const  double Point[4],
//                    *************         double *fieldArr ) 
// Given an input position/time vector 'Point', 
// this method must return the value of the field in "fieldArr".
//
// A field must also specify whether it changes a track's energy:
//                    DoesFieldChangeEnergy() 
//                    *********************
// A field must co-work with a corresponding Equation of Motion, to
// enable the integration of a particle's position, momentum and, optionally, 
// spin.  For this a field and its equation of motion must follow the
// same convention for the order of field components in the array "fieldArr"
// -------------------------------------------------------------------

#ifndef GUVVECTORFIELD_HH
#define GUVVECTORFIELD_HH

#include <vector>
#include "base/Vector3D.h"
#include "base/SOA3D.h"
#include "base/Global.h"
#include "backend/Backend.h"

#include "GUVField.h"


class GUVVectorField : public GUVField
{
  public: 

      //Vector interface with specialization
      virtual void GetFieldValue( const vecgeom::Vector3D<typename vecgeom::kVc::precision_v>  &Position,
                                        vecgeom::Vector3D<typename vecgeom::kVcFloat::precision_v>  &FieldValue ) = 0;

      inline
      GUVVectorField( int NumberOfComponents, bool changesEnergy );
      inline
      GUVVectorField( const GUVVectorField &);
      virtual ~GUVVectorField();

      GUVVectorField& operator = (const GUVVectorField &p); // Useful ?
      
      virtual GUVVectorField* Clone() const;
};

inline GUVVectorField::GUVVectorField( int numberOfComponents, bool changesEnergy )
   : GUVField(numberOfComponents, changesEnergy)
{
}


inline GUVVectorField::GUVVectorField( const GUVVectorField &field) 
  : GUVField(field)
{
}
#endif /* GUVFVECTORIELD_HH */
