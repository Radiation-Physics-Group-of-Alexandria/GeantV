//===----------------------------------------------------------------------===//
/**
 * @file GUVVectorField.h
 * @brief  Abstract field class for Geant-V prototype
 * @author John Apostolakis
 */
//===----------------------------------------------------------------------===//

//
// class GUVVectorField
//
// Class description:
//
// Abstract class for any kind of Field.
// It allows a vector field to be defined by implementing the inquiry function
// interface.
//
// The key method is GetFieldValue( vecgeom::Vector3D<double> const &Position,
//                   *************  vecgeom::Vector3D<float>        &FieldValue)
//                                  
// Given an input position 3-vector 'Position',
// this method must return the estimated value of the field in "FieldValue".
//
// A field must also specify whether it changes a track's energy:
//                    DoesFieldChangeEnergy()
//                    *********************
// A field must co-work with a corresponding Equation of Motion, to
// enable the integration of a particle's position, momentum and, 
// optionally, polarisation.
// -------------------------------------------------------------------

#ifndef GUVVECTORFIELD_HH
#define GUVVECTORFIELD_HH

// #include <vector>

#include "base/Global.h"     //  Defines Force_Inline, VectorBackend, .. 
#include "backend/Backend.h"
#include "base/Vector3D.h"
// #include "base/SOA3D.h"

// #include "GUVField.h"

/*
namespace vecfield {

VECGEOM_DEVICE_FORWARD_DECLARE(class GUVVectorField;);
VECGEOM_DEVICE_DECLARE_CONV(class, GUVVectorField);

inline namespace VECFIELD_IMPL_NAMESPACE {
*/

// using VectorBackend = vecCore::backend::VcVector;
// using VectorType = VectorBackend::double;

class GUVVectorField //  : public GUVField
{
   public:
      using Double_v = vecgeom::VectorBackend::Double_v;
      using Float_v  = vecgeom::VectorBackend::Float_v;
      // using Vector3D = vecgeom::VectorBackend::Float_v;
   
      inline
      GUVVectorField( int NumberOfComponents, bool changesEnergy );
      inline
      GUVVectorField( const GUVVectorField &);
      virtual ~GUVVectorField();

      VECGEOM_CUDA_HEADER_BOTH
      virtual void GetFieldValue( vecgeom::Vector3D<double> const &Position,
                                  vecgeom::Vector3D<float>        &FieldValue ) = 0;

      //Vector interface - chosen backend type (defined globally)
      virtual void GetFieldValueSIMD( vecgeom::Vector3D<Double_v> const &Position,
                                      vecgeom::Vector3D<Float_v>      &FieldValue ) = 0;

      // a helper tramponline to dispatch to DistanceToOut if type is not scalar
      template <typename S, typename T>
      VECGEOM_FORCE_INLINE
      VECGEOM_CUDA_HEADER_BOTH
      void GetFieldValue(vecgeom::Vector3D<S> const &position, vecgeom::Vector3D<T> &fieldValue)
      {
         return GetFieldValueSIMD(position, fieldValue);
      }

      VECGEOM_CUDA_HEADER_BOTH
      bool DoesFieldChangeEnergy() const { return fChangesEnergy; }
      VECGEOM_CUDA_HEADER_BOTH
      int  GetNumberOfComponents() const { return fNumberOfComponents; } 

      GUVVectorField& operator = (const GUVVectorField &p); // Useful ?
      
      virtual GUVVectorField* Clone() const;

  private:
      const int  fNumberOfComponents; 
      bool       fChangesEnergy; 
};

inline GUVVectorField::GUVVectorField( int numberOfComponents, bool changesEnergy )
   : fNumberOfComponents(numberOfComponents),
     fChangesEnergy(changesEnergy)
     //GUVField(numberOfComponents, changesEnergy)
{
  std::cout<<"-- entered GUVVectorField  constructor--"<<std::endl;
}


inline GUVVectorField::GUVVectorField( const GUVVectorField &field) 
  :  fNumberOfComponents(field.fNumberOfComponents)
    //GUVField(field)
{
  fChangesEnergy= field.fChangesEnergy;
}

/*
} // End inline namespace

} // End global namespace
*/

#endif /* GUVFVECTORIELD_HH */
