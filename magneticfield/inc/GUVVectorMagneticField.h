#ifndef GUVVectorMagneticField_H
#define GUVVectorMagneticField_H

#include "GUVVectorField.h"
#include "GUVField.h"
#include "GUVMagneticField.h"

/*
namespace vecfield {

VECGEOM_DEVICE_FORWARD_DECLARE(class GUVVectorField;);
VECGEOM_DEVICE_DECLARE_CONV(class, GUVVectorField);

inline namespace VECFIELD_IMPL_NAMESPACE {
*/

class GUVVectorMagneticField :  public GUVVectorField 
{
    // typedef typename vecgeom::kVc::precision_v      Double_v;
    // typedef typename vecgeom::kVcFloat::precision_v Float_v;
  
  public:
    static constexpr int   sNumFieldComponents= 3;
    static constexpr bool  sFieldChangesEnergy= false;
  
    GUVVectorMagneticField():
     GUVVectorField( sNumFieldComponents, sFieldChangesEnergy) 
    {std::cout<<"--- GUVVectorMagneticField def c-tor called ---"<<std::endl;}

    virtual ~GUVVectorMagneticField(){}; 

    /***
    void  GetFieldValue( const Double_v  Point[4],     // The old interface
                               Double_v* Field );
    ****/
    VECGEOM_CUDA_HEADER_BOTH
    virtual void GetFieldValue( vecgeom::Vector3D<double> const &Position,
                                vecgeom::Vector3D<float>        &FieldValue ) = 0;

    virtual void GetFieldValueSIMD( vecgeom::Vector3D<Double_v>  const &Position, 
                                    vecgeom::Vector3D<Float_v>         &FieldValue ) = 0;

    GUVVectorMagneticField& operator = (const GUVVectorMagneticField &p);
    //  Copy 'standard' components ...
};

/*******
void
GUVVectorMagneticField::GetFieldValue( const Double_v  Point[4],     // The old interface
                                             Double_v* FieldArr )
{
   // typedef typename vecgeom::kVc::precision_v Double_v;
   // typedef typename vecgeom::kVcFloat::precision_v Float_v;
   
   vecgeom::Vector3D<Double_v> PositionV3D( Point[0], Point[1], Point[2]);
   vecgeom::Vector3D<Float_v>  Field_v3f;
   this->GetFieldValue( PositionV3D, Field_v3f );
   FieldArr[0]= (Double_v) Field_v3f.x();
   FieldArr[1]= (Double_v) Field_v3f.y();
   FieldArr[2]= (Double_v) Field_v3f.z();
}
 ******/

/*
} // End inline namespace

} // End global namespace
*/

#endif
