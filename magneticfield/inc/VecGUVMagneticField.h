#ifndef VecGUVMagneticField_H
#define VecGUVMagneticField_H

#include "VecGUVField.h"


class VecGUVMagneticField :  public VecGUVField
{
  public:
    static constexpr int   fNumFieldComponents= 3;
    static constexpr bool  fFieldChangesEnergy= false;
  
    VecGUVMagneticField():  
     VecGUVField( fNumFieldComponents, fFieldChangesEnergy) 
    {
      // std::cout<<"--- VecGUVMagneticField entered here ---"<<std::endl;
    }

    virtual ~VecGUVMagneticField(){}; 

    template <class Backend>
    void  GetFieldValue( const typename Backend::precision_v  Point[4],     // The old interface
                               typename Backend::precision_v* Field );

    template <class Backend>
    void GetFieldValue( const vecgeom::Vector3D<typename Backend::precision_v> &Position, 
                              vecgeom::Vector3D<typename Backend::precision_v> &FieldValue );

    VecGUVMagneticField& operator = (const VecGUVMagneticField &p);
    //  Copy 'standard' components ...
};


template <class Backend>
void
VecGUVMagneticField::GetFieldValue( const typename Backend::precision_v  Point[4], // The old interface
                                          typename Backend::precision_v* FieldArr )
{
  typedef typename Backend::precision_v Double_v;

  vecgeom::Vector3D<Double_v> PositionV3D( Point[0], Point[1], Point[2]);
  vecgeom::Vector3D<Double_v>  Field_v3f;
  this->GetFieldValue( PositionV3D, Field_v3f );
  FieldArr[0]= Field_v3f.x();
  FieldArr[1]= Field_v3f.y();
  FieldArr[2]= Field_v3f.z();
}
#endif
