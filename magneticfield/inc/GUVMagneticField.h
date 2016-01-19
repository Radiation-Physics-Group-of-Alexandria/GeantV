#ifndef GUVMagneticField_H
#define GUVMagneticField_H

#include "GUVField.h"

class GUVMagneticField : public GUVField

{
  public:
    static constexpr int   fNumFieldComponents=3;
    static constexpr bool  fFieldChangesEnergy= false;
    GUVMagneticField():  GUVField( fNumFieldComponents, fFieldChangesEnergy) {}

    virtual ~GUVMagneticField(); // {}

    void  GetFieldValue( const double  Point[4],     // The old interface
                               double* Field );

    virtual void  GetFieldValue( const vecgeom::Vector3D<double> &Position, 
                                       vecgeom::Vector3D<float>  &FieldValue ) = 0;

    /*
     * The expected vector interface is: 
     *
     * virtual void GetFieldValue( const vecgeom::Vector3D<typename Backend::precision_v> &Position, 
     *                               vecgeom::Vector3D<typename Backend::precision_v> &FieldValue ) = 0;
     */

    // virtual GUVField* Clone() const;
    //   Concrete subclasses can (should?) implement it!

    GUVMagneticField& operator = (const GUVMagneticField &p);
    //  Copy 'standard' components ...
};
#endif
