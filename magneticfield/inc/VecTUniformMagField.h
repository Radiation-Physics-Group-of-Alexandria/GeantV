//
//  First version:      (Josh) - GSoC 2014 project
//  Current version:  J. Apostolakis

#ifndef VecTUniformMagField_H
#define VecTUniformMagField_H

#include "VecGUVMagneticField.h"
#include <iostream>

#include "base/Vector3D.h"

#include "Constants.h"  


class VecTUniformMagField : public VecGUVMagneticField
{
  public:  // with description

    VecTUniformMagField(const vecgeom::Vector3D<double>& FieldVector )
       : VecGUVMagneticField() //NumberOfComponents(3)
        // A field with value equal to FieldVector.
    {
      fFieldComponents = FieldVector;
    }

    VecTUniformMagField(double vField,
                             double vTheta,
                             double vPhi  );

    // virtual
    ~VecTUniformMagField() {}

    VecTUniformMagField(const VecTUniformMagField &p)   // : G4MagneticField(p)
    {
      fFieldComponents = p.fFieldComponents;
    }

    VecTUniformMagField& operator = (const VecTUniformMagField &p)
        // Copy constructor and assignment operator.
    {
      if (&p == this) return *this;
      // for (int i=0; i<3; i++) fFieldComponents[i] = p.fFieldComponents[i];
      fFieldComponents = p.fFieldComponents;
      return *this;
    }

    template <class Backend>
    void GetFieldValue( const vecgeom::Vector3D<typename Backend::precision_v> &, // Position,
                              vecgeom::Vector3D<typename Backend::precision_v> &FieldValue )
    {
      for (int i=0; i<3; i++) FieldValue[i] = fFieldComponents[i];
     // FieldValue= fFieldComponents;
    }

    void SetFieldValue(const vecgeom::Vector3D<double>& fieldValue)
    {
      fFieldComponents= fieldValue;
    }

    vecgeom::Vector3D<double> GetConstantFieldValue() const
    {
      return fFieldComponents;
    }
    // Return the field value

    // virtual
    VecTUniformMagField* Clone() const
    {
      return new VecTUniformMagField( *this );
    }

    VecTUniformMagField* CloneOrSafeSelf( bool /*Safe = 0*/ )
    // {  Safe= true; return this; }  //  Class is thread-safe, can use 'self' instead of clone
    // { Safe= false; return new VecTUniformMagField( this ); }  // Check ...
    { /*Safe= false;*/ return Clone(); }  // Check ...

    VecTUniformMagField* CloneOrSafeSelf( bool* pSafe )
    {
      if( pSafe ) *pSafe= true;
      return this; // ->CloneOrSafeSelf(*pSafe);
    }
    //  Class is thread-safe, can use 'self' instead of clone

  private:
    vecgeom::Vector3D<double> fFieldComponents;
};


VecTUniformMagField::VecTUniformMagField(double vField,
                                         double vTheta,
                                         double vPhi     )
{
  if ( (vField<0) || (vTheta<0) || (vTheta>Constants::pi) || (vPhi<0) || (vPhi>Constants::twopi) )
  {
     // Exception("VecTUniformMagField::VecTUniformMagField()",
     //     "GeomField0002", FatalException, "Invalid parameters.") ;
     std::cerr << "ERROR in VecTUniformMagField::VecTUniformMagField()"
               << "Invalid parameter(s): expect " << std::endl;
     std::cerr << " - Theta angle: Value = " << vTheta
               << "  Expected between 0 <= theta <= pi = " << Constants::pi << std::endl;
     std::cerr << " - Phi   angle: Value = " << vPhi
               << "  Expected between 0 <=  phi  <= 2*pi = " << Constants::twopi << std::endl;
     std::cerr << " - Magnitude vField: Value = " << vField
               << "  Expected vField > 0 " << Constants::twopi << std::endl;
  }

  //std::sin and cos should work since vTheta etc are always scalar, but what the heck
  fFieldComponents.Set( vField*sin(vTheta)*cos(vPhi),
                        vField*sin(vTheta)*sin(vPhi),
                        vField*cos(vTheta)                  );
}
#endif
