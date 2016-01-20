//
//  First version:      (Josh) - GSoC 2014 project
//  Current version:  J. Apostolakis

#ifndef TVectorUniformMagField_H
#define TVectorUniformMagField_H

#include "GUVMagneticField.h"
#include "GUVVectorMagneticField.h"
#include "GUVVectorField.h"

#include "TUniformMagField.h"
#include <iostream>

#include "base/Vector3D.h"

#include "Constants.h"  


class TVectorUniformMagField : public GUVVectorMagneticField
{  
    typedef typename vecgeom::kVc::precision_v      Double_v;
    typedef typename vecgeom::kVcFloat::precision_v Float_v;
    public:  

        TVectorUniformMagField(const vecgeom::Vector3D<Float_v>& FieldVector )
           : GUVVectorMagneticField() 
        {
           for (int i=0; i<3; i++) fFieldComponents[i] = /*(Double_v)*/ FieldVector[i];
        }

        TVectorUniformMagField(Double_v vField,
                               Double_v vTheta,
                               Double_v vPhi  );

        // virtual
        ~TVectorUniformMagField() {}

        TVectorUniformMagField(const TVectorUniformMagField &p)   // : G4MagneticField(p)
          // : TUniformMagField(p)
        {
           fFieldComponents = p.fFieldComponents;
        }

        TVectorUniformMagField& operator = (const TVectorUniformMagField &p)
            // Copy constructor and assignment operator.
        {
           if (&p == this) return *this;
           // for (int i=0; i<3; i++) fFieldComponents[i] = p.fFieldComponents[i];
           fFieldComponents = p.fFieldComponents;
           return *this;
        }

        // virtual
        void GetFieldValue( const vecgeom::Vector3D<Double_v> &, // Position,
                                  vecgeom::Vector3D<Float_v> &FieldValue )
        {
           for (int i=0; i<3; i++) FieldValue[i] = fFieldComponents[i];
           // FieldValue= fFieldComponents;
        }

        void SetFieldValue(const vecgeom::Vector3D<Float_v>& fieldValue)
        {
           for (int i=0; i<3; i++) fFieldComponents[i] = fieldValue[i];
           // fFieldComponents= fieldValue;
        }

        vecgeom::Vector3D<Float_v> GetConstantFieldValue() const
        {
           return fFieldComponents;
        }

        //below 3 functions similar as TUniformMagField
        TVectorUniformMagField* Clone() const
        {
          //discuss
          //uncommenting the line below gives error
          //why?
          //Resolved. 
          //error was because of inheriting GUVVectorMagField from GUVField ultimately which had a virtual function with scalar signature
           return new TVectorUniformMagField( *this );
        }

        TVectorUniformMagField* CloneOrSafeSelf( bool /*Safe = 0*/ )
        { /*Safe= false;*/ return Clone(); }  

        TVectorUniformMagField* CloneOrSafeSelf( bool* pSafe )
        {
           if( pSafe ) *pSafe= true;
           return this; 
        }
      

    private:
        // vecgeom::Vector3D<float> fFieldComponents;
        vecgeom::Vector3D<Float_v> fFieldComponents;

};

// This function is not working now. Uncomment and compile to see error
// discuss
TVectorUniformMagField::TVectorUniformMagField(typename vecgeom::kVc::precision_v vField,
                                               typename vecgeom::kVc::precision_v vTheta,
                                               typename vecgeom::kVc::precision_v vPhi   )
  // : TUniformMagField(vField[0], vTheta[0], vPhi[0])
{
   // if ( (vField<0) || (vTheta<0) || (vTheta>Constants::pi) || (vPhi<0) || (vPhi>Constants::twopi) )

   bool fieldLessThanZero = !(vecgeom::IsEmpty(vField<0));
   bool thetaLessThanZero = !(vecgeom::IsEmpty(vTheta<0));
   bool thetaMoreThanPi   = !(vecgeom::IsEmpty(vTheta>Constants::pi));
   bool phiLessThanZero   = !(vecgeom::IsEmpty(vPhi<0));
   bool phiMoreThanPi     = !(vecgeom::IsEmpty(vPhi>Constants::pi));
   if ( fieldLessThanZero || thetaLessThanZero || thetaMoreThanPi || phiLessThanZero || phiMoreThanPi )
   {
      std::cerr << "ERROR in TVectorUniformMagField::TVectorUniformMagField()"
                << "Invalid parameter(s): expect " << std::endl;
      std::cerr << " - Theta angle: Value = " << vTheta
                << "  Expected between 0 <= theta <= pi = " << Constants::pi << std::endl;
      std::cerr << " - Phi   angle: Value = " << vPhi
                << "  Expected between 0 <=  phi  <= 2*pi = " << Constants::twopi << std::endl;
      std::cerr << " - Magnitude vField: Value = " << vField
                << "  Expected vField > 0 " << Constants::twopi << std::endl;
   }

  
  fFieldComponents.Set((Float_v) (vField*Vc::sin(vTheta)*Vc::cos(vPhi)),
                        (Float_v) (vField*Vc::sin(vTheta)*Vc::sin(vPhi)),
                        (Float_v) (vField*Vc::cos(vTheta))              );
}
#endif
