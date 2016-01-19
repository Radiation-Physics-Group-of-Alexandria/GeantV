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


class TVectorUniformMagField : public TUniformMagField
{
    typedef typename vecgeom::kVc::precision_v Double_v;
    typedef typename vecgeom::kVcFloat::precision_v Float_v;
    public:  

/*        TVectorUniformMagField(const vecgeom::Vector3D<Float_v>& FieldVector )
           : TUniformMagField() //, GUVMagneticField() 
        {
           for (int i=0; i<3; i++) fFieldComponents[i] = (Double_v) FieldVector[i];
        }*/

        TVectorUniformMagField(Double_v vField,
                               Double_v vTheta,
                               Double_v vPhi  );

        // virtual
        ~TVectorUniformMagField() {}

        TVectorUniformMagField(const TVectorUniformMagField &p)   // : G4MagneticField(p)
          : TUniformMagField(p)
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
                                  vecgeom::Vector3D<Double_v> &FieldValue )
        {
           FieldValue= fFieldComponents;
        }

        void SetFieldValue(const vecgeom::Vector3D<Double_v>& fieldValue)
        {
           fFieldComponents= fieldValue;
        }

        vecgeom::Vector3D<Double_v> GetConstantFieldValue() const
        {
           return fFieldComponents;
        }
      

    private:
        vecgeom::Vector3D<Double_v> fFieldComponents;
};

TVectorUniformMagField::TVectorUniformMagField(typename vecgeom::kVc::precision_v vField,
                                               typename vecgeom::kVc::precision_v vTheta,
                                               typename vecgeom::kVc::precision_v vPhi   )
  : TUniformMagField(vField[0], vTheta[0], vPhi[0])
{
   // if ( (vField<0) || (vTheta<0) || (vTheta>Constants::pi) || (vPhi<0) || (vPhi>Constants::twopi) )

   typedef typename vecgeom::kVc::bool_v Bool_v;
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


   fFieldComponents.Set( vField*Vc::sin(vTheta)*Vc::cos(vPhi),
                         vField*Vc::sin(vTheta)*Vc::sin(vPhi),
                         vField*Vc::cos(vTheta)                );
}
#endif
