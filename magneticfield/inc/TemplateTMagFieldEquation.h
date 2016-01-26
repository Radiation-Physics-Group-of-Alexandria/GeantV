// Approach is derived from the Geant4 class G4MagFieldEquation
//

#include <cmath>

#include "TemplateGUVEquationOfMotion.h"

#include "base/Global.h"

//  Ensure that equation Right Hand Side is inlined - may be compiler dependend
#define INLINERHS 1

#ifdef  INLINERHS
#define REALLY_INLINE   inline __attribute__((always_inline)) 
#else
#define REALLY_INLINE   inline
#endif

#ifndef TEMPLATETMAGFIELDEQUATION_H
#define TEMPLATETMAGFIELDEQUATION_H  1

// #include <vector>
#include "base/Vector3D.h"

#include "Units.h"
#include "Constants.h"
//  Update to GeantV units ASAP

template
<class Backend, class Field, unsigned int Size>
class TemplateTMagFieldEquation :  public TemplateGUVEquationOfMotion<Backend>
{
   public:
     typedef typename Backend::precision_v      Double_v;
     typedef Field T_Field;

     static const unsigned int  N   = Size;
     static constexpr double fCof   = Constants::c_light;  

     TemplateTMagFieldEquation(T_Field* pF) : TemplateGUVEquationOfMotion<Backend>(pF) { fPtrField = pF; }

     TemplateTMagFieldEquation(const TemplateTMagFieldEquation& );
     ~TemplateTMagFieldEquation()  {}  // Was virtual - but now no inheritance

     TemplateTMagFieldEquation<Backend,Field,Size>* Clone() const;
     TemplateTMagFieldEquation<Backend,Field,Size>* CloneOrSafeSelf(bool& safe);
     TemplateTMagFieldEquation<Backend,Field,Size>* CloneOrSafeSelf(bool* safe=0);
     
     REALLY_INLINE  
     void GetFieldValue(const Double_v Point[4],
                              Double_v Value[]) const
     {
        fPtrField->T_Field::GetFieldValue(Point, Value);
     }

     inline 
     void RightHandSide(const Double_v y[], 
                        const Double_v charge, 
                              Double_v dydx[] ) const;

     REALLY_INLINE
     void TEvaluateRhsGivenB( const Double_v y[],
                              const vecgeom::Vector3D<Double_v> B,  // Was double B[3],
                              const Double_v charge,
                                    Double_v dydx[] ) const;

     // virtual
     void EvaluateRhsGivenB( const Double_v y[],
                             const vecgeom::Vector3D<Double_v> B,  // Was const double B[3],
                             const Double_v charge,
                                   Double_v dydx[] ) const
     { TEvaluateRhsGivenB( y, B, charge, dydx); }

     REALLY_INLINE
     void FieldFromY(const Double_v y[], 
                     const Double_v charge,
                           Double_v Bfield[] ) const;

     REALLY_INLINE
     void FieldFromY(const Double_v y[], 
                    /* const Double_v charge,*/
                           vecgeom::Vector3D<Double_v> &Bfield ) const;

     REALLY_INLINE
     void PrintInputFieldAndDyDx(const Double_v y[],  
                                 const Double_v charge,  
                                       Double_v dydx[] ) const;

     //Discuss
     //Function needed? probably not if we don't care about particleCharge
     //or not take any input but just do InformReady?
     REALLY_INLINE
     void InitializeCharge(double particleCharge) final 
      {  fParticleCharge= particleCharge; TemplateGUVEquationOfMotion<Backend>::InformReady();  }

      //should get this func. from inheritance
      void InvalidateParameters() final { TemplateGUVEquationOfMotion<Backend>::InformDone();}

   private:
     enum { G4maximum_number_of_field_components = 24 };
     T_Field *fPtrField;
     double   fParticleCharge;
};

template
<class Backend, class Field, unsigned int Size>
   TemplateTMagFieldEquation<Backend,Field,Size>
   ::TemplateTMagFieldEquation(const TemplateTMagFieldEquation& right)
   :  TemplateGUVEquationOfMotion<Backend>( (TemplateGUVField<Backend>*) 0 ),
      fPtrField( right.fPtrField->CloneOrSafeSelf( (bool *)0 ) )
      // fPtrField( new Field(right.fPtrField) )
{
   // G4bool threadSafe;
   // fPtrField = right.fPtrField->CloneOrSafeSelf( &threadSafe );

   // std::cout <<  "TemplateTMagFieldEquation - copy constructor called." << std::endl;
   TemplateGUVEquationOfMotion<Backend>::SetFieldObj( fPtrField ); //  Also stored in base class ... for now
}

template
<class Backend, class Field, unsigned int Size>
   TemplateTMagFieldEquation<Backend,Field,Size>*
   TemplateTMagFieldEquation<Backend,Field,Size>
   ::CloneOrSafeSelf(bool& safe)
{
   TemplateTMagFieldEquation<Backend,Field,Size>* equation;
   Field* pField=
      fPtrField->CloneOrSafeSelf(safe);

   std::cerr << " #TemplateTMagFieldEquation<Backend,Field,Size>::CloneOrSafeSelf(bool& safe) called# " << std::endl;

   equation = new TemplateTMagFieldEquation( pField );
   safe= false;

}



template
<class Backend, class Field, unsigned int Size>
REALLY_INLINE
   void  TemplateTMagFieldEquation<Backend,Field, Size>
   ::TEvaluateRhsGivenB( const typename Backend::precision_v y[],
                         const vecgeom::Vector3D<typename Backend::precision_v> B, //Bfloat, 
                         const typename Backend::precision_v charge,
                               typename Backend::precision_v dydx[]  ) const
{
  
    typedef typename Backend::precision_v Double_v;
    
    Double_v momentum_mag_square = y[3]*y[3] + y[4]*y[4] + y[5]*y[5];
    Double_v inv_momentum_magnitude = 1. / vecgeom::VECGEOM_IMPL_NAMESPACE::Sqrt( momentum_mag_square );

    // vecgeom::Vector3D<Double_v> B( (Double_v) Bfloat[0], (Double_v) Bfloat[1], (Double_v) Bfloat[2] );

    Double_v cof = fParticleCharge*fCof*inv_momentum_magnitude;

    dydx[0] = y[3]*inv_momentum_magnitude;       //  (d/ds)x = Vx/V
    dydx[1] = y[4]*inv_momentum_magnitude;       //  (d/ds)y = Vy/V
    dydx[2] = y[5]*inv_momentum_magnitude;       //  (d/ds)z = Vz/V

    dydx[3] = cof*(y[4]*B[2] - y[5]*B[1]) ;  // Ax = a*(Vy*Bz - Vz*By)
    dydx[4] = cof*(y[5]*B[0] - y[3]*B[2]) ;  // Ay = a*(Vz*Bx - Vx*Bz)
    dydx[5] = cof*(y[3]*B[1] - y[4]*B[0]) ;  // Az = a*(Vx*By - Vy*Bx)

    return ;
}


template
<class Backend, class Field, unsigned int Size>
REALLY_INLINE
void
TemplateTMagFieldEquation<Backend,Field,Size>
   ::FieldFromY(const typename Backend::precision_v y[], 
                const typename Backend::precision_v charge,  
                      typename Backend::precision_v Bfield[3] ) const
{
    // double  Bfield[3];  //G4maximum_number_of_field_components];
    typename Backend::precision_v PositionAndTime[4];
    PositionAndTime[0] = y[0];
    PositionAndTime[1] = y[1];
    PositionAndTime[2] = y[2];
    PositionAndTime[3] = 0;        // Time

    GetFieldValue(PositionAndTime, Bfield) ;
}

template
<class Backend, class Field, unsigned int Size>
REALLY_INLINE
void
TemplateTMagFieldEquation<Backend,Field,Size>
   ::FieldFromY(const typename Backend::precision_v                                  y[],  
                               vecgeom::Vector3D<typename Backend::precision_v> &Bfield ) const
{
    vecgeom::Vector3D<typename Backend::precision_v> Position( y[0], y[1], y[2] );

    fPtrField->T_Field::GetFieldValue( Position, Bfield );
}


template
<class Backend, class Field, unsigned int Size>
REALLY_INLINE
void
TemplateTMagFieldEquation<Backend,Field,Size>
   ::RightHandSide(const typename Backend::precision_v y[], 
                   const typename Backend::precision_v charge, 
                         typename Backend::precision_v dydx[] ) const
{
    vecgeom::Vector3D<Double_v> BfieldVec;

    FieldFromY( y, BfieldVec );
    TEvaluateRhsGivenB( y, BfieldVec, charge, dydx );
}

#include <iostream>   // For debuging only
using std::cout;
using std::endl;

template
<class Backend, class Field, unsigned int Size>
REALLY_INLINE
void
TemplateTMagFieldEquation<Backend,Field,Size>
   ::PrintInputFieldAndDyDx(const typename Backend::precision_v y[], 
                            const typename Backend::precision_v charge, 
                                  typename Backend::precision_v dydx[] ) const
{

    RightHandSide(y, dydx);

    // Obtain the field value
    typedef typename Backend::precision_v Double_v;
    Double_v  Bfield[3];  //G4maximum_number_of_field_components];
    FieldFromY( y, charge, Bfield );
    TEvaluateRhsGivenB(y, Bfield, charge, dydx);

    cout.precision(8);

    // cout.setf (std::ios_base::fixed);
    // cout << " Position = " << y[0] << " " << y[1] << " " << y[3] << endl;
    // cout.unsetf(std::ios_base::fixed);
    cout << "\n# Input & B field \n";
    cout.setf (std::ios_base::scientific);
    cout << " Position = " << y[0] << " " << y[1] << " " << y[2] << endl;
    cout << " Momentum = " << y[3] << " " << y[4] << " " << y[5] << endl;
    cout << " B-field  = " << Bfield[0] << " " << Bfield[1] << " " << Bfield[2] << endl;
    cout.unsetf(std::ios_base::scientific);

    cout << "\n# 'Force' from B field \n";
    cout.setf (std::ios_base::fixed);
    cout << " dy/dx [0-2] (=dX/ds) = " << dydx[0]   << " " << dydx[1]   << " " << dydx[2] << endl;
    cout << " dy/dx [3-5] (=dP/ds) = " << dydx[3]   << " " << dydx[4]   << " " << dydx[5] << endl;
    cout.unsetf(std::ios_base::fixed);
}
#endif  // TMAGFIELDEQUATION_H
