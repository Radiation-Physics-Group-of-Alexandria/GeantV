//
// class TemplateGUVEquationOfMotion
//
// Class description:
//
// Abstract Base Class for the right hand size of the equation of
// motion of a particle in a field.

// History:
// - Created. J.Apostolakis     Dec 2014/Jan 2015
// -------------------------------------------------------------------

#ifndef TemplateGUVEquationOfMotion_H
#define TemplateGUVEquationOfMotion_H

#include <cassert>
#include <iostream>

// #include <vector>
#include "base/Vector3D.h"
#include <ostream>
#include "TemplateGUVField.h"   // required in inline method implementations

// #define DEBUGAnanya

//vecgeom::VECGEOM_DEVICE_FORWARD_DECLARE( template <typename Type> class TemplateGUVEquationOfMotion;)

template <class Backend>
class TemplateGUVEquationOfMotion //: public GUVEquationOfMotion 
{
  
  typedef typename Backend::precision_v Double_v;

  public:  // with description

     TemplateGUVEquationOfMotion( TemplateGUVField<Backend> *Field, unsigned int verbose=0 );

     virtual ~TemplateGUVEquationOfMotion();
       // Constructor and virtual destructor. No operations, just checks

     virtual void EvaluateRhsGivenB( const Double_v yVec[],
                                     const vecgeom::Vector3D<Double_v> B,  // Was double B[3],
                                           Double_v charge,
                                           Double_v dydx[]   ) const = 0;
       // Given the value of the  field "B", this function 
       // calculates the value of the derivative dydx.
       // --------------------------------------------------------
       // This is the _only_ function a subclass must define.
       // The other two functions use Rhs_givenB.

      
      virtual void InitializeCharge(double particleCharge)=0;
       // Must be called to correctly initialise and provide charge
      virtual void InvalidateParameters()=0;
      
      inline void InformReady(); // All parameters have been set (charge+)
      inline void InformDone();  // Invalidate charge, other parameters
      inline void CheckInitialization() const; // Ensure initialization
      inline void CheckDone() const;

     // virtual void SetChargeMomentumMass(double particleCharge,
     //                                    double MomentumXc,
     //                                    double MassXc2) = 0;
     //   // Set the charge, momentum and mass of the current particle
     //   // --> used to set the equation's coefficients ...

      inline void RightHandSide( const  Double_v y[],
                                        Double_v charge ,
                                        Double_v dydx[] ) const;
       // This calculates the value of the derivative dydx at y.
       // It is the usual enquiry function.
       // ---------------------------
       // It uses the virtual function EvaluateRhsGivenB

     void EvaluateRhsReturnB( const Double_v y[],
                                    Double_v dydx[],
                                    Double_v charge,
                                    vecgeom::Vector3D<Double_v> &Field ) const;
       // Same as RHS above, but also returns the value of B.
       // Should be made the new default ? after putting dydx & B in a class.

     void GetFieldValue( const Double_v Point[4],
                               Double_v Field[] )  const;
       // Obtain only the field - the stepper assumes it is pure Magnetic.
       // Not protected, because GUVRKG3_Stepper uses it directly.
     inline
     void GetFieldValue( const Double_v                     Point[4],
                               vecgeom::Vector3D<Double_v>  &FieldValue ) const;

     inline
     void GetFieldValue( const vecgeom::Vector3D<Double_v> &Position,
                               vecgeom::Vector3D<Double_v> &FieldValue ) const;

     const TemplateGUVField<Backend>* GetFieldObj() const {return fField;}
           TemplateGUVField<Backend>* GetFieldObj()       {return fField;}

     void SetFieldObj(TemplateGUVField<Backend>* pField){fField=pField;}

     bool         Initialised() const { return fInitialised; } 
     unsigned int GetId() const       { return fEquationId; }
     static unsigned int GetNumCreated() { return fNumObjectsCreated; }
     static unsigned int GetNumLive() { return fNumObjectsCreated - fNumObjectsDeleted; }

     template <class Backend_>
     friend std::ostream&
             operator<<( std::ostream& os, const TemplateGUVEquationOfMotion<Backend_>& eq);

  public:
     static const unsigned int idxTime=3;  // Convention for location of time 't' in vector

  private:
     static unsigned int fNumObjectsCreated;
     static unsigned int fNumObjectsDeleted;
     // const int GUVmaximum_number_of_field_components = 24;
     enum { GUVmaximum_number_of_field_components = 24 } ;

     TemplateGUVField<Backend> *     fField;
     unsigned int   fEquationId;  //
     unsigned short fVerbose;
     bool           fInitialised;


};

template <class Backend>
unsigned int TemplateGUVEquationOfMotion<Backend>::fNumObjectsCreated=0;

template <class Backend>
unsigned int TemplateGUVEquationOfMotion<Backend>::fNumObjectsDeleted=0;

template <class Backend>
inline
TemplateGUVEquationOfMotion<Backend>::TemplateGUVEquationOfMotion(TemplateGUVField<Backend>* pField, unsigned int verbose)
   : // GUVEquationOfMotion(pField, verbose),
     fField(pField), 
     fEquationId(fNumObjectsCreated++),
     fVerbose(verbose), 
     fInitialised(false)
{
   if(fVerbose)
   {
     std::cout << " Called TemplateGUVEquationOfMotion::InformDone() " << std::endl;
     // std::cout << *this << std::endl;
   }
  #ifdef DEBUGAnanya
     std::cout<<"----Entered constructor of TemplateGUVEquationOfMotion "<<std::endl;
  #endif

}

template <class Backend>
inline
void TemplateGUVEquationOfMotion<Backend>::InformReady() // was Initialize()
{
   fInitialised= true;
}

template <class Backend>
inline
void TemplateGUVEquationOfMotion<Backend>::InformDone()  // was Clear() and before Finished();
{
   if(fVerbose)
   {
     std::cout << " Called TemplateGUVEquationOfMotion::InformDone() " << std::endl;
     // std::cout << *this << std::endl;
   }
   assert( fInitialised );
   fInitialised= false;
}

template <class Backend>
inline
void TemplateGUVEquationOfMotion<Backend>::GetFieldValue( const typename Backend::precision_v Point[4],
                                                                typename Backend::precision_v Field[] ) const
{
   vecgeom::Vector3D<typename Backend::precision_v> Position( Point[0], Point[1], Point[2] );
   vecgeom::Vector3D<typename Backend::precision_v> FieldVec;
   fField-> GetFieldValue( Position, FieldVec );
   Field[0] = (Double_v) FieldVec[0];
   Field[1] = (Double_v) FieldVec[1];
   Field[2] = (Double_v) FieldVec[2];
}

template <class Backend>
inline
void TemplateGUVEquationOfMotion<Backend>::GetFieldValue( const Double_v Point[4],
                                                                vecgeom::Vector3D<typename Backend::precision_v>  &FieldValue ) const
{
   vecgeom::Vector3D<Double_v> Position( Point[0], Point[1], Point[2] );
   fField-> GetFieldValue( Position, FieldValue );
}

template <class Backend>
inline
void TemplateGUVEquationOfMotion<Backend>::GetFieldValue( const vecgeom::Vector3D<typename Backend::precision_v> &Position,
                                                                vecgeom::Vector3D<typename Backend::precision_v> &FieldValue ) const
{
   fField-> GetFieldValue( Position, FieldValue );
}

template <class Backend>
inline
void
TemplateGUVEquationOfMotion<Backend>::RightHandSide( const typename Backend::precision_v y[],
                                                           typename Backend::precision_v charge,
                                                           typename Backend::precision_v dydx[] ) const
{
   using ThreeVectorD = vecgeom::Vector3D<typename Backend::precision_v>;
   CheckInitialization();

   ThreeVectorD  Field_3vf;
   ThreeVectorD  Position( y[0], y[1], y[2] );

   GetFieldValue    ( Position, Field_3vf );
   EvaluateRhsGivenB( y, Field_3vf, charge, dydx );
   // std::cout<<"\n----Field_3vf is: "<<Field_3vf[0]<<std::endl;
   // std::cout<<"----Field_3vf is: "<<Field_3vf[1]<<std::endl;
   // std::cout<<"----Field_3vf is: "<<Field_3vf[2]<<std::endl;
/*
   #ifdef DEBUGAnanya
    std::cout<<"\n----Field_3vf is: "<<Field_3vf[0]<<std::endl;
    std::cout<<"----dydx is: "<<dydx[0]<<std::endl;
   #endif */
}

#include <iostream>

template <class Backend>
void TemplateGUVEquationOfMotion<Backend>::CheckInitialization() const
{
#ifdef GUVERBOSE
   if( fVerbose && !fInitialised ){
      std::cerr << "TemplateGUVEquationOfMotion is not Initialised" << std::endl;
   }
#endif
   // assert( fInitialised ); //Ananya:temp
}

template <class Backend>
void TemplateGUVEquationOfMotion<Backend>::CheckDone() const
{
#ifdef GUVERBOSE
   if( fVerbose && fInitialised ){
      std::cerr << "TemplateGUVEquationOfMotion was NOT told it is Done!" << std::endl;
   }
#endif
   assert( !fInitialised );
}

template <class Backend>
TemplateGUVEquationOfMotion<Backend>::~TemplateGUVEquationOfMotion()
{
  CheckDone();
  fNumObjectsDeleted++;
}

template <class Backend>
void
TemplateGUVEquationOfMotion<Backend>::EvaluateRhsReturnB( const typename Backend::precision_v  y[],
                                                                typename Backend::precision_v  dydx[],
                                                                typename Backend::precision_v  charge,
                                                                vecgeom::Vector3D<typename Backend::precision_v> &Field ) const
{
   typedef typename Backend::precision_v Double_v;
   Double_v  PositionAndTime[4];
   PositionAndTime[0] = y[0];
   PositionAndTime[1] = y[1];
   PositionAndTime[2] = y[2];
   // PositionAndTime[3] = y[7];

   GetFieldValue( PositionAndTime, Field) ;
   EvaluateRhsGivenB( y, Field, charge, dydx );
}

template <class Backend>
std::ostream&  operator<<( std::ostream& os, const TemplateGUVEquationOfMotion<Backend>& eq)
{
   os << " Equation of Motion # " << eq.GetId()
      << "   field ptr= "  << eq.GetFieldObj() << "  Initialised= " << eq.Initialised()
      << std::endl;
   // os << "  Total # of E-of-M = " << TemplateGUVEquationOfMotion<Backend>::GetNumCreated()
   //    << " live= " << TemplateGUVEquationOfMotion<Backend>::GetNumLive() << std::endl;

   return os;
}

// template <typename T>
// std::ostream& operator<<(std::ostream& os, Vector3D<T> const &vec) {
//   os << "(" << vec[0] << ", " << vec[1] << ", " << vec[2] << ")";
//   return os;
// }


#endif /* TemplateGUVEquationOfMotion_DEF */
