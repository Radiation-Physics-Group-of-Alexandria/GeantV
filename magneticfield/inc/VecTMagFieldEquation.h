// Approach is derived from the Geant4 class G4MagFieldEquation
//

#include <cmath>

#include "base/Global.h"

//  Ensure that equation Right Hand Side is inlined - may be compiler dependend
#define INLINERHS 1

#ifdef  INLINERHS
#define REALLY_INLINE   inline __attribute__((always_inline)) 
#else
#define REALLY_INLINE   inline
#endif

#ifndef VECTMAGFIELDEQUATION_H
#define VECTMAGFIELDEQUATION_H  1

// #include <vector>
#include "base/Vector3D.h"

#include "Units.h"
#include "Constants.h"
//  Update to GeantV units ASAP

#include <cassert>
#include <iostream>

#include <ostream>
#include "VecGUVField.h"   // required in inline method implementations

#define DEBUGAnanya

template
<class T_Field, unsigned int Size>
class VecTMagFieldEquation // :  public VecGUVEquationOfMotion
{
  public:
    static const unsigned int  N   = Size;
    static constexpr double fCof   = Constants::c_light;  

    VecTMagFieldEquation(T_Field* pF, unsigned int verbose=0);

    VecTMagFieldEquation(const VecTMagFieldEquation& );
     ~VecTMagFieldEquation()  { CheckDone(); fNumObjectsDeleted++; }  // Was virtual - but now no inheritance

    VecTMagFieldEquation<T_Field,Size>* Clone() const;
    VecTMagFieldEquation<T_Field,Size>* CloneOrSafeSelf(bool& safe);
    VecTMagFieldEquation<T_Field,Size>* CloneOrSafeSelf(bool* safe=0);
     
    template <class Backend>
    REALLY_INLINE  
    void GetFieldValue(const typename Backend::precision_v Point[4],
                             typename Backend::precision_v Value[]) const;

    template <class Backend>
    inline
    void GetFieldValue( const typename Backend::precision_v                     Point[4],
                              vecgeom::Vector3D<typename Backend::precision_v>  &FieldValue ) const;


    template <class Backend>
    inline
    void GetFieldValue( const vecgeom::Vector3D<typename Backend::precision_v> &Position,
                              vecgeom::Vector3D<typename Backend::precision_v> &FieldValue ) const;

    template <class Backend>
    inline 
    void RightHandSide(const typename Backend::precision_v y[], 
                       const typename Backend::precision_v charge, 
                             typename Backend::precision_v dydx[]) const;
    
    template <class Backend>
    void RightHandSide(const typename Backend::precision_v y[],  
                             typename Backend::precision_v dydx[]) const
    { typename Backend::precision_v charge  = -1.0;
      RightHandSide(y, charge, dydx);  }; //Ananya
      //added this function to get RightHandSide functions compatible irrespecitive of 
      //whether charge is given in input or not. 
      //Assumed that in final version, charge will be included everywhere.

    template <class Backend>
    REALLY_INLINE
    void TEvaluateRhsGivenB( const typename Backend::precision_v y[],
                             const vecgeom::Vector3D<typename Backend::precision_v> B,  // Was double B[3],
                             const typename Backend::precision_v charge= -1.,
                                   typename Backend::precision_v dydx[]= 0. ) const;

    // virtual
    template <class Backend>
    void EvaluateRhsGivenB( const typename Backend::precision_v y[],
                            const vecgeom::Vector3D<typename Backend::precision_v> B,  // Was const double B[3],
                            const typename Backend::precision_v charge= -1,
                                  typename Backend::precision_v dydx[]= 0. ) const
    { TEvaluateRhsGivenB( y, B, charge, dydx); }

    template <class Backend>
    void EvaluateRhsReturnB( const typename Backend::precision_v y[],
                                   typename Backend::precision_v dydx[],
                                   typename Backend::precision_v charge,
                                   vecgeom::Vector3D<typename Backend::precision_v> &Field ) const;
     
    template <class Backend>
    REALLY_INLINE
    void FieldFromY(const typename Backend::precision_v y[], 
                    const typename Backend::precision_v charge,
                          typename Backend::precision_v Bfield[] ) const;

    template <class Backend>
    REALLY_INLINE
    void FieldFromY(const typename Backend::precision_v y[], 
                   /* const typename Backend::precision_v charge,*/
                          vecgeom::Vector3D<typename Backend::precision_v> &Bfield ) const;

    template <class Backend>
    REALLY_INLINE
    void PrintInputFieldAndDyDx(const typename Backend::precision_v y[],  
                                const typename Backend::precision_v charge,  
                                      typename Backend::precision_v dydx[] ) const;

    void InvalidateParameters() { InformDone();}


    inline void InformReady(); // All parameters have been set (charge+)
    inline void InformDone();  // Invalidate charge, other parameters
    inline void CheckInitialization() const; // Ensure initialization
    inline void CheckDone() const;

    const T_Field* GetFieldObj() const {return fPtrField;}
          T_Field* GetFieldObj()       {return fPtrField;}

    bool         Initialised() const { return fInitialised; } 
    unsigned int GetId() const       { return fEquationId; }
    static unsigned int GetNumCreated() { return fNumObjectsCreated; }
    static unsigned int GetNumLive() { return fNumObjectsCreated - fNumObjectsDeleted; }

    void SetFieldObj(T_Field* pField){fPtrField = pField;}

    template <class T_Field_, unsigned int Size_>
    friend std::ostream&
           operator<<( std::ostream& os, const VecTMagFieldEquation<T_Field_, Size_>& eq);


  public:
    static const unsigned int idxTime=3;  // Convention for location of time 't' in vector

  private:
    static unsigned int fNumObjectsCreated;
    static unsigned int fNumObjectsDeleted;

    enum { G4maximum_number_of_field_components = 24 };

    T_Field *fPtrField;

    unsigned int   fEquationId;  //
    unsigned short fVerbose;
    bool           fInitialised;
};

template <class T_Field, unsigned int Size>
unsigned int VecTMagFieldEquation<T_Field,Size>::fNumObjectsCreated=0;

template <class T_Field, unsigned int Size>
unsigned int VecTMagFieldEquation<T_Field,Size>::fNumObjectsDeleted=0;

template <class T_Field, unsigned int Size>
inline
VecTMagFieldEquation<T_Field,Size>
  ::VecTMagFieldEquation(T_Field* pField, unsigned int verbose)
  : // GUVEquationOfMotion(pField, verbose),
    fPtrField(pField), 
    fEquationId(fNumObjectsCreated++),
    fVerbose(verbose), 
    fInitialised(false)
{
  if(fVerbose)
  {
    std::cout << " Called VecTMagFieldEquation::InformDone() " << std::endl;
    // std::cout << *this << std::endl;
  }
 #ifdef DEBUGAnanya
  std::cout<<"----Entered constructor of VecTMagFieldEquation "<<std::endl;
 #endif

}

template
<class T_Field, unsigned int Size>
  VecTMagFieldEquation<T_Field,Size>
  ::VecTMagFieldEquation(const VecTMagFieldEquation& right)
  :  // VecGUVEquationOfMotion( (VecGUVField*) 0 ), // Commented for now, discuss what to do : Ananya
    fPtrField( right.fPtrField->CloneOrSafeSelf( (bool *)0 ) )
    // fPtrField( new T_Field(right.fPtrField) )
{
  // G4bool threadSafe;
  // fPtrField = right.fPtrField->CloneOrSafeSelf( &threadSafe );
  // std::cout <<  "VecTMagFieldEquation - copy constructor called." << std::endl;
  SetFieldObj( fPtrField ); //  Also stored in base class ... for now
}

template
<class T_Field, unsigned int Size>
  VecTMagFieldEquation<T_Field,Size>*
  VecTMagFieldEquation<T_Field,Size>
  ::CloneOrSafeSelf(bool& safe)
{
  VecTMagFieldEquation<T_Field,Size>* equation;
  T_Field* pField=
    fPtrField->CloneOrSafeSelf(safe);

  std::cerr << " #VecTMagFieldEquation<T_Field,Size>::CloneOrSafeSelf(bool& safe) called# " << std::endl;

  equation = new VecTMagFieldEquation( pField );
  safe= false;
}

template <class T_Field, unsigned int Size>
inline
void VecTMagFieldEquation<T_Field,Size>
  ::InformReady() // was Initialize()
{
  fInitialised= true;
}


template <class T_Field, unsigned int Size>
inline
void VecTMagFieldEquation<T_Field,Size>
  ::InformDone()  // was Clear() and before Finished();
{
  if(fVerbose)
  {
    std::cout << " Called VecTMagFieldEquation<T_Field,Size>::InformDone() " << std::endl;
    // std::cout << *this << std::endl;
  }
  assert( fInitialised );
  fInitialised= false;
}

template <class T_Field, unsigned int Size>
void VecTMagFieldEquation<T_Field,Size>
  ::CheckInitialization() const
{
#ifdef GUVERBOSE
  if( fVerbose && !fInitialised ){
    std::cerr << "VecTMagFieldEquation<T_Field,Size> is not Initialised" << std::endl;
  }
#endif
  // assert( fInitialised ); //Ananya:temp
}

template <class T_Field, unsigned int Size>
void VecTMagFieldEquation<T_Field,Size>
  ::CheckDone() const
{
#ifdef GUVERBOSE
  if( fVerbose && fInitialised ){
    std::cerr << "VecTMagFieldEquation<T_Field,Size> was NOT told it is Done!" << std::endl;
  }
#endif
  assert( !fInitialised );
}

template <class T_Field, unsigned int Size>
template <class Backend>
REALLY_INLINE
void VecTMagFieldEquation<T_Field,Size>
  ::GetFieldValue( const typename Backend::precision_v Point[4],
                         typename Backend::precision_v Field[] ) const
{
   vecgeom::Vector3D<typename Backend::precision_v> Position( Point[0], Point[1], Point[2] );
   vecgeom::Vector3D<typename Backend::precision_v> FieldVec;
   fPtrField-> T_Field::GetFieldValue( Position, FieldVec );
   Field[0] = FieldVec[0];
   Field[1] = FieldVec[1];
   Field[2] = FieldVec[2];
}

template <class T_Field, unsigned int Size>
template <class Backend>
inline
void VecTMagFieldEquation<T_Field,Size>
  ::GetFieldValue( const typename Backend::precision_v Point[4],
                         vecgeom::Vector3D<typename Backend::precision_v>  &FieldValue ) const
{
   vecgeom::Vector3D<typename Backend::precision_v> Position( Point[0], Point[1], Point[2] );
   fPtrField-> GetFieldValue<Backend>( Position, FieldValue );
}

template <class T_Field, unsigned int Size>
template <class Backend>
inline
void VecTMagFieldEquation<T_Field,Size>
  ::GetFieldValue( const vecgeom::Vector3D<typename Backend::precision_v> &Position,
                         vecgeom::Vector3D<typename Backend::precision_v> &FieldValue ) const
{
   fPtrField-> GetFieldValue<Backend>( Position, FieldValue );
}

template <class T_Field, unsigned int Size>
template <class Backend>
REALLY_INLINE
void  VecTMagFieldEquation<T_Field,Size>
   ::TEvaluateRhsGivenB( const typename Backend::precision_v y[],
                         const vecgeom::Vector3D<typename Backend::precision_v> B, //Bfloat, 
                         const typename Backend::precision_v charge,
                               typename Backend::precision_v dydx[]  ) const
{
  
  typedef typename Backend::precision_v Double_v;
  
  Double_v momentum_mag_square = y[3]*y[3] + y[4]*y[4] + y[5]*y[5];
  Double_v inv_momentum_magnitude = 1. / vecgeom::VECGEOM_IMPL_NAMESPACE::Sqrt( momentum_mag_square );
/*
  #ifdef DEBUGAnanya
    std::cout<<"\n----y is: "<<y[3]<<" "<<y[4]<<" " <<y[5]<<std::endl;
    std::cout<<"----inv_momentum is: "<<inv_momentum_magnitude<<std::endl;
    std::cout<<"----momentum is: "<< momentum_mag_square <<std::endl;
  #endif*/


  // std::cout<<"\n\n\n AM I BEING CALLED SOMEHOW?"<<std::endl;
  // vecgeom::Vector3D<Double_v> B( (Double_v) Bfloat[0], (Double_v) Bfloat[1], (Double_v) Bfloat[2] );

  Double_v cof = charge*fCof*inv_momentum_magnitude;

  dydx[0] = y[3]*inv_momentum_magnitude;       //  (d/ds)x = Vx/V
  dydx[1] = y[4]*inv_momentum_magnitude;       //  (d/ds)y = Vy/V
  dydx[2] = y[5]*inv_momentum_magnitude;       //  (d/ds)z = Vz/V

  dydx[3] = cof*(y[4]*B[2] - y[5]*B[1]) ;  // Ax = a*(Vy*Bz - Vz*By)
  dydx[4] = cof*(y[5]*B[0] - y[3]*B[2]) ;  // Ay = a*(Vz*Bx - Vx*Bz)
  dydx[5] = cof*(y[3]*B[1] - y[4]*B[0]) ;  // Az = a*(Vx*By - Vy*Bx)

  return ;
}


template <class T_Field, unsigned int Size>
template <class Backend>
REALLY_INLINE
void
VecTMagFieldEquation<T_Field,Size>
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

  GetFieldValue<Backend>(PositionAndTime, Bfield) ;
}

template <class T_Field, unsigned int Size>
template <class Backend>
REALLY_INLINE
void
VecTMagFieldEquation<T_Field,Size>
 ::FieldFromY(const typename Backend::precision_v                                  y[],  
                             vecgeom::Vector3D<typename Backend::precision_v> &Bfield ) const
{
  vecgeom::Vector3D<typename Backend::precision_v> Position( y[0], y[1], y[2] );

  fPtrField->T_Field::GetFieldValue( Position, Bfield );
}


template <class T_Field, unsigned int Size>
template <class Backend>
REALLY_INLINE
void
VecTMagFieldEquation<T_Field,Size>
 ::RightHandSide(const typename Backend::precision_v y[], 
                 const typename Backend::precision_v charge, 
                       typename Backend::precision_v dydx[] ) const
{
  vecgeom::Vector3D<typename Backend::precision_v> BfieldVec;

  FieldFromY<Backend>( y, BfieldVec );
  TEvaluateRhsGivenB<Backend>( y, BfieldVec, charge, dydx );
}

#include <iostream>   // For debuging only
using std::cout;
using std::endl;

template <class T_Field, unsigned int Size>
template <class Backend>
REALLY_INLINE
void
VecTMagFieldEquation<T_Field,Size>
 ::PrintInputFieldAndDyDx(const typename Backend::precision_v y[], 
                          const typename Backend::precision_v charge, 
                                typename Backend::precision_v dydx[] ) const
{
  RightHandSide<Backend>(y, dydx);

  // Obtain the field value
  typedef typename Backend::precision_v Double_v;
  Double_v  Bfield[3];  //G4maximum_number_of_field_components];
  FieldFromY<Backend>( y, charge, Bfield );
  TEvaluateRhsGivenB<Backend>(y, Bfield, charge, dydx);

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

template <class T_Field, unsigned int Size>
template <class Backend>
void
VecTMagFieldEquation<T_Field,Size>
  ::EvaluateRhsReturnB( const typename Backend::precision_v  y[],
                              typename Backend::precision_v  dydx[],
                              typename Backend::precision_v  charge,
                              vecgeom::Vector3D<typename Backend::precision_v> &Field ) const
{
  typedef typename Backend::precision_v Double_v;
  Double_v  PositionAndTime[4];
  PositionAndTime[0] = y[0];
  PositionAndTime[1] = y[1];
  PositionAndTime[2] = y[2];
  PositionAndTime[3] = y[7];  

  GetFieldValue<Backend>( PositionAndTime, Field) ;
  EvaluateRhsGivenB<Backend>( y, Field, charge, dydx );
}

template <class T_Field, unsigned int Size>
std::ostream&  operator<<( std::ostream& os, const VecTMagFieldEquation<T_Field,Size>& eq)
{
  os << " Equation of Motion # " << eq.GetId()
     << "   field ptr= "  << eq.GetFieldObj() << "  Initialised= " << eq.Initialised()
     << std::endl;
  // os << "  Total # of E-of-M = " << VecGUVEquationOfMotion<Backend>::GetNumCreated()
  //    << " live= " << VecGUVEquationOfMotion<Backend>::GetNumLive() << std::endl;

  return os;
}

#endif  // VECTMAGFIELDEQUATION_H
