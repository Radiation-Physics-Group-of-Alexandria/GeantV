#ifndef _FIELDTRACK_H_
#define _FIELDTRACK_H_

/*---------------------
Data structure in place of GUFieldTrack to be used 
for input and output stream arrays of AccurateAdvance 
in IntegrationDriver. Functions DumpToArray and LoadFromArray
can be removed if PosMomVector is made public data member.
Same goes for SetCurveLength and GetCurveLength functions.
----------------*/
#include <iostream>
#include "base/Vector3D.h"   // VecGeom/base/Vector3D.h 

struct FieldTrack{
  
  typedef vecgeom::Vector3D<double> ThreeVector;

private: 
  double fDistanceAlongCurve = 0.0;

public:

  FieldTrack()= default ;
  FieldTrack( const ThreeVector& pPosition,
              const ThreeVector& pMomentum); // {};

  FieldTrack( const FieldTrack& rStVec );
  inline FieldTrack& operator = ( const FieldTrack & rStVec );

  //data members
  double PosMomVector[6];

  //And functions 
  void DumpToArray(double valArr[]){ //12 from ncompSVEC as in both TemplateGUIntegrationDriver
    for (int i = 0; i < 6; ++i)        //and GUFieldTrack function
    {
      valArr[i] = PosMomVector[i];
    }
  }

  void LoadFromArray(const double valArr[], int noVarsIntegrated = 6){
    for (int i = 0; i < noVarsIntegrated; ++i)
    {
      PosMomVector[i] = valArr[i];
    }
  }

  void SetCurveLength(double len){
    fDistanceAlongCurve = len;
  }

  double GetCurveLength(){
    return fDistanceAlongCurve;
  }

friend std::ostream&
          operator<<( std::ostream& os, const FieldTrack& fieldTrack)
          {
            os<< " ( ";
            os<< " X= "<< fieldTrack.PosMomVector[0]<<" "
                       << fieldTrack.PosMomVector[1]<<" "
                       << fieldTrack.PosMomVector[2]<<" "; //Position
            os<< " P= "<< fieldTrack.PosMomVector[3]<<" "
                       << fieldTrack.PosMomVector[4]<<" "
                       << fieldTrack.PosMomVector[5]<<" "; //Momentum
            os<< " ) ";

            return os;
          }

};

FieldTrack::FieldTrack( const vecgeom::Vector3D<double>& pPosition,
                        const vecgeom::Vector3D<double>& pMomentum)
{
  for (int i = 0; i < 3; ++i)
    {
      PosMomVector[i]   = pPosition[i];
      PosMomVector[i+3] = pMomentum[i];
    }  
}

inline
FieldTrack::FieldTrack( const FieldTrack& rStVec)
{
  for (int i = 0; i < 6; ++i)
  {
    PosMomVector[i] = rStVec.PosMomVector[i];
  }
}

inline
FieldTrack & FieldTrack::operator = ( const FieldTrack& rStVec)
{
  if (&rStVec == this) return *this;
  
  for (int i = 0; i < 6; ++i)
  {
    PosMomVector[i] = rStVec.PosMomVector[i];
  }
  return *this;

}

#endif 