//
//
#include "base/Vector3D.h"

#include "GUVEquationOfMotion.h"

#include "TMagFieldEquation.h"
#include "TUniformMagField.h"

using ThreeVector_f = vecgeom::Vector3D<float>;
using ThreeVector_d = vecgeom::Vector3D<double>;

GUVEquationOfMotion*  CreateFieldAndEquation(ThreeVector_f constFieldValue);
bool  TestEquation(GUVEquationOfMotion* );

constexpr unsigned int gNposmom= 6; // Position 3-vec + Momentum 3-vec

ThreeVector_f  FieldValue(0.0, 0.0, 1.0);

int
main( int, char** )
{
  GUVEquationOfMotion* eq = CreateFieldAndEquation( FieldValue );
  TestEquation(eq);

  return 1;
}

// GUVField*       pField;

GUVEquationOfMotion* CreateFieldAndEquation(ThreeVector_f FieldValue)
{
  TUniformMagField*     ConstBfield = new TUniformMagField( FieldValue );
  // typedef typename TMagFieldEquation<TUniformMagField, gNposmom>  EquationType;
  using EquationType = TMagFieldEquation<TUniformMagField, gNposmom>;
  
  // GUVEquationOfMotion*  magEquation= new TMagFieldEquation<TUniformMagField, gNposmom>(ConstBfield);
  GUVEquationOfMotion*  magEquation = new EquationType(ConstBfield);
  
  // pField = ConstBfield; 
  return magEquation;
}

int gVerbose= 1;

bool TestEquation(GUVEquationOfMotion* equation)
{
  constexpr double perMillion = 1e-6;
  bool   hasError = false;  // Return value
  
  ThreeVector_d PositionVec( 1., 2.,  3.);  // initial
  ThreeVector_d MomentumVec( 0., 0.1, 1.);
  ThreeVector_f FieldVec( 0., 0., 1.);  // Magnetic field value (constant)

  double PositionTime[4] = { PositionVec.x(), PositionVec.y(), PositionVec.z(), 0.0};
  // double PositionTime[4] = { 1., 2., 3., 4.};
  PositionTime[0] = PositionVec.x();
  PositionTime[1] = PositionVec.y();
  PositionTime[2] = PositionVec.z();

  // double magField[3];

  double dydx[gNposmom];
  double PositionMomentum[gNposmom];

  double charge= -1;

  PositionMomentum[0] = PositionVec[0];
  PositionMomentum[1] = PositionVec[1];
  PositionMomentum[2] = PositionVec[2];
  PositionMomentum[3] = MomentumVec[0];
  PositionMomentum[4] = MomentumVec[1];
  PositionMomentum[5] = MomentumVec[2];

  // double FieldArr[3]= { FieldVec.x(), FieldVec.y(), FieldVec.z() };
  
  equation->InitializeCharge( charge );
  equation->EvaluateRhsGivenB( PositionTime, FieldVec, /* charge, */ dydx );

  ThreeVector_d  ForceVec( dydx[3], dydx[4], dydx[5]);

  // Check result
  double MdotF = MomentumVec.Dot(ForceVec);
  double BdotF = FieldVec.Dot(ForceVec);

  double momentumMag = MomentumVec.Mag();
  double fieldMag =   FieldVec.Mag();
  double ForceMag =   ForceVec.Mag();

  if( ForceMag != momentumMag * fieldMag ) {
     std::cerr << "ERROR: Force magnitude is not equal to momentum * field."  << std::endl;     
  }
     
  assert( ForceMag != momentumMag * fieldMag );  // Must add coefficient !!
  
  if( std::fabs(MdotF) > perMillion * MomentumVec.Mag() * ForceVec.Mag() )
  { 
     std::cerr << "ERROR: Force due to magnetic field is not perpendicular to momentum!!"  << std::endl;
     hasError= true;
  }
  else if ( gVerbose )
  {
     std::cout << " Success:  Good (near zero) dot product momentum . force " << std::endl;
  }
  if( std::fabs(BdotF) > perMillion * FieldVec.Mag() * ForceVec.Mag() )
  { 
    std::cerr << "ERROR: Force due to magnetic field is not perpendicular to B field!"
              << std::endl; 
    std::cerr << " Vectors:  BField   Force " << std::endl;
    for ( int i = 0; i < 3; i ++ )
      std::cerr << "   [" << i << "] " << FieldVec[i] << " " << ForceVec[i] << std::endl;

    hasError= true;
  }
  else if ( gVerbose )
  {
    std::cout << " Success:  Good (near zero) dot product magnetic-field . force " << std::endl;
  }

  return hasError;
}
