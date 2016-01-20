
#include "base/Vector3D.h"

#include "GUVEquationOfMotion.h"

#include "GUVVectorEquationOfMotion.h"
#include "GUVVectorField.h"
#include "TVectorMagFieldEquation.h"
#include "GUVMagneticField.h"
#include "GUVVectorMagneticField.h"
#include "TVectorUniformMagField.h"

#include "GUVField.h"
#include "TMagFieldEquation.h"
#include "TUniformMagField.h"

using namespace std;

typedef typename vecgeom::kVc::precision_v Double_v;
typedef typename vecgeom::kVcFloat::precision_v Float_v;
typedef typename vecgeom::kVc::int_v Int_v;
typedef typename vecgeom::kVc::bool_v Bool_v;

using ThreeVector_f = vecgeom::Vector3D<float>;
using ThreeVector_d = vecgeom::Vector3D<double>;

using ThreeVectorSimd_f = vecgeom::Vector3D<Float_v>;
using ThreeVectorSimd_d = vecgeom::Vector3D<Double_v>;


/*int main(){

  ThreeVectorSimd_f FieldValue2(0.0, 0.0, 1.0);

  TVectorUniformMagField* ConstBfield = new TVectorUniformMagField( FieldValue2 );
  cout<<"here 0.5"<<endl;
  using EquationType = TVectorMagFieldEquation<TVectorUniformMagField, gNposmom>;

  return 0;
}*/

GUVVectorEquationOfMotion*  CreateFieldAndEquation(ThreeVectorSimd_f constFieldValue);
bool  TestEquation(GUVVectorEquationOfMotion* );

constexpr unsigned int gNposmom= 6; // Position 3-vec + Momentum 3-vec

ThreeVectorSimd_f  FieldValue(0.0, 0.0, 1.0);

/*int
main( int, char** )
{
  std::cout<<"Entered main"<<std::endl;

  GUVVectorEquationOfMotion* eq = CreateFieldAndEquation( FieldValue );

  std::cout<<"eq created "<<std::endl;
  TestEquation(eq);

  return 1;
}*/
int main(int argc, char *argv[]){
  int flag = atoi(argv[1]);

  if(flag) {

    ThreeVectorSimd_f FieldValue2(0.0, 0.0, 1.0);

    TVectorUniformMagField* ConstBfield = new TVectorUniformMagField( FieldValue2 );
    cout<<"here 0.5"<<endl;
    // using EquationType = TVectorMagFieldEquation<TVectorUniformMagField, gNposmom>;
  }

  else{
    std::cout<<"Entered main"<<std::endl;

    GUVVectorEquationOfMotion* eq = CreateFieldAndEquation( FieldValue );

    std::cout<<"eq created "<<std::endl;
    TestEquation(eq);
  }

  return 0;
}


GUVVectorEquationOfMotion* CreateFieldAndEquation(ThreeVectorSimd_f FieldValue)
{
  std::cout<<"Creating field and equation..."<<std::endl;

  TVectorUniformMagField* ConstBfield = new TVectorUniformMagField( FieldValue );
  cout<<"here 0.5"<<endl;
  using EquationType = TVectorMagFieldEquation<TVectorUniformMagField, gNposmom>;
  cout<<"here1"<<endl;
  // GUVVectorEquationOfMotion*  magEquation= new TVectorMagFieldEquation<TUniformMagField, gNposmom>(ConstBfield);
  GUVVectorEquationOfMotion*  magEquation = new EquationType(ConstBfield);
  
  // pField = ConstBfield; 
  return magEquation;
}

int gVerbose= 1;

bool TestEquation(GUVVectorEquationOfMotion* equation)
{
  constexpr double perMillion = 1e-6;
  bool   hasError = false;  // Return value
  
  ThreeVectorSimd_d PositionVec( 1., 2.,  3.);  // initial
  ThreeVectorSimd_d MomentumVec( 0., 0.1, 1.);
  ThreeVectorSimd_f FieldVec( 0., 0., 1.);  // Magnetic field value (constant)

  Double_v PositionTime[4] = { PositionVec.x(), PositionVec.y(), PositionVec.z(), 0.0};

  Double_v dydx[gNposmom];
  Double_v PositionMomentum[gNposmom];

  double charge= -1;  
  Double_v Charge(-1);

  PositionMomentum[0] = PositionVec[0];
  PositionMomentum[1] = PositionVec[1];
  PositionMomentum[2] = PositionVec[2];
  PositionMomentum[3] = MomentumVec[0];
  PositionMomentum[4] = MomentumVec[1];
  PositionMomentum[5] = MomentumVec[2];

  
  equation->InitializeCharge( charge );

  equation->EvaluateRhsGivenB( PositionTime, FieldVec,  Charge,   dydx );

  ThreeVectorSimd_d  ForceVec( dydx[3], dydx[4], dydx[5]);

  // Check result
  Double_v MdotF = MomentumVec.Dot(ForceVec);
  ThreeVectorSimd_d doubleFieldVec;
  for(int i=0; i<3;i++) doubleFieldVec[i] = (Double_v) FieldVec[i];

  Double_v BdotF = doubleFieldVec.Dot(ForceVec);

  Double_v momentumMag = MomentumVec.Mag();
  Double_v fieldMag =   doubleFieldVec.Mag();
  Double_v ForceMag =   ForceVec.Mag();

  bool cond1 = !( vecgeom::IsEmpty( ForceMag != momentumMag*fieldMag ) );
  bool cond2 = !( vecgeom::IsEmpty( Vc::abs(MdotF) > perMillion * MomentumVec.Mag()    * ForceVec.Mag() ) );
  bool cond3 = !( vecgeom::IsEmpty( Vc::abs(BdotF) > perMillion * doubleFieldVec.Mag() * ForceVec.Mag() ) );


  if( cond1 ) {
     std::cerr << "ERROR: Force magnitude is not equal to momentum * field."  << std::endl;     
  }
     
  assert( cond1 );  // Must add coefficient !!
  
  if( cond2 )
  { 
     std::cerr << "ERROR: Force due to magnetic field is not perpendicular to momentum!!"  << std::endl;
     hasError= true;
  }
  else if ( gVerbose )
  {
     std::cout << " Success:  Good (near zero) dot product momentum . force " << std::endl;
  }

  if( cond3 )
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
  //return true;
}
