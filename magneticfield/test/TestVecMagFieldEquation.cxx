//
//
#include "base/Vector3D.h"

#include "VecGUVMagneticField.h"
#include "VecGUVField.h"
#include "VecTMagFieldEquation.h"
#include "VecTUniformMagField.h"

//To be able to call vecgeom::kScalar/kVc/kCuda
#include <Vc/Vc>
// #include "backend/vc/Backend.h"
// #include "backend/vcfloat/Backend.h"
// #include "backend/scalarfloat/Backend.h"


using ThreeVector_d = vecgeom::Vector3D<double>;

typedef vecgeom::Vector3D<typename vecgeom::kVc::precision_v> ThreeVectorSimd_d;

template <class T_Field, unsigned int Size>
VecTMagFieldEquation<T_Field,Size>*  CreateFieldAndEquation(ThreeVector_d constFieldValue);

template <class Backend, class T_Field, unsigned int Size>
bool TestEquation(VecTMagFieldEquation<T_Field,Size>* );

constexpr unsigned int gNposmom= 6; // Position 3-vec + Momentum 3-vec
ThreeVector_d  FieldValue(0.0, 0.0, 1.0);
int gVerbose= 1;

int
main( int, char** )
{
   
  auto eq = CreateFieldAndEquation<VecTUniformMagField, gNposmom>( FieldValue );

  TestEquation<vecgeom::kVc, VecTUniformMagField, gNposmom>(eq);

  return 1;
}

template <class Type>
Type CustomAbs(Type T ){
  return vecgeom::VECGEOM_IMPL_NAMESPACE::Abs(T);
}

template <class T_Field, unsigned int Size>
VecTMagFieldEquation<T_Field,Size> * CreateFieldAndEquation(ThreeVector_d FieldValue)
{
  VecTUniformMagField*   ConstBfield = new VecTUniformMagField( FieldValue );
  
  using EquationType = VecTMagFieldEquation <T_Field, Size>;

  VecTMagFieldEquation<T_Field,Size>* magEquation = new EquationType(ConstBfield);
  
  return magEquation;
}



template <class Backend, class T_Field, unsigned int Size>
bool TestEquation(VecTMagFieldEquation<T_Field,Size>* equation)
{
  constexpr double perMillion = 1e-6;
  bool   hasError = false;  // Return value
  
  // vecgeom::Vector3D<typename Backend::precision_v> 
  vecgeom::Vector3D<typename Backend::precision_v> PositionVec( 1., 2.,  3.);  // initial
  vecgeom::Vector3D<typename Backend::precision_v> MomentumVec( 0., 0.1, 1.);
  vecgeom::Vector3D<typename Backend::precision_v> FieldVec   ( 0., 0.,  1.);  // Magnetic field value (constant)

  typedef typename Backend::precision_v Double;
  typedef typename Backend::bool_v Bool_v;


  Double dydx[gNposmom];
  Double PositionMomentum[gNposmom];

  // double charge= -1;  
  Double Charge(-1);

  PositionMomentum[0] = PositionVec[0];
  PositionMomentum[1] = PositionVec[1];
  PositionMomentum[2] = PositionVec[2];
  PositionMomentum[3] = MomentumVec[0];
  PositionMomentum[4] = MomentumVec[1];
  PositionMomentum[5] = MomentumVec[2];

  
  // equation->InitializeCharge( charge );
  equation->InformReady();
  equation->template TEvaluateRhsGivenB<Backend>( PositionMomentum, FieldVec, Charge, dydx );

  vecgeom::Vector3D<typename Backend::precision_v>  ForceVec( dydx[3], dydx[4], dydx[5]);

  // Check result
  Double MdotF = MomentumVec.Dot(ForceVec);
  vecgeom::Vector3D<typename Backend::precision_v>  doubleFieldVec;
  // for(int i=0; i<3;i++) doubleFieldVec[i] = (Double) FieldVec[i];

  Double BdotF = FieldVec.Dot(ForceVec);

  Double momentumMag = MomentumVec.Mag();
  Double fieldMag =   FieldVec.Mag();
  Double ForceMag =   ForceVec.Mag();

  cout<<"\n";
  cout<<" momentumMag: "<<momentumMag<<endl;
  cout<<" fieldMag:    "<<fieldMag<<endl;
  cout<<" ForceMag:    "<<ForceMag<<endl;
  cout<<" abs(BdotF):  "<<CustomAbs(BdotF)<<endl;
  cout<<" abs(MdotF):  "<<CustomAbs(MdotF)<<endl;

  bool cond1 = !( vecgeom::IsEmpty( ForceMag != momentumMag*fieldMag ) );
  bool cond2 = !( vecgeom::IsEmpty( CustomAbs(MdotF) > perMillion * MomentumVec.Mag()    * ForceVec.Mag() ) );
  bool cond3 = !( vecgeom::IsEmpty( CustomAbs(BdotF) > perMillion * FieldVec.Mag() * ForceVec.Mag() ) );

  Bool_v print1 = (ForceMag != momentumMag*fieldMag) ;
  Bool_v print2 = CustomAbs(MdotF) > perMillion * MomentumVec.Mag()    * ForceVec.Mag();
  Bool_v print3 = CustomAbs(BdotF) > perMillion * FieldVec.Mag() * ForceVec.Mag();

  cout<<"\nPrinting Bool_v "<<endl;
  cout<< print1 <<" and "<<cond1<<endl;
  cout<< print2 <<" and "<<cond2<<endl;
  cout<< print3 <<" and "<<cond3<<endl;
  cout<<"\n"<<endl;


  if( cond1 ) {
     std::cerr << "ERROR: Force magnitude is not equal to momentum * field."  << std::endl;  
     std::cout<< ForceMag<<" vs "<< momentumMag*fieldMag<<"\n"<<endl;   
  }
     
  assert( cond1 );  // Must add coefficient !!
  
  if( cond2 )
  { 
     std::cerr << "ERROR: Force due to magnetic field is not perpendicular to momentum!!"  << std::endl;
     hasError= true;
     cout<<CustomAbs(MdotF)<<" vs "<< perMillion * MomentumVec.Mag()    * ForceVec.Mag() <<"\n"<<endl;
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

    cout<<"\n"<<CustomAbs(BdotF)<<" vs "<< perMillion * FieldVec.Mag()  * ForceVec.Mag() <<endl;
  }
  else if ( gVerbose )
  {
    std::cout << " Success:  Good (near zero) dot product magnetic-field . force " << std::endl;
  }

  return hasError;
  //return true;
}
