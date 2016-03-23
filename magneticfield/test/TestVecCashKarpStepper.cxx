//
//  Check if VecGUTCashKarpRKF45 has been designed well
// 
//  Based on testStepperFixed.cc
//    was the work of Somnath Banerjee in GSoC 2015
//
#include <iomanip>


#include "Units.h"

// using fieldUnits::meter;
using fieldUnits::millimeter;   
using fieldUnits::second;  
using fieldUnits::eplus;  
using fieldUnits::tesla;
using fieldUnits::degree;

#include <Vc/Vc>
#include "base/Vector3D.h"

#include "VecTUniformMagField.h"
#include "VecTMagFieldEquation.h"
#include "VecGUTCashKarpRKF45.h"

using namespace std;

#define DEBUGAnanya

int main()
{
  constexpr unsigned int Nposmom= 6; // Position 3-vec + Momentum 3-vec

  using Backend = vecgeom::kVc ;
  typedef typename Backend::precision_v Double_v;
  // typedef vecgeom::Vector3D<Double_v> ThreeVectorSimd;
  typedef vecgeom::Vector3D<double> ThreeVector_d;

  using  GvEquationType=  VecTMagFieldEquation<VecTUniformMagField, Nposmom>;
 
  /* -----------------------------SETTINGS-------------------------------- */
  
  /* Parameters of test
   - Modify values  */
  
  int no_of_steps = 250;         // No. of Steps for the stepper
  // int stepper_no =  5;         // Choose stepper no., for refernce see above
  double step_len_mm = 200.;    // meant as millimeter;  //Step length 
  double z_field_in = DBL_MAX;
  
  double step_len = step_len_mm * fieldUnits::millimeter;
  
  //Choice of output coordinates
  int
  columns[] =
  {
     1 , 1 , 1 ,  // position  x, y, z 
     1 , 1 , 1 ,  // momentum  x, y, z
     0 , 0 , 0 ,  // dydx pos  x, y, z
     1 , 1 , 0    // dydx mom  x, y, z
  }; //Variables in yOut[] & dydx[] we want to display - 0 for No, 1 for yes
  
  const unsigned int nwdf= 12;  // Width for float/double
  
  //Set coordinates here
  double
     x_pos = 0.,                 //pos - position  : input unit = mm
     y_pos = 0.,
     z_pos = 0.;
  double   
     x_mom = 0.,                 //mom - momentum  : input unit = GeV / c
     y_mom = 1.,
     z_mom = 1.;
  double          
     x_field = 0.,               //Uniform Magnetic Field (x,y,z)
     y_field = 0.;               //  Tesla // *tesla ;
  double z_field;
  // vecgeom::MaskedAssign(z_field_in < DBL_MAX, z_field_in, &z_field);

  if( z_field_in < DBL_MAX )
     z_field = z_field_in;
  else
     z_field = -1.0;  //  Tesla // *tesla ;

  // Field
  auto gvField= new VecTUniformMagField( fieldUnits::tesla * ThreeVector_d(x_field, y_field, z_field) );

  cout << "#  Initial  Field strength (GeantV) = "
       << x_field << " , " << y_field << " , " << z_field 
     // << (1.0/fieldUnits::tesla) * gvField->GetValue()->X() << ",  "
       << " Tesla " << endl;
  cout << "#  Initial  momentum * c = " << x_mom << " , " << y_mom << " , " << z_mom << " GeV " << endl;
  //Create an Equation :
 
  auto gvEquation = new GvEquationType(gvField);
     // TemplateFieldEquationFactory<Backend>::CreateMagEquation<VecTUniformMagField>(gvField);

  /*-------------------------PREPARING STEPPER-----------------------------*/
  
  //Create a stepper :

  auto myStepper = new VecGUTCashKarpRKF45<GvEquationType,Nposmom>(gvEquation);

  // myStepper->InitializeCharge( particleCharge );

  //Initialising coordinates
  const double mmGVf = fieldUnits::millimeter;
  const double ppGVf = fieldUnits::GeV ;  //   it is really  momentum * c_light
                                       //   Else it must be divided by fieldUnits::c_light;
  // const double ppGVf = fieldUnits::GeV / Constants::c_light;     // OLD

  Double_v yIn[] = {x_pos * mmGVf, y_pos * mmGVf ,z_pos * mmGVf,
                  x_mom * ppGVf ,y_mom * ppGVf ,z_mom * ppGVf};

  #ifdef DEBUGAnanya
    cout<<yIn[0]<<endl;
  #endif 
  
  std::cout << "# step_len_mm = " << step_len_mm;
  
  
  //Empty buckets for results
  Double_v dydx[8] = {0.,0.,0.,0.,0.,0.,0.,0.},  // 2 extra safety buffer
         yout[8] = {0.,0.,0.,0.,0.,0.,0.,0.},
         yerr[8] = {0.,0.,0.,0.,0.,0.,0.,0.};
  
  /*-----------------------END PREPARING STEPPER---------------------------*/




  /*---------------------------------------------------*/
  //        -> First Print the (commented) title header
  cout<<"\n#";
  cout<<setw(5)<<"Step";
  for (int i=0; i<6;i++)
      if (columns[i])
      {
         cout << setw(nwdf-2)<< "yOut[" << i << "]";
      }
  for (int i=0; i<6;i++)
      if (columns[i+6])
      {
         cout << setw(nwdf-2)<< "dydx[" << i << "]";
      }    
  cout<<setw(nwdf)<<"tan-1(y/x)";
  cout<<"\n";
  
  // Units
  const char *nameUnitLength= "mm";
  const char *nameUnitMomentum= "GeV/c";
  cout<<setw(6)<<"#Numbr";
  for (int i=0; i<6;i++)
      if (columns[i])
      {
         const char* nameUnit = ( i<3 ) ? nameUnitLength : nameUnitMomentum ; 
         cout << setw(nwdf)<< nameUnit;
      }    
  cout<<"\n";
  

  cout.precision(3);
  
  /*----------------NOW STEPPING-----------------*/
  no_of_steps = 1;
  for(int j=0; j<no_of_steps; j++)
  {
      cout<<setw(6)<<j ;           //Printing Step number
      Double_v charge(-1.);

      myStepper->RightHandSideVIS<Backend>(yIn, charge, dydx);               //compute dydx - to supply the stepper

      if( j > 0 )  // Let's print the initial points!
      {
         myStepper->StepWithErrorEstimate<Backend>(yIn,dydx,step_len,yout,yerr);   //Call the 'trial' stepper

      }
      //-> Then print the data
      cout.setf (ios_base::fixed);
      cout.precision(4);
      for(int i=0; i<3;i++)
          if(columns[i]){
             cout<<setw(nwdf)<< yout[i] / mmGVf ;
             cout.precision(3);
          }

      cout.unsetf (ios_base::fixed);        
      cout.setf (ios_base::scientific);
      for(int i=3; i<6;i++)
          if(columns[i]){
             cout<<setw(nwdf) << yout[i] / ppGVf ;
          }
      cout.unsetf (ios_base::scientific);
      
      for(int i=0; i<6;i++)   // Print auxiliary components
      {
         double unitGVf=1;  
         // double unitRef=1;
         // if( i < 3 )             // length / length

         if( i >= 3 ){
            unitGVf = ppGVf / mmGVf; //  dp / ds
         }
         
         if( i == 0 ){            // length / length                      
            // cout.unsetf (ios_base::scientific);
            cout.setf (ios_base::fixed);
         }else if( i == 3 ){
            cout.unsetf (ios_base::fixed);              
            cout.setf (ios_base::scientific);
         }
         
         if(columns[i+6])
         {
             cout<<setw(nwdf)<< dydx[i] / unitGVf ;
         }
      }
      cout.unsetf(ios_base::scientific);
      if( j > 0 )  // Step 0 did not move -- printed the starting values
      {
         // cout.unsetf(ios_base::scientific);
         cout.setf (ios_base::fixed);                         
         cout.precision(2);
         cout<<setw(nwdf) << atan2(yout[1],yout[0])/degree;
         
         //Copy yout into yIn
         for(int i=0; i<6;i++){
            yIn[i] = yout[i];

         }
      }
      
      cout<<"\n";
  }
  
  /*-----------------END-STEPPING------------------*/

  /*------ Clean up ------*/
  myStepper->InformDone(); 

  delete myStepper;
  delete gvField;
  
  cout<<"\n\n#-------------End of output-----------------\n";
}
