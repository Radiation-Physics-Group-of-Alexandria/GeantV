//
//  Compare the output of different steppers
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

#include "TemplateTUniformMagField.h"
#include "TemplateTMagFieldEquation.h"
#include "TemplateFieldEquationFactory.h"
#include "TemplateGUVIntegrationStepper.h"
#include "TemplateGUTCashKarpRKF45.h"

#include "TemplateGUIntegrationDriver.h"
#include "FieldTrack.h"

using namespace std;

// #define DEBUGAnanya

int main(int argc, char *args[])
{
    constexpr unsigned int Nposmom= 6; // Position 3-vec + Momentum 3-vec

    using Backend = vecgeom::kVc ;
    typedef typename Backend::precision_v Double;
    typedef vecgeom::Vector3D<Double> ThreeVectorSimd;
    typedef vecgeom::Vector3D<double> ThreeVector_d;

    using  GvEquationType=  TemplateTMagFieldEquation<Backend, TemplateTUniformMagField<Backend>, Nposmom>;
   
    /* -----------------------------SETTINGS-------------------------------- */
    
    /* Parameters of test
     - Modify values  */
    
    int no_of_steps = 250;         // No. of Steps for the stepper
    int stepper_no =  5;         // Choose stepper no., for refernce see above
    double step_len_mm = 200.;    // meant as millimeter;  //Step length 
    double z_field_in = DBL_MAX;
    
    //Checking for command line values :
    if(argc>1)
        stepper_no = atoi(args[1]);
    if(argc > 2)
       step_len_mm = (float)(stof(args[2]));   // *mm);
    if(argc > 3)
        no_of_steps = atoi(args[3]);
    if(argc > 4)
       z_field_in = (float) (stof(args[4]));     // tesla
    double step_len = step_len_mm * fieldUnits::millimeter;
    
    //Set Charge etc.
    double particleCharge = +1.0;      // in e+ units
    
    //Choice of output coordinates
    int
    columns[] =
    {
       1 , 1 , 1 ,  // position  x, y, z 
       1 , 1 , 1 ,  // momentum  x, y, z
       0 , 0 , 0 ,  // dydx pos  x, y, z
       1 , 1 , 0    // dydx mom  x, y, z
    }; //Variables in yOut[] & dydx[] we want to display - 0 for No, 1 for yes

    bool printErr= 0;   // Print the predicted Error
    bool printSep = 0;  // separator  '|'
    bool printInp = 0;  // print the input values
    
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

     #ifdef DEBUGAnanya
      cout<<"----Just before making TemplateTUniformMagField"<<endl;
     #endif 

    // Field
    auto gvField= new TemplateTUniformMagField<Backend>( fieldUnits::tesla * ThreeVector_d(x_field, y_field, z_field) );
    #ifdef DEBUGAnanya
     cout<<"----TemplateTUniformMagField Object constructed"<<endl;
    #endif
    cout << "#  Initial  Field strength (GeantV) = "
         << x_field << " , " << y_field << " , " << z_field 
       // << (1.0/fieldUnits::tesla) * gvField->GetValue()->X() << ",  "
         << " Tesla " << endl;
    cout << "#  Initial  momentum * c = " << x_mom << " , " << y_mom << " , " << z_mom << " GeV " << endl;
    //Create an Equation :
    #ifdef DEBUGAnanya
      cout<<"----Just before making EquationFactory"<<endl;
    #endif 
    auto gvEquation =
       TemplateFieldEquationFactory<Backend>::CreateMagEquation<TemplateTUniformMagField<Backend> >(gvField);
    #ifdef DEBUGAnanya
       cout<<"----EquationFactory made "<<endl;
    #endif 

    // gvEquation->InitializeCharge( particleCharge );  // Send it via Stepper instead    

    /*-------------------------PREPARING STEPPER-----------------------------*/
    
    //Create a stepper :

  #ifdef DEBUGAnanya
     cout<<"---- "<<endl;
     cout<<"---- Making TemplateGUTCashKarpRKF45"<<endl;
  #endif   
    // TemplateGUTCashKarpRKF45<Backend,GvEquationType,Nposmom> myStepper2(gvEquation);
  #ifdef DEBUGAnanya
    cout<<"---- constructed TemplateGUTCashKarpRKF45"<<endl;
  #endif
/*    TemplateGUVIntegrationStepper<Backend> *myStepper;
    myStepper = &myStepper2;*/

    TemplateGUVIntegrationStepper<Backend> *myStepper = new TemplateGUTCashKarpRKF45<Backend,GvEquationType,Nposmom>(gvEquation);

    myStepper->InitializeCharge( particleCharge );

    //Initialising coordinates
    const double mmGVf = fieldUnits::millimeter;
    const double ppGVf = fieldUnits::GeV ;  //   it is really  momentum * c_light
                                         //   Else it must be divided by fieldUnits::c_light;
    // const double ppGVf = fieldUnits::GeV / Constants::c_light;     // OLD

    Double yIn[] = {x_pos * mmGVf, y_pos * mmGVf ,z_pos * mmGVf,
                    x_mom * ppGVf ,y_mom * ppGVf ,z_mom * ppGVf};

/*    Double X_pos, Y_pos, Z_pos, X_mom, Y_mom, Z_mom;

    for (int i = 0; i < 4; ++i)
    {
      X_pos[i] = i-1.;
      Y_pos[i] = i-1.;
      Z_pos[i] = i-1.;
      X_mom[i] = i-1.;
      Y_mom[i] = i+1.-1.;
      Z_mom[i] = i+1.-1.;

    }
    cout<<"New X position is: "<<X_pos<<endl;

    Double yIn[] = {X_pos * mmGVf, Y_pos * mmGVf ,Z_pos * mmGVf,
                    X_mom * ppGVf ,Y_mom * ppGVf ,Z_mom * ppGVf};*/


    #ifdef DEBUGAnanya
      cout<<yIn[0]<<endl;
    #endif 
       

    auto gvEquation2 = new GvEquationType(gvField);
                   // new TMagFieldEquation<TUniformMagField, Nposmom>(gvField);
    // gvEquation2->InitializeCharge( particleCharge ); // Let's make sure
    
    // Should be able to share the Equation -- eventually
    // For now, it checks that it was Done() -- and fails an assert


    std::cout << "# step_len_mm = " << step_len_mm;
    
    
    //Empty buckets for results
    Double dydx[8] = {0.,0.,0.,0.,0.,0.,0.,0.},  // 2 extra safety buffer
           yout[8] = {0.,0.,0.,0.,0.,0.,0.,0.},
           yerr[8] = {0.,0.,0.,0.,0.,0.,0.,0.};
    
    /*-----------------------END PREPARING STEPPER---------------------------*/


    //=======================Test part for Integration driver====================
    auto testDriver = new TemplateGUIntegrationDriver<Backend>(0.2, myStepper);

    const ThreeVectorSimd  startPosition( yIn[0], yIn[1], yIn[2]);
    const ThreeVectorSimd  startMomentum( yIn[3], yIn[4], yIn[5]);
    TemplateGUFieldTrack<Backend>  yTrackIn(  startPosition, startMomentum );  // yStart
    TemplateGUFieldTrack<Backend>  yTrackOut( startPosition, startMomentum );  // yStart
    Double total_step = 0.;

    typedef typename Backend::bool_v Bool;
    Bool goodAdvance(true);
    double epsTol = 1.0e-5;

    // goodAdvance = testDriver->AccurateAdvance( yTrackIn, total_step, epsTol, yTrackOut );

    int nTracks = 16;
    FieldTrack yInput[nTracks], yOutput[nTracks];
    double hstep[nTracks];
    bool   succeeded[nTracks];
    testDriver->AccurateAdvance( yInput, hstep, epsTol, yOutput, nTracks, succeeded );

    cout<<"\n----Given hstep was: "<<endl;
    for (int i = 0; i < nTracks; ++i)
    {
      cout<<hstep[i]<<" ";
    }
    cout<<endl;



    //========================End testing IntegrationDriver=======================


    /*------ Clean up ------*/
    myStepper->InformDone(); 
    delete myStepper; 
    delete gvField;

    // deleting IntegrationDriver
    delete testDriver;
    
    cout<<"\n\n#-------------End of output-----------------\n";
}
