//
//  Test GUIntegrationDriver 
//   * compares with the output of a reference stepper (high accuracy) - ok for small steps
// 
//  Based on testStepperFixed.cc
//    which was started from the work of Somnath Banerjee in GSoC 2015
//

#include <iomanip>

// #include "G4UniformMagField.hh"
// #include "G4SystemOfUnits.hh"
#include "Units.h"

// using fieldUnits::meter;
using fieldUnits::millimeter;   
using fieldUnits::second;  
using fieldUnits::eplus;  
using fieldUnits::tesla;
using fieldUnits::degree;

// #include "G4ios.hh"

#include "ThreeVector.h"

#include "TUniformMagField.h"

#include "TMagFieldEquation.h"

#include "GUVIntegrationStepper.h"
#include "TClassicalRK4.h"
#include "GUTCashKarpRKF45.h"

#include "TSimpleRunge.h"
#include "GUExactHelixStepper.h"
// #include "GULineSection.hh"

// #include "BogackiShampine23.hh"
// #include "DormandPrince745.hh"
// #include "BogackiShampine45.h"
// #include "G4SimpleHeum.hh"

#include "GUFieldTrack.h"
#include "GUIntegrationDriver.h"

// #define  COMPARE_TO_G4  1

#ifdef COMPARE_TO_G4
//#include <system.h>
//#include "G4Types.h"
#include "G4UniformMagField.hh"
#include "G4ExactHelixStepper.hh"
#include "G4ClassicalRK4.hh"
#endif

using namespace std;
// using namespace CLHEP;

//Version 2.0 - Includes direct comparison with ExactHelix

/* Stepper No.
 0. GUExactHelix
 4. TClassicalRK4
 5. TCashKarpRKF45 

Potential expansion:
 2. G4SimpleHeum
 3. BogackiShampine23
 6. BogackiShampine45
 7. DormandPrince745
 */

const unsigned int Nposmom= 6; // Position 3-vec + Momentum 3-vec

int main(int argc, char *args[])
{
    using  GvEquationType=  TMagFieldEquation<TUniformMagField, Nposmom>;
   
    /* -----------------------------SETTINGS-------------------------------- */
    
    /* Parameters of test
     - Modify values  */
    
    int no_of_steps = 20;         // No. of Steps for the stepper
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
    double particleCharge = +1.0,      // in e+ units
       spin=0.0,                       // ignore the spin
       magneticMoment= 0.0,            // ignore the magnetic moment
       mass = 1;
    
    //Choice of output coordinates
    int
    columns[] =
    {
       1 , 1 , 1 ,  // position  x, y, z 
       0 , 0 , 0 ,  // momentum  x, y, z
       1 , 1 , 1 ,  // dydx pos  x, y, z
       1 , 1 , 1    // dydx mom  x, y, z
    }; //Variables in yOut[] & dydx[] we want to display - 0 for No, 1 for yes

    bool printErr= 0,   // Print the predicted Error
         printRef= 0,   // Print the reference Solution 
         printDiff= 1;  // Print the difference 
    bool printSep = 1;  // separator  '|'
    bool printOk  = 1;  // "Ok" or "Bad"
    bool printRatio = 0;  //  yDiff / yAverage
    bool printInp = 0;  // print the input values
    bool printInpX= 0;  // print the input values for Ref 

    const int    epsDigits = 5;       //  Digits - for printing
    const double errAmplif = 1.0e5;   //  Amplification for printing of error
    
    const unsigned int wd= 12;  // Width for float/double
    const unsigned int wdf = wd + epsDigits - 3;
    
    //Set coordinates here
    double
       x_pos = 0.,                 //pos - position  : input unit = mm
       y_pos = 0.,
       z_pos = 0.;
    double   
       x_mom = 0.0,                 //mom - momentum  : input unit = GeV / c
       y_mom = 1.0,
       z_mom = 1.0;
    double          
       x_field = 0.,               //Uniform Magnetic Field (x,y,z)
       y_field = 0.;               //  Tesla // *tesla ;
    double z_field;
    if( z_field_in < DBL_MAX )
       z_field = z_field_in;
    else
       z_field = -1.0;  //  Tesla // *tesla ;

    // Field
    auto gvField= new TUniformMagField( fieldUnits::tesla * ThreeVector(x_field, y_field, z_field) );

    cout << "#  Initial  Field strength (GeantV) = "
         << x_field << " , " << y_field << " , " << z_field 
       // << (1.0/fieldUnits::tesla) * gvField->GetValue()->X() << ",  "
         << " Tesla " << endl;
    cout << "#  Initial  momentum * c = " << x_mom << " , " << y_mom << " , " << z_mom << " GeV " << endl;
    //Create an Equation :
    auto gvEquation =
       new GvEquationType(gvField);
       // new TMagFieldEquation<TUniformMagField, Nposmom>(gvField);
    // gvEquation->InitializeCharge( particleCharge );

    /*-------------------------PREPARING STEPPER-----------------------------*/
    
    //Create a stepper :
    GUVIntegrationStepper *myStepper; // , *exactStepper;
    // G4MagIntegrationStepper *g4refStepper;    

    cout << "#  Chosen the   ";
    //Choose the stepper based on the command line argument
    switch(stepper_no){
       // case 0: myStepper = new GUExactHelixStepper(gvEquation);
       //  cout << " GUExactHelix stepper. " << endl;
       //  break;
      case 1: myStepper = new TSimpleRunge<GvEquationType,Nposmom>(gvEquation);
         cout << " TSimpleRunge stepper. " << endl;
         break;         
         // case 2: myStepper = new G4SimpleHeum(gvEquation);   break;
       // case 3: myStepper = new BogackiShampine23(gvEquation); break;
      case 4: myStepper = new TClassicalRK4<GvEquationType,Nposmom>(gvEquation);
         cout << " TClassicalRK4 stepper. " << endl;         
         break;
      case 5: myStepper = new GUTCashKarpRKF45<GvEquationType,Nposmom>(gvEquation);
         cout << " GUTCashKarpRKF45 stepper. " << endl;                  
         break;         
       // case 6: myStepper = new BogackiShampine45(gvEquation); break;
       // case 7: myStepper = new DormandPrince745(gvEquation);  break;
      default : myStepper = 0 ;
    }

    const double hminimum = 0.0001 * millimeter;  // Minimum step = 0.1 microns
    const double epsTol    = 1.0e-5;              // Relative error tolerance
    int   statisticsVerbosity= 1;

    auto integrDriver= new GUIntegrationDriver( hminimum,
                                                myStepper,
                                                Nposmom,
                                                statisticsVerbosity); 
    // myStepper->InitializeCharge( particleCharge );
    integrDriver->InitializeCharge( particleCharge );
    
    //Initialising coordinates
    const double mmGVf = fieldUnits::millimeter;
    const double ppGVf = fieldUnits::GeV ;  //   it is really  momentum * c_light
                                         //   Else it must be divided by fieldUnits::c_light;
    // const double ppGVf = fieldUnits::GeV / Constants::c_light;     // OLD

    // double yIn[] = {x_pos,y_pos,z_pos,x_mom,y_mom,z_mom};
    double yIn[] = {x_pos * mmGVf, y_pos * mmGVf ,z_pos * mmGVf,
                    x_mom * ppGVf ,y_mom * ppGVf ,z_mom * ppGVf};
    
#if COMPARE_TO_G4
    const double mmG4 = CLHEP::millimeter;
    const double ppG4 = CLHEP::GeV ;  //  In G4 too 'p' means p*c -- so no division  / CLHEP::c_light;

    const double mmRef = mmG4;  // Unit for reference of lenght   - milli-meter
    const double ppRef = ppG4;  // Unit for reference of momentum - GeV / c^2
    
    // G4double yInX[] = {x_pos * mmG4, y_pos * mmG4 ,z_pos * mmG4,
    //                   x_mom * ppG4 ,y_mom * ppG4 ,z_mom * ppG4};

    G4UniformMagField myField( CLHEP::tesla * G4ThreeVector(x_field, y_field, z_field));
    cout << " Field strength = " << (1.0/CLHEP::tesla) * myField << " Tesla " 
         << endl;
    
    G4Mag_UsualEqRhs *g4Equation = new G4Mag_UsualEqRhs(&myField);
    
    G4ChargeState chargeState(particleCharge,             // The charge can change (dynamic)
                              spin=0.0,
                              magneticMoment=0.0);
    
    g4Equation->SetChargeMomentumMass( chargeState,
                                       G4ThreeVector(x_mom, y_mom, z_mom).mag(), //momentum magnitude
                                       mass);  // unused
//  auto g4exactStepper= = new G4ExactHelixStepper(g4Equation);
    auto g4exactStepper= = new G4ClassicalRK4(g4Equation);
    
    auto exactStepper = g4ExactStepperGV;
#else
    // double yInX[] = {x_pos * mmGVf, y_pos * mmGVf ,z_pos * mmGVf,
    //                 x_mom * ppGVf ,y_mom * ppGVf ,z_mom * ppGVf};    
    const double mmRef = mmGVf; // Unit for reference of lenght   - milli-meter
    const double ppRef = ppGVf; // Unit for reference of momentum - GeV / c^2
    
    auto gvEquation2 = new GvEquationType(gvField);
                   // new TMagFieldEquation<TUniformMagField, Nposmom>(gvField);
    // gvEquation2->InitializeCharge( particleCharge ); // Let's make sure
    
    // Should be able to share the Equation -- eventually
    // For now, it checks that it was Done() -- and fails an assert

    //Creating the baseline stepper
    auto exactStepperGV =
        new GUTCashKarpRKF45<GvEquationType,Nposmom>(gvEquation2);       
        // new TClassicalRK4<GvEquationType,Nposmom>(gvEquation2);
    cout << "#  Reference stepper is: GUTCashKarpRKF45 <GvEquationType,Nposmom>(gvEquation2);" << endl;
    // cout << "#  Reference stepper is: TClassicalRK4 <GvEquationType,Nposmom>(gvEquation2);" << endl;    
       // new TSimpleRunge<GvEquationType,Nposmom>(gvEquation2);    
       // new GUExactHelixStepper(gvEquation2);

    // Configure Stepper for current particle
    exactStepperGV->InitializeCharge( particleCharge ); // Passes to Equation, is cached by stepper
    // gvEquation2->InitializeCharge( particleCharge ); //  Different way - in case this works
    
    auto exactStepper = exactStepperGV;
#endif
    std::cout << "# step_len_mm = " << step_len_mm;
    std::cout << " mmRef= " << mmRef << "   ppRef= " << ppRef << std::endl;
    
    double yInX[] = {x_pos * mmRef, y_pos * mmRef ,z_pos * mmRef,
                     x_mom * ppRef ,y_mom * ppRef ,z_mom * ppRef};

    double stepLengthRef = step_len_mm * mmRef;
    
    //Empty buckets for results
    double dydx[8] = {0.,0.,0.,0.,0.,0.,0.,0.},  // 2 extra safety buffer
        dydxRef[8] = {0.,0.,0.,0.,0.,0.,0.,0.},
          yout [8] = {0.,0.,0.,0.,0.,0.,0.,0.},
          youtX[8] = {0.,0.,0.,0.,0.,0.,0.,0.},
          yerr [8] = {0.,0.,0.,0.,0.,0.,0.,0.},
          yerrX[8] = {0.,0.,0.,0.,0.,0.,0.,0.},
          yDiff[8] = {0.,0.,0.,0.,0.,0.,0.,0.},
          yAver[8] = {0.,0.,0.,0.,0.,0.,0.,0.};          
    /*-----------------------END PREPARING STEPPER---------------------------*/


    /*---------------------------------------------------*/
    //        -> First Print the (commented) title header
    cout<<"\n#";
    cout<<setw(5)<<"Step";
    cout<<setw(8)<<" StepLen";    
    for (int i=0; i<6;i++)
        if (columns[i])
        {
           if( printSep ) cout << " | " ;  // Separator
           if( printInp )
              cout << setw(wd-2)<< "yIn[" << i << "]";
           if( printInp )
              cout << setw(wd-2)<< "yInX[" << i << "]";           
           cout << setw(wd-2)<< "yOut[" << i << "]";
           if( printRef )
              cout << setw(wd-2) << "yOutX[" << i << "]";
           if( printErr )
              cout << setw(wd-1) << "yErr[" << i << "]";
           if( printDiff )
              cout << setw(wd) << " yOut-yOutX[" << i << "]";
           if( printOk  ) cout << " Cmp";
           if( printRatio ) cout << setw(wd+2) << " Ratio ";           
        }
    for (int i=0; i<6;i++)
        if (columns[i+6])
        {
           if( printSep ) cout << " | " ;  // Separator
           if( printInp )                               
              cout << setw(wd-2)<< "yIn[" << i << "]";           
           cout << setw(wd-2)<< "dydx[" << i << "]";
           if( printRef )                    
              cout << setw(wd-2) << "dydxRef[" << i << "]";
           // if( printErr ) cout << setw(wd-2) << "yErr[" << i << "]";
           // if( printDiff ) cout << setw(wdf-2) << "dydx-Ref[" << i << "]";
        }    
    cout<<setw(wd)<<"tan-1(y/x)";
    cout<<"\n";
    
    // Units
    const char *nameUnitLength= "mm";
    const char *nameUnitMomentum= "GeV/c";
    cout<<setw(6)<<"#Numbr";
    cout<<setw(8)<<nameUnitLength;
    for (int i=0; i<6;i++)
        if (columns[i])
        {
           if( printSep ) cout << " | " ;  // Separator                      
           const char* nameUnit = ( i<3 ) ? nameUnitLength : nameUnitMomentum ; 
           cout << setw(wd)<< nameUnit;
           if( printRef )  cout << setw(wd) << nameUnit;  // Reference output
           if( printErr ){
              cout << setw(4) << nameUnit;  // Estim. error
              cout << " * " << setw(wd-7) << errAmplif << " ";
           }
           if( printDiff ) cout << setw(wdf) << nameUnit;    // Diff  new - ref
           if( printOk  )  cout << "    ";
           if( printRatio ) cout << setw(wd+2) << "/epsTol ";                      
        }    
    cout<<"\n";
    
    //-> Then print the data

    
    // cout.setf (ios_base::scientific);
    // cout.setf (ios_base::scientific);    
    cout.precision(3);

    ThreeVector  startPosition( yIn[0], yIn[1], yIn[2]);
    ThreeVector  startMomentum( yIn[3], yIn[4], yIn[5]);

    GUFieldTrack yStart( startPosition, startMomentum ); 
    double total_step =0;
    /*----------------NOW STEPPING-----------------*/
    
    for(int j=0; j<no_of_steps; j++)
    {
        cout<<setw(6)<<j ;           //Printing Step number

        // myStepper->RightHandSide(yIn, dydx);               //compute dydx - to supply the stepper
        exactStepper->RightHandSide(yInX, dydxRef);   //compute the value of dydx for the exact stepper

        // Driver begins at the start !
        GUFieldTrack  yTrack( startPosition, startMomentum );  // yStart
        
        if( j > 0 )  // Let's print the initial points!
        {
           total_step += step_len;
           
           bool goodAdvance=
              integrDriver->AccurateAdvance( yTrack, total_step, epsTol); // , hInitial );
           // *****************************
           
           // Output is in yTrack
           
           ThreeVector PositionOut = yTrack.GetPosition(); 
           yout[0]= PositionOut.x();
           yout[1]= PositionOut.y();
           yout[2]= PositionOut.z();           
           ThreeVector MomentumOut = yTrack.GetMomentum(); 
           yout[3]= MomentumOut.x();
           yout[4]= MomentumOut.y();
           yout[5]= MomentumOut.z();           
          
           // myStepper->StepWithErrorEstimate(yIn,dydx,step_len,yout,yerr);   //Call the 'trial' stepper

           // Compare with a high-quality stepper -- expect it works well for this step size (check!)
           //   This builds on its previous step to create the solution !
#ifdef COMPARE_TO_G4        
           g4ExactStepper->Stepper(yInX,dydxRef,stepLengthRef,youtX,yerrX); //call the reference stepper
#else
           exactStepperGV->StepWithErrorEstimate(yInX,dydxRef,stepLengthRef,youtX,yerrX); //call the reference stepper
#endif
        }
        //-> Then print the data
        cout.setf (ios_base::fixed);

        cout.precision(1);
        cout << setw(8) << total_step << " ";
        
        // Report Position        
        cout.precision(4);
        for(int i=0; i<3;i++)
            if(columns[i]){
               if( printSep ) cout << " | " ;  // Separator
               if( printInp ) cout << setw(wd-2)<< yIn[i] / mmGVf;
               if( printInpX ) cout << setw(wd-2)<< yInX[i] / mmRef;
               cout<<setw(wd)<< yout[i] / mmGVf ;
               // cout.precision(3);
               if( printRef ) cout<<setw(wd)<< youtX[i] / mmRef; // Reference Solution
               if( printErr ) cout<<setw(wd)<< errAmplif * (yerrX[i] / mmGVf) ;
               if( printDiff ){
                  cout.precision(epsDigits);
                  cout<<setw(wdf)<< yDiff[i];
               }               
               if( printOk || printRatio ){
                  yDiff[i] = yout[i] / mmGVf - youtX[i] / mmRef ;
                  yAver[i] = yout[i] / mmGVf + youtX[i] / mmRef ;
                  bool good =  fabs(yDiff[i]) <  0.5 * epsTol * fabs(yAver[i]);
                  if(printOk) cout<< " " << setw(3)<< (good ? "Ok " : "Bad" );
                  double ratio = (yAver[i] != 0.0) ? yDiff[i] / yAver[i] : 0.0;
                  if( printRatio ) cout << " " << setw(wd) << ratio/epsTol << " ";
               }
            }
        cout.unsetf (ios_base::fixed);

        // Report momentum 
        cout.setf (ios_base::scientific);
        for(int i=3; i<6;i++)
            if(columns[i]){
               if( printSep ) cout << " | " ;  // Separator
               if( printInp ) cout << setw(wd-1)<< yIn[i] / ppGVf << " ";
               if( printInpX ) cout << setw(wd-1)<< yInX[i] / ppRef << " ";
               cout<<setw(wd) << yout[i] / ppGVf ;
               if( printRef ) cout<<setw(wd)<< youtX[i] / ppRef; // Reference Solution
               if( printErr ) cout<< " " << setw(wd)<< errAmplif * yerrX[i] / ppGVf ;
               if( printDiff ){
                  cout.precision(epsDigits);
                  cout<<setw(wdf)<< yDiff[i];
               }                 
               if( printOk || printRatio ){
                  yDiff[i] = yout[i] / ppGVf - youtX[i] / ppRef ;
                  yAver[i] = yout[i] / ppGVf + youtX[i] / ppRef ;
                  bool good =  fabs(yDiff[i]) <  0.5 * epsTol * fabs(yAver[i]);
                  if( printOk ) cout<< " " << setw(3)<< (good ? "Ok " : "Bad" );
                  double ratio = (yAver[i] != 0.0) ? yDiff[i] / yAver[i] : 0.0;
                  if( printRatio ){
                     cout.unsetf (ios_base::scientific);
                     cout.setf (ios_base::fixed);                     
                     cout << " " << setw(wd) << ratio/epsTol << " ";
                     cout.unsetf (ios_base::fixed);
                     cout.setf (ios_base::scientific);                     
                  }
               }
            }
        cout.unsetf (ios_base::scientific);
        
        for(int i=0; i<6;i++)   // Print auxiliary components
        {
           double unitGVf=1;  
           double unitRef=1;
           // if( i < 3 )             // length / length

           if( i >= 3 ){
              unitGVf = ppGVf / mmGVf; //  dp / ds
              unitRef = ppRef / mmRef; //  dp / ds
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
               // cout << " dy/dx["<<i<<"] = ";
               if( printSep ) cout << " | " ;  // Separator
               // cout<<setw(wd)<< dydx[i] / unitGVf ; -- Not available from Driver
               // if( printRef ) 
               cout<<setw(wd)<< dydxRef[i] / unitRef; // Reference Solution
               cout.precision(epsDigits);               
               if( 0 ) {  //  ( printDiff ){  -- Not meaningful
                  cout<<setw(wdf)<< ( dydx[i] / unitGVf )
                                   - ( dydxRef[i] / unitRef );
               }
           }
           // if( i == 2 )     { cout.unsetf(ios_base::fixed);      }
           // else if ( i==5 ) { cout.unsetf(ios_base::scientific); }
        }
        cout.unsetf(ios_base::scientific);
        if( j > 0 )  // Step 0 did not move -- printed the starting values
        {
           // cout.unsetf(ios_base::scientific);
           cout.setf (ios_base::fixed);                         
           cout.precision(2);
           cout<<setw(wd) << atan2(yout[1],yout[0])/degree;
                             // atan(yout[1]/yout[0])/degree;
           if( printRef ) 
             cout<<setw(wd) << atan2(youtX[1],youtX[0])/degree;
           
           //Copy yout into yIn
           for(int i=0; i<6;i++){
              // yIn[i] = yout[i];
              yInX[i] = youtX[i];
           }
        }
        
        cout<<"\n";
    }
    
    /*-----------------END-STEPPING------------------*/

    /*------ Clean up ------*/
    gvEquation->InformDone();    

#ifndef COMPARE_TO_G4
    gvEquation2->InformDone();    
#endif
    delete myStepper;
    delete exactStepper;
    delete gvEquation;
    delete gvEquation2;    
    delete gvField;
    
    cout<<"\n\n#-------------End of output-----------------\n";
    
}
