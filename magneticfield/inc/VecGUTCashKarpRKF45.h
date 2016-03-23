//
// Runge-Kutta Stepper using Cash Karp's RK tableau
//
// Adapted from 'GUTCashKarpRKF45' by Qieshen Xie, GSoC 2014
//         (derived from G4CashKarpRKF45)
//
// First version:  John Apostolakis,  4 Nov 2015
//
#ifndef VECGUTCASHKARPRKF45_H
#define VECGUTCASHKARPRKF45_H


// #include "TMagErrorStepper.h" //for sake of GUIntegrationNms::NumVars
// or add as below
namespace GUIntegrationNms
{
   constexpr unsigned int NumVarBase  = 8;  //
}

#include "AlignedBase.h"
#include "VecTMagFieldEquation.h"

#define INLINERHS 1
#define RemoveAuxStepper

#ifdef INLINERHS
#define REALLY_INLINE   inline __attribute__((always_inline))
#else
#define REALLY_INLINE   inline
#endif


template
<class T_Equation, unsigned int Nvar>
class VecGUTCashKarpRKF45 : public AlignedBase
{
  public:
    static constexpr unsigned int sOrderMethod= 4;
    static constexpr unsigned int sNstore = (GUIntegrationNms::NumVarBase > Nvar) ? GUIntegrationNms::NumVarBase : Nvar;
    // std::max( GUIntegrationNms::NumVarBase,  Nvar);
    // static const IntegratorCorrection = 1./((1<<4)-1);
    inline double IntegratorCorrection() { return 1./((1<<4)-1); }

    inline
    VecGUTCashKarpRKF45( T_Equation *EqRhs,
                              unsigned int numStateVariables=0,
                              bool primary=true);

    VecGUTCashKarpRKF45( const VecGUTCashKarpRKF45& );
    
    virtual ~VecGUTCashKarpRKF45();

    VecGUTCashKarpRKF45<T_Equation,Nvar>* Clone() const;

    template <class Backend>
    REALLY_INLINE
    void StepWithErrorEstimate(const typename Backend::precision_v* yInput,    // Consider __restrict__
                               const typename Backend::precision_v*  dydx,
                                     typename Backend::precision_v   Step,
                                     typename Backend::precision_v*  yOut,
                                     typename Backend::precision_v*  yErr) ;

    template <class Backend>
    REALLY_INLINE
    void RightHandSideInl(typename Backend::precision_v y[], typename Backend::precision_v dydx[]) 
    {fEquation_Rhs->T_Equation::template RightHandSide<Backend>(y, dydx);}

    template <class Backend>
    inline 
    void RightHandSideVIS( const typename Backend::precision_v y[], 
                                 typename Backend::precision_v charge, 
                                 typename Backend::precision_v dydx[] ); 

    void SetEquationOfMotion(T_Equation* equation);

    inline T_Equation *GetEquationOfMotion() { return  fEquation_Rhs; }
    inline const T_Equation *GetEquationOfMotion() const { return  fEquation_Rhs; }  
    
    void InformDone() { GetEquationOfMotion()->InformDone();}
    
    inline unsigned int  GetNumberOfVariables() const;
    // Get the number of variables that the stepper will integrate over.

    inline unsigned int  GetNumberOfStateVariables() const;
    // Get the number of variables of state variables (>= above, integration)
    
    unsigned int IntegratorOrder() const { return fIntegrationOrder; };

  private:

  
    VecGUTCashKarpRKF45& operator=(const VecGUTCashKarpRKF45&) = delete;
    // private assignment operator.

  private:
    // 'Invariant' during integration - the pointers must not change
    // -----------
    T_Equation* fEquation_Rhs; // corresponds to fAbstrEquation of GUVIntegrationStepper
    bool        fOwnTheEquation;
    #ifndef RemoveAuxStepper
    VecGUTCashKarpRKF45* fAuxStepper;
    #endif

/*    // State -- intermediate values used during RK step
    // -----        
    typename Backend::precision_v ak2[sNstore];
    typename Backend::precision_v ak3[sNstore];
    typename Backend::precision_v ak4[sNstore];
    typename Backend::precision_v ak5[sNstore];
    typename Backend::precision_v ak6[sNstore];
    typename Backend::precision_v ak7[sNstore];
    typename Backend::precision_v yTemp[sNstore];
    typename Backend::precision_v yIn[sNstore];
    // scratch space

    // State -- values used for subsequent call to DistChord
    // -----
    
    typename Backend::precision_v  fLastStepLength;

    typename Backend::precision_v* fLastInitialVector;
    typename Backend::precision_v* fLastFinalVector;
    typename Backend::precision_v* fLastDyDx;
*/
    const unsigned int fIntegrationOrder; // RK or similar order - if any. Else 0
    const unsigned int fNoIntegrationVariables; // # of Variables in integration
    const unsigned int fNoStateVariables;       // # required for FieldTrack

};

template <class T_Equation, unsigned int Nvar>
inline
unsigned int VecGUTCashKarpRKF45<T_Equation,Nvar>
  ::GetNumberOfVariables() const
{
  return fNoIntegrationVariables;
}

template <class T_Equation, unsigned int Nvar>
inline
unsigned int VecGUTCashKarpRKF45<T_Equation,Nvar>
  ::GetNumberOfStateVariables() const
{
  return fNoStateVariables;
}

template <class T_Equation, unsigned int Nvar>
inline
VecGUTCashKarpRKF45<T_Equation,Nvar>::
  VecGUTCashKarpRKF45( T_Equation *EqRhs,
                       unsigned int numStateVariables,
                       bool primary                  )
  : fEquation_Rhs(EqRhs),
    fOwnTheEquation(primary),
   #ifndef RemoveAuxStepper
    fAuxStepper(0),
   #endif
    // fLastStepLength(0.),
    fIntegrationOrder(sOrderMethod),
    fNoIntegrationVariables(Nvar),
    fNoStateVariables((numStateVariables>0) ? numStateVariables : sNstore)
{
 #if 0
  std::cout<<"\n----Entered constructor of VecGUTCashKarpRKF45 "<<std::endl;
  std::cout<<"----In VecGUTCashKarpRKF45 constructor, Nvar is: "<<Nvar<<std::endl;
 #endif
  assert( dynamic_cast<VecTMagFieldEquation*>(EqRhs) != 0 );

  assert( (numStateVariables == 0) || (numStateVariables >= Nvar) );
 
/*  typedef typename Backend::precision_v Double_v;

  fLastInitialVector = new Double_v[sNstore] ;
  fLastFinalVector   = new Double_v[sNstore] ;
  fLastDyDx          = new Double_v[sNstore];*/

 #ifndef RemoveAuxStepper
  if( primary )
  {
    // Reuse the Equation of motion in the Auxiliary Stepper      
    VecGUTCashKarpRKF45<T_Equation,Nvar> fAuxStepper2(EqRhs, numStateVariables, false);
    fAuxStepper = &fAuxStepper2;
    // fAuxStepper = new VecGUTCashKarpRKF45<T_Equation,Nvar>(EqRhs, numStateVariables, false);
    #ifdef DEBUGAnanya
    std::cout<<"----VecGUTCashKarpRKF45 Auxiliary Stepper being made"<<std::endl;
    #endif 
  }
 #endif 

  #if 0
   std::cout<<"----end of constructor of VecGUTCashKarpRKF45"<<std::endl;
  #endif
}


template <class T_Equation, unsigned int Nvar>
void VecGUTCashKarpRKF45<T_Equation,Nvar>::
  SetEquationOfMotion(T_Equation* equation)
{
  if (equation != 0)
  {
    fEquation_Rhs= equation;
  }
}

//  Copy - Constructor
// 
template <class T_Equation,unsigned int Nvar>
inline
VecGUTCashKarpRKF45<T_Equation,Nvar>::
  VecGUTCashKarpRKF45( const VecGUTCashKarpRKF45& right )
  : fEquation_Rhs( (T_Equation*) 0 ),
    fOwnTheEquation(true),

    #ifndef RemoveAuxStepper
    fAuxStepper(0),   //  May overwrite below
    #endif 

   fIntegrationOrder(sOrderMethod),
   fNoIntegrationVariables(Nvar),
   fNoStateVariables(right.GetNumberOfStateVariables())

    // fLastStepLength(0.)
    // fPrimary( right.fPrimary )
{
/*  #ifdef DEBUGAnanya
   std::cout<<"----Entered constructor of VecGUTCashKarpRKF45 "<<std::endl;
 #endif*/
  // if( primary )
  SetEquationOfMotion( new T_Equation( *(right.fEquation_Rhs)) );
  fOwnTheEquation=true;
   // fEquation_Rhs= right.GetEquationOfMotion()->Clone());
  
  assert( dynamic_cast<TemplateGUVEquationOfMotion<Backend>*>(fEquation_Rhs) != 0 );  
  assert( this->GetNumberOfStateVariables() >= Nvar);

#if 1
  // if( verbose )
  std::cout << " VecGUTCashKarpRKF45 - copy constructor: " << std::endl
            << " Nvar = " << Nvar << " Nstore= " << sNstore 
            << " Own-the-Equation = " << fOwnTheEquation << std::endl;
#endif   

#ifndef RemoveAuxStepper
  if( right.fAuxStepper )
  {
    // Reuse the Equation of motion in the Auxiliary Stepper
    fAuxStepper = new VecGUTCashKarpRKF45(fEquation_Rhs, this->GetNumberOfStateVariables(), false);
    std::cout<<"Auxiliary stepper made"<<std::endl;
  }
#endif
}



template <class T_Equation,unsigned int Nvar>
REALLY_INLINE
VecGUTCashKarpRKF45<T_Equation,Nvar>::~VecGUTCashKarpRKF45()
{
#if 0
  std::cout<<"----- CashKarp destructor0"<<std::endl;
#endif

/*   delete[] fLastInitialVector;
   delete[] fLastFinalVector;
   delete[] fLastDyDx;*/

  #ifndef RemoveAuxStepper
  std::cout<<fAuxStepper<<std::endl;
  delete fAuxStepper;
  #endif 


  if( fOwnTheEquation )
    delete fEquation_Rhs; // Expect to own the equation, except if auxiliary (then sharing the equation)

#if 0
  std::cout<<"----- CashKarp destructor1"<<std::endl;
#endif

}


template <class T_Equation, unsigned int Nvar>
VecGUTCashKarpRKF45<T_Equation,Nvar>* 
VecGUTCashKarpRKF45<T_Equation,Nvar>::Clone() const
{
  // return new VecGUTCashKarpRKF45( *this );
  return new VecGUTCashKarpRKF45<T_Equation,Nvar>( *this );   
}


template <class T_Equation, unsigned int Nvar>
template <class Backend> 
inline void 
VecGUTCashKarpRKF45<T_Equation,Nvar>::
  RightHandSideVIS( const  typename Backend::precision_v y[], 
                           typename Backend::precision_v charge,
                           typename Backend::precision_v dydx[] )
{
   fEquation_Rhs->T_Equation::template RightHandSide<Backend>(y, charge, dydx);
}

template <class T_Equation, unsigned int Nvar>
template <class Backend>
inline void
VecGUTCashKarpRKF45<T_Equation,Nvar>::
  StepWithErrorEstimate(const typename Backend::precision_v  yInput[],    
                        const typename Backend::precision_v  dydx[],
                              typename Backend::precision_v   Step,
                              typename Backend::precision_v  yOut[],
                              typename Backend::precision_v  yErr[])
{
  typedef typename Backend::precision_v Double_v;

  Double_v ak2[sNstore];
  Double_v ak3[sNstore];
  Double_v ak4[sNstore];
  Double_v ak5[sNstore];
  Double_v ak6[sNstore];
  // Double_v ak7[sNstore];
  Double_v yTemp[sNstore];
  Double_v yIn[sNstore];

  Double_v  fLastInitialVector[sNstore];
  Double_v  fLastFinalVector[sNstore];
  Double_v  fLastDyDx[sNstore];

  Double_v  fLastStepLength;

  // const int nvar = 6 ;
  // const double a2 = 0.2 , a3 = 0.3 , a4 = 0.6 , a5 = 1.0 , a6 = 0.875;
  unsigned int i;

  const double  b21 = 0.2 ,
        b31 = 3.0/40.0 , b32 = 9.0/40.0 ,
        b41 = 0.3 , b42 = -0.9 , b43 = 1.2 ,

        b51 = -11.0/54.0 , b52 = 2.5 , b53 = -70.0/27.0 ,
        b54 = 35.0/27.0 ,

        b61 = 1631.0/55296.0 , b62 =   175.0/512.0 ,
        b63 =  575.0/13824.0 , b64 = 44275.0/110592.0 ,
        b65 =  253.0/4096.0 ,

        c1 = 37.0/378.0 , c3 = 250.0/621.0 , c4 = 125.0/594.0 ,
        c6 = 512.0/1771.0 ,
        dc5 = -277.0/14336.0 ;

  const double dc1 = c1 - 2825.0/27648.0 ,  
        dc3 = c3 - 18575.0/48384.0 ,
        dc4 = c4 - 13525.0/55296.0 , 
        dc6 = c6 - 0.25 ;

  // Initialise time to t0, needed when it is not updated by the integration.
  //       [ Note: Only for time dependent fields (usually electric) 
  //                 is it neccessary to integrate the time.] 
  //yOut[7] = yTemp[7]   = yIn[7]; 

  //  Saving yInput because yInput and yOut can be aliases for same array
  for(i=0;i<Nvar;i++) 
  {
    yIn[i]=yInput[i];
  }
  // RightHandSideInl(yIn, dydx) ;              // 1st Step

  for(i=0;i<Nvar;i++) 
  {
    yTemp[i] = yIn[i] + b21*Step*dydx[i] ;
  }
  this->template RightHandSideInl<Backend>(yTemp, ak2) ;              // 2nd Step

  for(i=0;i<Nvar;i++)
  {
    yTemp[i] = yIn[i] + Step*(b31*dydx[i] + b32*ak2[i]) ;
  }
  this->template RightHandSideInl<Backend>(yTemp, ak3) ;              // 3rd Step

  for(i=0;i<Nvar;i++)
  {
    yTemp[i] = yIn[i] + Step*(b41*dydx[i] + b42*ak2[i] + b43*ak3[i]) ;
  }
  this->template RightHandSideInl<Backend>(yTemp, ak4) ;              // 4th Step

  for(i=0;i<Nvar;i++)
  {
    yTemp[i] = yIn[i] + Step*(b51*dydx[i] + b52*ak2[i] + b53*ak3[i] +
            b54*ak4[i]) ;
  }
  this->template RightHandSideInl<Backend>(yTemp, ak5) ;              // 5th Step

  for(i=0;i<Nvar;i++)
  {
    yTemp[i] = yIn[i] + Step*(b61*dydx[i] + b62*ak2[i] + b63*ak3[i] +
            b64*ak4[i] + b65*ak5[i]) ;
  }
  this->template RightHandSideInl<Backend>(yTemp, ak6) ;              // 6th Step

  for(i=0;i<Nvar;i++)
  {
    // Accumulate increments with proper weights

    yOut[i] = yIn[i] + Step*(c1*dydx[i] + c3*ak3[i] + c4*ak4[i] + c6*ak6[i]) ;
  }
  for(i=0;i<Nvar;i++)
  {
    // Estimate error as difference between 4th and
    // 5th order methods

    yErr[i] = Step*(dc1*dydx[i] + dc3*ak3[i] + dc4*ak4[i] +
            dc5*ak5[i] + dc6*ak6[i]) ;
  #ifdef DEBUGAnanya
    // std::cout<< "----In Stepper, yerrr is: "<<yErr[i]<<std::endl;
  #endif 
  }
  for(i=0;i<Nvar;i++)
  {
    // Store Input and Final values, for possible use in calculating chord
    fLastInitialVector[i] = yIn[i] ;
    fLastFinalVector[i]   = yOut[i];
    fLastDyDx[i]          = dydx[i];
  }
  fLastStepLength =Step;

  return ;
}



#endif /*VecGUTCashKARP_RKF45 */
