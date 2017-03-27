
#ifndef GammaModelWrapper_H
#define GammaModelWrapper_H

#include <string>
#include <vector>
#include <cassert>

#include "base/VecPhys.h"
#include "base/SystemOfUnits.h"   //  The OLD units - CLHEP namespace

// #include "Particle.h"
#include "Electron.h"
#include "Positron.h"
#include "Gamma.h"

// #include "VecCore/Backend/Scalar.h"
#include "GeantTaskData.h"

#include "PhysicsData.h"  // Type of Physics data in 'task data'

#include "PhysicalConstants.h"

#include "EMModel.h"
// #include "EmModelBase.h"   // VecPhys
#include "LightTrack.h"
#include "Element.h"

#include "GUTrack.h"

#include "Material.h"
#include "MaterialCuts.h"

namespace geantphysics {

// forward declarations
// class GeantTaskData;   
// class Material;
// class MaterialCuts;
class MaterialProperties;
class Element;
class Particle;
class PhysicsParameters;
class EMElementSelector;

/**
 * @brief   Wrapper class for gamma physics models.
 * @class   GammaModelWrapper
 * @author  J. Apostolakis S.Y. Jun
 * @date    March 2017
 */

template <class PhysModelType, bool isConversion= false>
class GammaModelWrapper : public EMModel 
{
  public:
  /**
   * @brief Constructor
  GammaModelWrapper(const std::string&  name,
     @param[in]  vecPhysModel,
     @param[in]  emittedType    int:  code of the particle emitted by the interaction
     @param[in]  gammaSurvives  boolean: does the projectile survive ?
   */
  GammaModelWrapper(const std::string&  name,
                    PhysModelType*      vecPhysModel,
                    int                 emittedType,
                    bool                gammaSurvives  // Static information -- for checking
     );  

  /**
   * @brief DTR
   */
  ~GammaModelWrapper() {}

//
// The following 5 virtual methods might be implemented by the derived electromagnetic models.
//
  // will set the physics parameter object that belongs to the regions(s) where this model is active
  /**
   * @brief The base class implementation Initialize() must be called in the first
   *        line because to set the PhysicsParameters member variable for the 
   *        regions(s) where this model is active.
   *        This method is called when its EMPhysicsProcess, is initialized by the
   *        EmModelManager member of the EMPhysicsProcess.
   */
  void  Initialize() override final;

  /**
   * @brief Method to compute stopping power relevant to energy loss models only.
   *
   */
  double ComputeDEDX(const MaterialCuts * /*matcut*/, double /*kinenergy*/,
                     const Particle * /*particle*/,   bool   /*istotal=false*/
                    ) override final { return 1.0; }

  /**
   * @brief Method to compute macroscopic cross section for a given MaterialCuts, 
   *        Partcile, kinetic energy.
   *
   * Mandatory method called at initialization to build the lambda tables through 
   * the corresponding PhysicsProcess by the PhysicsManagerPerParticle object.
   *
   * @param[in] matcut      Pointer to the MaterialCuts object
   * @param[in] kinenergy   Kinetic energy of the Particle
   * @param[in] particle    Pointer to the Particle object
   * @return    Macroscopic cross section computed by the given electromagnetic model in internal [1/length] units for
   *            the given Particle, MaterialCuts/Material and particle kinetic energy combination.
   */
  double ComputeMacroscopicXSection(const MaterialCuts * /*matcut*/, double /*kinenergy*/,
                                    const Particle * /*particle*/) override final;

  /**
   * @brief Compute atomic cross section for given Element, MaterialCuts, Particle, kinetic energy.
   *
   * Method is called at initialization to build the lambda tables through the corresponding 
   *  PhysicsProcess by
   *  - the PhysicsManagerPerParticle object, and 
   *  - the EMElementSelector object, if this was requested in the derived class
   *     Initialize() method.
   * 
   * Note, that the MaterialCut object, that corresponds to the material-region to
   *    which the current element belongs to is also provided to include any 
   *    material, or cut dependences. However, this infomation does not necessary
   *    to use by each GammaModelWrapper-s (e.g. those models that 
   *    describes discrete interaction that doesn't depend on production
   *    threshold and material dependent properties won't use this extra information).
   *
   * @param[in] elem        Pointer to the Element object relating request
   * @param[in] matcut      Pointer to the MaterialCuts object 
   * @param[in] kinenergy   Kinetic energy of the Partcile
   * @param[in] particle    Pointer to the Partcile object
   * @return    Atomic cross section computed by the given electromagnetic model in
   *            internal [length^2] units for the given ELement, Particle, 
   *            MaterialCuts/Material and particle kinetic energy combination.
   */
  double ComputeXSectionPerAtom(const Element * /*elem*/,
                                const MaterialCuts * /*matcut*/,
                                double /*kinenergy*/,
                                const Particle * /*particle*/) override final;

  int    SampleSecondaries(LightTrack & /*track*/, std::vector<LightTrack> & /*sectracks*/,
                                   Geant::GeantTaskData * /*td*/) override final;

  /**
   * @brief Method to obtain minimum primary particle kinetic energy at which the discrete part
   *         (if any) of the interaction can happen i.e. kinematic minimum for the discrete 
   *         part of the interaction.
   *
   * All discrete models which have production-cut dependence *must* implement this method. 
   *    It is used e.g. to build the energy grid of the target element selector (if it was
   *    by the derived model) for providing run-time functionality to select target element
   *    for the discrete interaction.
   *
   * @param[in] matcut      Pointer to the MaterialCuts object (for interaction)
   * @param[in] particle    Pointer to the Particle object 
   * @return    Minimum primary particle kinetic energy (at which the discrete interaction
   *            can happen) in internal [energy] units (for the given Particle, 
   *            MaterialCuts/Material combination).
   */
  double MinimumPrimaryEnergy(const MaterialCuts * /*matcut*/, const Particle * /*part*/) const final
  { return  ( isConversion ? 2.0 * geant::kElectronMassC2 : 0.0 ); }

protected:
  // Rotate outgoing particle's momentum  'finalMomentum' using projectile's initial vector
  void RotateToLabFrame(       vecgeom::Vector3D<double> &outgoingVector,
                         const vecgeom::Vector3D<double> &initialVec ) const;

  // Some short cuts will be used in methods below, for speed - knowing 
  //    the type of primary and the implementation of the interactions
  GUTrack    ConvertLightToGUTrack(LightTrack &lightTrack) const;
  LightTrack ConvertGUtoLightTrack( GUTrack      &guTrack,
                                    LTrackStatus  status,  
                                    int           materialCutCoupleIndex,
                                    double        time     = -1.0,
                                    double        weight   =  1.0,                                   
                                    double        eDeposit =  0.0,
                                    ExtraInfo*    aExtraInfo = 0
                                   ) const ;

private:
  PhysModelType*  fVecPhysModel;  
  // double          fModelminimumPrimaryEnergy= 0.0;
  const int       fEmittedType;    // Particle type of emitted particle (GeantV id)
  const bool      fGammaSurvives;  // Is projectile gamma expected to survive interaction ? 
};

// *********************  Implementation *******************************************
// *********************************************************************************

template <class PhysModelType, bool isConversion>
   GammaModelWrapper<PhysModelType, isConversion>::
   GammaModelWrapper(const std::string&  name,
                    PhysModelType*      vecPhysModel,
                     int                 emittedType,
                     bool                gammaSurvives                     
      )
     : EMModel( name ),
       fVecPhysModel(vecPhysModel),
       fEmittedType(emittedType),
       fGammaSurvives(gammaSurvives)
{
   std::cout << "Constructor called for Gamma Model Wrapper with process name= "
             << name << "  physics model ptr= " << vecPhysModel
             << std::endl;
}

// **************************************************************************

template <class PhysModelType, bool isConversion>         
void
   GammaModelWrapper<PhysModelType,isConversion>::Initialize()
{
  std::cout << "GammaModelWrapper method Initialize() called for model "
            << this->GetName() << std::endl;

  //  fVecPhysModel->template Initialize<vecCore::backend::Scalar>Initialization();
  fVecPhysModel->Initialization();
}

// **************************************************************************

template <class PhysModelType, bool isConversion>         
GUTrack
   GammaModelWrapper<PhysModelType,isConversion>::
      ConvertLightToGUTrack(LightTrack &lightTrack) const
{
  GUTrack guTrack;
  const int gammaCode= Gamma::Definition()->GetInternalCode();

  if (lightTrack.GetGVcode() != gammaCode ) {
     std::cerr << " ERROR: Expected light track with gamma (code= " << gammaCode << " ) "
               << " .  Track has GV code = " << lightTrack.GetGVcode() << std::endl;
  }
  // assert(lightTrack.GetGVcode() == gammaCode );
  
  guTrack.status= 1;                            // lightTrack.GetTrackStatus();
  guTrack.particleType= lightTrack.GetGVcode();
  guTrack.id= lightTrack.GetTrackIndex();       // counter
  guTrack.parentId= 0;                          // id of parent
  guTrack.proc= lightTrack.GetProcessIndex();   // index of physics process
  // guTrack.x= 0.0;
  // guTrack.y= 0.0;
  // guTrack.z= 0.0;

  double eKin= lightTrack.GetKinE();
  // double momentum= std::sqrt(eKin*(eKin+2*restMass));
  double momentum= eKin;  // Shortcut for gamma :  E = E_kin = p * c
  guTrack.px= momentum * lightTrack.GetDirX();
  guTrack.py= momentum * lightTrack.GetDirY();
  guTrack.pz= momentum * lightTrack.GetDirZ();
  guTrack.E= eKin;
  guTrack.q= 0;  // Shortcut for photon 
  guTrack.nint= lightTrack.GetNumOfInteractionLegthLeft();
  //  guTrack.lambda= lightTrack.GetIntLen();   // interaction length
  guTrack.s= lightTrack.GetStepLength();    // step length ??

  assert(guTrack.particleType == gammaCode);
  
  return guTrack;
}

// **************************************************************************

inline
double
VectorToUnitAndMag( const vecgeom::Vector3D<double> & vectorIn,
                          vecgeom::Vector3D<double> & outputDir ) // const
{
  double magnitude = vectorIn.Mag();
  double invFinalMom = 1.0 / magnitude;
  outputDir = invFinalMom * vectorIn;

  return magnitude;
}

// **************************************************************************
inline bool EstimatedMass( double Energy, vecgeom::Vector3D<double> momentum )
{
   double xAbs = fabs(momentum.x());
   double yAbs = fabs(momentum.y());
   double zAbs = fabs(momentum.z());
   double maxXY= std::max( xAbs, yAbs );
   double minXY= std::min( xAbs, yAbs );

   double maxAll = std::max( maxXY, zAbs );
   double minUP  = std::min( maxXY, zAbs );
   double middle, minAll;
   if(  minUP > minXY ) { middle = minUP; minAll= minXY; }
   else                 { middle = minXY; minAll= minUP; }   
   // if(  zAbs > minXY ) { middle = zAbs;  minAll= minXY; }
   // else                { middle = minXY; minAll= zAbs; }
   assert ( minAll <= middle );
   assert ( middle <= maxAll );

#if 1   
   double diff = ( Energy - maxAll) * ( Energy + maxAll ) - middle*middle - minAll*minAll;
#else
   double diff1 = ( Energy - maxAll) * ( Energy + maxAll );
   double s1 = diff1 > 0.0 ?  sqrt(diff1) : 0.0; 
   double diff2 = ( s1 - middle ) * ( s1 + middle ); 
   double s2 = diff2 > 0.0 ?  sqrt(diff1) : 0.0;
   double diff = ( s2 - minAll ) * ( s2 + minAll );
#endif   
   diff = ( diff < 0.0 ) ? 0.0 : diff ;
   return ( sqrt(diff) );
}

inline bool GoodMass( double Energy, vecgeom::Vector3D<double> momentum, double mass )
{
   const double epsilonOK = 1.0e-5;
   
   double diff = ( Energy - mass ) * ( Energy - mass ) - momentum.Mag2();
   return ( diff  < epsilonOK * mass * mass );
}

// **************************************************************************

template <class PhysModelType, bool isConversion>
void
GammaModelWrapper<PhysModelType,isConversion>::
   RotateToLabFrame(       vecgeom::Vector3D<double> & outputVec,
                     const vecgeom::Vector3D<double> & initialDir ) const
{
  // using  ThreeVector= vecgeom::Vector3D<double>;   
  // double outDirX= invFinalMom * outputVec.x(), outDirY = invFinalMom * outputVec.y(), outDirZ= invFinalMom * outputVec.z();
  double outDirX= outputVec.x(), outDirY= outputVec.y(), outDirZ= outputVec.z();
  
  // double mag= initialDir.Mag();
  assert (  std::fabs(initialDir.Mag() - 1.0 ) < 1.0e-6 );
  // double invMag= mag > 0.0 ? 1.0 / mag : 1.0; 
  // ThreeVector dirInit= invMag * initialVec;

  EMModel::RotateToLabFrame(outDirX, outDirY, outDirZ,
                            -initialDir.x(), -initialDir.y(), -initialDir.z());                            
  //                         initialDir.x(), initialDir.y(), initialDir.z());

  outputVec.Set(outDirX, outDirY, outDirZ );  
  return;
}
   
// **************************************************************************

const double electronMass= 1000.0 * geant::kElectronMassC2;

// static atomic<int> numBadEPMass= 0;

inline
bool CheckEnergyMomentum( const vecgeom::Vector3D<double>& Momentum,
                          double             energy,
                          int                particleCode,
                          const std::string& name,
                          bool               verbose= false)
{
  const  int electronCode = Electron::Definition()->GetInternalCode();
  const  int positronCode = Positron::Definition()->GetInternalCode();  
  const  int    gammaCode =    Gamma::Definition()->GetInternalCode();
  bool   good= false;

  if( particleCode == gammaCode )
  {
     good= GoodMass( energy, Momentum, 0.0 );
     if( ! good ) { 
        if( verbose ) { // || (++numBadEPMass % 100 == 0) ){
           double pEratio = Momentum.Mag2() / (energy*energy);
           std::cerr << " ERROR> Incorrect p/E ratio for gamma, r = " << pEratio
                     << " r-1= "  << pEratio - 1.0 << std::endl;
             // << "  occurance # " << numBadEPMass << std::endl;
        }
     }
     assert( fabs( Momentum.Mag() / energy ) < 1.0e-4 * energy );
  }
  else if( particleCode == electronCode || particleCode == positronCode )
  {
     good = GoodMass( energy, Momentum, electronMass );
     if( verbose && ! good ) {
        double estimM= EstimatedMass( energy, Momentum );
        std::cerr << " ERROR> BAD mass from E/p - for secondary : "
                  << " Mass estimated = " << estimM
                  << " vs expected = " << electronMass
                  << "  rel. diff = " << (estimM/electronMass) - 1.0
                  << " momentum= " << Momentum << "  E= " << energy
                  << " label= " << name
                  << std::endl;
     }
     assert( GoodMass( energy, Momentum, electronMass ) );
  }
  return good;
}

// **************************************************************************

inline
bool CheckEPGamma( const GUTrack& aTrack, const std::string &label, bool verbose= false )
{
  using  ThreeVector= vecgeom::Vector3D<double>;     
  const  int  gammaCode =    Gamma::Definition()->GetInternalCode();
  const  ThreeVector  momentum(aTrack.px, aTrack.py, aTrack.pz);
  
  return CheckEnergyMomentum( momentum, aTrack.E, gammaCode, label, verbose);
}

// **************************************************************************

inline
bool CheckEPMass( const GUTrack& guTrack, const std::string &label, bool verbose )  
{
  vecgeom::Vector3D<double> momentum(guTrack.px, guTrack.py, guTrack.pz);
  return CheckEnergyMomentum( momentum, guTrack.E, guTrack.particleType,
                              label, verbose );
}
   
// **************************************************************************
inline
bool CheckConservationLaws( const typename vecgeom::Vector3D<double> initialMom,   double initialEn,
                            const typename vecgeom::Vector3D<double> finalProjMom, double finalE,
                            const typename vecgeom::Vector3D<double> secondaryMom, double secondaryE,
                            bool verbose = false )
{
  using  ThreeVector= vecgeom::Vector3D<double>;   
  double initialMomentumMag= initialMom.Mag();
  const double epsilon = 1.0e-5;
     
  ThreeVector totalMomentumOut= finalProjMom + secondaryMom;  
  ThreeVector diffMomentum=     initialMom - totalMomentumOut;
  double      diffEnergy  =     initialEn - guProjectile.E - guSecondary.E;

  bool bigDiffP = ( totalMomentumOut.Mag() > epsilon * initialMomentumMag );
  bool bigDiffE = ( std::fabs(diffEnergy) > epsilon * initialEn );

  bool bigDiff = bigDiffE || bigDiffP ;  
  
  if( bigDiff || verbose )
  {
     int oldPrc= cout.precision(10);     
     cout << " particle 1 - out " << finalProjMom << endl;
     cout << " particle 2 - out " << secondaryMom << endl;
     cout << " Sum P out        " << totalMomentumOut << endl;
     cout << " P     in         " << initialMom   << endl;     
     cout.precision(oldPrc);
     
     if( initialMomentumMag == 0.0 ) { initialMomentumMag= -0.1e-70; }
     cout << " DiffE = ";
     if (diffEnergy == 0.0) {cout << "0 ";} else { cout << diffEnergy; }
     cout << " (UnRotated) d-Momentum= " << diffMomentum
          << " mag(dP)= " << diffMomentum.Mag() << " relative: " <<  diffMomentum.Mag() / initialMomentumMag
          << " d|P|= " << initialMomentumMag - totalMomentumOut.Mag()
          << " relative diff " << totalMomentumOut.Mag() / initialMomentumMag - 1.0
          << endl;
  }   // DEBUG printout - END
}
   
// initialEn - guProjectile.E - guSecondary.E;

/*    
inline
bool CheckEPMass( const LightTrack& track, const std::string &label, bool verbose )  
{
  vecgeom::Vector3D<double> momentum(track.GetDirX(), track.GetDirY(), guTrack.GetDirZ());
  momentum *= 
  return CheckEnergyMomentum( momentum, guTrack.E, guTrack.particleType,
                       "Convert GU to Light Track", verbose );
}
*/
   
// **************************************************************************

   
template <class PhysModelType, bool isConversion>
LightTrack 
GammaModelWrapper<PhysModelType,isConversion>::ConvertGUtoLightTrack(
                                         GUTrack      &guTrack,
                                         LTrackStatus  status,
                                         int           materialCutCoupleIndex,
                                         double        time,
                                         double        weight,
                                         double        eDeposit,
                                         ExtraInfo*    aExtraInfo
                                         ) const
{
  const  int gammaCode = Gamma::Definition()->GetInternalCode();

  double aMass = (guTrack.id == gammaCode) ? 0.0: geant::kElectronMassC2;

  // Basic check of relation between E, p & mass
  vecgeom::Vector3D<double> momentum(guTrack.px, guTrack.py, guTrack.pz);
  CheckEnergyMomentum( momentum, guTrack.E, guTrack.particleType,
                       "Convert GU to Light Track", true ); // verbose
  
  double momentumMag = ( guTrack.particleType != gammaCode ) ? momentum.Mag() : guTrack.E;
  double invMag= 1.0 / momentumMag;
  momentum *= invMag;
  
  return LightTrack( 
             status,                // LTrackStatus::kAlive
             guTrack.particleType,  // aGVcode,
             guTrack.parentId,      // --> needed by PostStepDoIt of PhysicsProcessHandler ( used to get => x,y,z,..com )
             materialCutCoupleIndex,
             0, // aProcessIndex ==> this particle is new, no process is chosen for it
             0, // aTargetZ,
             0, // aTargetN,
             momentum.x(), momentum.y(), momentum.z(), // invMag * guTrack.px, invMag * guTrack.py, invMag * guTrack.pz,
             guTrack.E,
             aMass,
             time,
             weight,
             guTrack.s,
             eDeposit,
             guTrack.nint,    // aNintLen,
             guTrack.lambda,  // aIntLen,
             aExtraInfo
         );
}

// **************************************************************************
   
template <class PhysModelType, bool isConversion>
double
GammaModelWrapper<PhysModelType,isConversion>::
    ComputeXSectionPerAtom(const  Element      *elem,
                           const  MaterialCuts * /*matcut*/ ,
                           double               kinEnergy,
                           const Particle *     /*particle*/ )
{
   int  zElement =  (int) (elem-> GetZ()) ;
   double xsec;

   static bool firstCall= true;
   if( firstCall ) {
     std::cout << " GammaModelWrapper - method ComputeXSectionPerAtom() called for model " << this->GetName() << std::endl;
     firstCall= false;
   }

   std::cout << " GMWrapper CXSPA call - model " << this->GetName()
             << " Elem= " << elem->GetZ() << " Ekin= " << kinEnergy << std::endl;
   // Use the  CrossSectionKernel() method of
   //  template <class Backend>
   //  VECCORE_ATT_HOST_DEVICE typename Backend::Double_v
   //      CrossSectionKernel(typename Backend::Double_v          energyIn,
   //                         Index_v<typename Backend::Double_v> zElement);
   xsec= fVecPhysModel -> template CrossSectionKernel<vecCore::backend::Scalar> ( kinEnergy, zElement );

   // Use 
   //    VECCORE_ATT_HOST double GetG4CrossSection(int Z, double energyIn);
   // xsec= fVecPhysModel->GetG4CrossSection( zElement, kinEnergy);
           
   return xsec;
}

// **************************************************************************
   
template <class PhysModelType, bool isConversion>   
int
GammaModelWrapper<PhysModelType,isConversion>::
   SampleSecondaries(LightTrack &             projectile,
                     std::vector<LightTrack> & /*secondaryTracks*/,
                     Geant::GeantTaskData    *td )
{
  using std::cout;
  using std::endl;
  using  ThreeVector= vecgeom::Vector3D<double>;

  const  int electronCode = Electron::Definition()->GetInternalCode();
  const  int positronCode = Positron::Definition()->GetInternalCode();  
  const  int    gammaCode =    Gamma::Definition()->GetInternalCode();

  static bool firstCall= true;
  if( 0 ) { // firstCall ) {
    cout << " GammaModelWrapper - method SampleSecondaries() called for model " << this->GetName() << endl;     
    cout << " Gamma    'internal' code/type = " << electronCode << endl;
    cout << " Electron 'internal' code/type = " << gammaCode << endl;  
    cout << " Positron 'internal' code/type = " << positronCode << endl;
    firstCall= false;
  }
  assert( projectile.GetGVcode() == gammaCode );

  cout << " Gamma Wrapper - sampling model " << this->GetName() << endl;
  
  GUTrack  guSecondary;
  GUTrack  guProjectile= ConvertLightToGUTrack(projectile);
  
  int         matCutId= projectile.GetMaterialCutCoupleIndex();
  ExtraInfo *extraInfo= projectile.GetExtraInfo();
  int   ZtargetElement= projectile.GetTargetZ();
  double        weight= projectile.GetWeight();
  double          time= projectile.GetTime();

  double      initialEn =   guProjectile.E;
  ThreeVector initialMom= { guProjectile.px, guProjectile.py, guProjectile.pz };
  ThreeVector initialDir;   // double      initialMomMag=

  VectorToUnitAndMag( initialMom, initialDir );  
  
  int       idProjectile= projectile.GetTrackIndex();

  cout << " ****** 'Input' of Interaction :  ***************** " << endl;  
  cout << " **     Input projectile: " << guProjectile << endl;
  
  fVecPhysModel->template Interact<vecCore::backend::Scalar>( guProjectile, ZtargetElement, guSecondary);
  // **********  ******** ********

  cout << " ****** Results of Interact(ion) : **************** " << endl;
  cout << " ** 1. Changed projectile: " << guProjectile << endl;
  // cout << " ************************************************** " << endl;  
  cout << " ** 2.  Secondary :         " << guSecondary  << endl;
  // cout << " ************************************************** " << endl;

  ThreeVector finalMom=    { guProjectile.px, guProjectile.py, guProjectile.pz };
  ThreeVector secondaryMom= { guSecondary.px, guSecondary.py, guSecondary.pz };

#if 1 // CHECK_EACH_EN_MOMENTUM
  bool verbose= true;
  if( fGammaSurvives ){
     CheckEPGamma( guProjectile, "Surviving gamma", verbose );
     // CheckEPGamma( finalMom, guProjectile.E, "Surviving gamma", verbose );
  } else {
     // guProjectile.particleType= electronCode;
     // CheckEPMass( guProjectile, "Primary (Changed to e)" );
     CheckEnergyMomentum( finalMom,  guProjectile.E, electronCode, "Primary (Changed to e)", verbose );
  }
  CheckEnergyMomentum( secondaryMom, guSecondary.E, electronCode, "Secondary electron from gamma", verbose );
#else  
  // Check the relationship between E and p for output particles
  const double electronMass= 1000.0 * geant::kElectronMassC2; // geant::kElectronMassC2;
  if( fGammaSurvives ){
     double pEratio = finalMom.Mag() / guProjectile.E;
     if( fabs(pEratio - 1.0) > 1.0e-6 ) 
        cout << " p/E ratio (gamma-out) = " << pEratio << " r-1= "  << pEratio - 1.0 << endl;
     assert( fabs( finalMom.Mag() - guProjectile.E ) < 1.0e-4 * guProjectile.E );     
  } else {
     assert( GoodMass( guProjectile.E, finalMom, electronMass ) );
  }
  // Secondary is always an electron / positron
  // assert( GoodMass( guSecondary.E, secondaryMom, electronMass ) );
  if( ! GoodMass( guSecondary.E, secondaryMom, electronMass ) ){
     double estimM= EstimatedMass( guSecondary.E, secondaryMom );
     cout << " BAD mass from E/p - for secondary : "
          << " Mass estimated = " << estimM
          << " vs expected = " << electronMass
          << "  rel. diff = " << (estimM/electronMass) - 1.0
          << endl;
  }
#endif
#if  1   // DEBUG_CONSERVATION_LAWS
  //     DEBUG printout - START
  double initialMomentumMag= initialMom.Mag();
  {
     ThreeVector totalMomentumOut= finalMom + secondaryMom;  
     ThreeVector diffMomentum= initialMom - totalMomentumOut;
     double      diffEnergy  = initialEn - guProjectile.E - guSecondary.E;

     int oldPrc= cout.precision(10);
     cout << " particle 1 - out " << finalMom << endl;
     cout << " particle 2 - out " << secondaryMom << endl;
     cout << " Sum P out        " << totalMomentumOut << endl;
     cout << " P     in         " << initialMom   << endl;     
     cout.precision(oldPrc);
     
     if( initialMomentumMag == 0.0 ) { initialMomentumMag= -0.1e-70; }
     cout << " DiffE = ";
     if (diffEnergy == 0.0) {cout << "0 ";} else { cout << diffEnergy; }
     cout << " (UnRotated) d-Momentum= " << diffMomentum
          << " mag(dP)= " << diffMomentum.Mag() << " relative: " <<  diffMomentum.Mag() / initialMomentumMag
          << " d|P|= " << initialMomentumMag - totalMomentumOut.Mag()
          << " relative diff " << totalMomentumOut.Mag() / initialMomentumMag - 1.0
          << endl;
  }   // DEBUG printout - END
#endif

  // Rotate outgoing momentum directions - and re-check conservation
  ThreeVector finalMomDir;
  double finalMomMag = VectorToUnitAndMag( finalMom, finalMomDir );
  
  ThreeVector secondaryMomDir;
  double secondaryMomMag = VectorToUnitAndMag( secondaryMom, secondaryMomDir );  

  // rotate back to lab frame - 1st argument is changed in place
  // RotateToLabFrame(     finalMomDir, initialDir );
  // RotateToLabFrame( secondaryMomDir, initialDir );

  finalMom=     finalMomMag     * finalMomDir;
  secondaryMom= secondaryMomMag * secondaryMomDir;  

#if  1 // DEBUG_CONSERVATION_LAWS  
  {
     ThreeVector totalMomentumOutRot= finalMom + secondaryMom;     
     ThreeVector diffMomentumRot= initialMom - totalMomentumOutRot;
     double      diffEnergyRot  = initialEn - guProjectile.E - guSecondary.E;

     cout << " Sum P out (Rot)  " << totalMomentumOutRot << endl;
     
     cout << " DiffE = ";
     if (diffEnergyRot == 0.0) {cout << "0 ";} else { cout << diffEnergyRot; }
     cout << " (Rotated)   d-Momentum= " << diffMomentumRot
          << " mag(d-P)= " << diffMomentumRot.Mag() << " relative: " <<  diffMomentumRot.Mag() /  initialMomentumMag
          << " d|P|= " << initialMomentumMag - totalMomentumOutRot.Mag()
          << " relative diff " << totalMomentumOutRot.Mag() / initialMomentumMag - 1.0
          << endl;
  }
  
  if( projectile.GetGVcode() != gammaCode ){
     cout << " - projectile changed from gamma (code= " << gammaCode << " ) to code "
          << projectile.GetGVcode() << endl;
  }
#endif
  
  // To DO:  if the models do not respect the production thresholds (they don't),
  //         we can (should?) apply the production threshold for ongoing gammas & electrons.
  //         Ongoing positrons should be cut only if the production threshold is above 2 * massElectron
  
  // status = guProjectile.status; // GetStatus();
  // if( status == GUTrack::kAlive ){ ...
  int numberSecondaries= isConversion ? 2 : 1;  // Conversion creates two new tracks - others 1.

  //  It seems that the list "secondaryTracks" passed in argument is ignored ... e.g. in SeltzerBergerBremmModel
  //   let's copy what that does instead ...
  //
  PhysicsData* physicsData= td->fPhysicsData;
  physicsData->PrepareForSecondaries( numberSecondaries );

  std::vector<LightTrack>& sectracks = physicsData->GetListOfSecondaries();  
  int secIndx   = physicsData->GetNumUsedSecondaries();

  guSecondary.parentId= idProjectile;
  
  if( isConversion ) {
     // Deal with the HACK in Conversion - where the projectile photon is turned into a lepton!!

     // Lets assume ongoing track is an electron
     guProjectile.particleType= electronCode;
     guSecondary.particleType= positronCode;
     guProjectile.parentId= idProjectile;

     LightTrack convertedTrack = ConvertGUtoLightTrack( guProjectile,
                                                      LTrackStatus::kAlive,
                                                      // guProjectile.particleType, // PID  right now
                                                      matCutId,
                                                      time,
                                                      weight,
                                                      0.0,    // eDeposit = 0 for discrete process
                                                      extraInfo
        );
     // secondaryTracks.push_back( convertedTrack );  //  A copy is made into the vector !?
     sectracks[secIndx] = convertedTrack;
     secIndx++;

     // The original photon is killed
     projectile.SetTrackStatus( LTrackStatus::kKill /*ed*/ );  // Defined in LightTrack.h
  }
  if( fGammaSurvives ){
     CheckEPGamma( guProjectile, "Surviving gamma", true );
     
     projectile.SetKinE( guProjectile.E );
     // projectile.SetDirection( guProjectile.GetDirection() );
     // double px=guProjectile.px, py=guProjectile.py, pz=guProjectile.pz;
     // double invMomentum = 1.0 / std::sqrt( px*px + py*py + pz*pz ); // ( pMag2 );
     // projectile.SetDirection( px * invMomentum, py * invMomentum, pz * invMomentum);
     projectile.SetDirection( finalMom.x(), finalMom.y(), finalMom.z() );
     // projectile.SetDirection( finalMom );
  }

  guSecondary.particleType= fEmittedType;
  LightTrack  lightSecondary= ConvertGUtoLightTrack( guSecondary,
                                                     LTrackStatus::kAlive,
                                                     // fEmittedType, // electronCode, // guSecondary.particleType,
                                                     matCutId,
                                                     time,
                                                     1.0,
                                                     0.0,    // eDeposit = 0 for discrete process
                                                     extraInfo );
  // secondaryTracks.push_back(lightSecondary);

  lightSecondary.SetTrackIndex( idProjectile );  // Expects parent index for new tracks !!

  cout << endl;
  cout << ">>> Secondary (light) is " << lightSecondary << endl;
  
  sectracks[secIndx] = lightSecondary;
  secIndx++;
  
  /* sectracks[secIndx].SetDirX(gamDirX); */
  /* sectracks[secIndx].SetDirY(gamDirY); */
  /* sectracks[secIndx].SetDirZ(gamDirZ); */
  /* sectracks[secIndx].SetKinE(gammaEnergy); */
  /* sectracks[secIndx].SetGVcode(fSecondaryInternalCode);  // gamma GV code */
  /* sectracks[secIndx].SetMass(0.0); */
  /* sectracks[secIndx].SetTrackIndex(track.GetTrackIndex()); */
  
  // void Interact(GUTrack &projectile, const int targetElement, GUTrack &secondary);

  cout << " GammaModelWrapper - end of method SampleSecondaries() called for model " << this->GetName() << endl;
  
  return numberSecondaries;
}

/**
 * @brief Method to compute macroscopic cross section for a given MaterialCuts, Particle, kinetic energy.
 *
 * Mandatory method for discrete processes
 * Called at initialization to build the lambda tables by the PhysicsManagerPerParticle
 * (via its registered PhysicsProcess object. )
 *
 * @param[in] matcut      Pointer to the MaterialCuts
 * @param[in] kinenergy   Kinetic energy of the Particle
 * @param[in] particle    Pointer to the Particle object
 * @return    Macroscopic cross section computed in internal [1/length] units 
 *            for the given combination of MaterialCuts/Material, particle type
 *            and particle kinetic energy.
 */
template <class PhysModelType, bool isConversion>   
double
GammaModelWrapper<PhysModelType,isConversion>::  
   ComputeMacroscopicXSection(const  MaterialCuts *matcut,
                              double kinEnergy,
                              const  Particle * /*particle*/)
{
   double xsec = 0.0;
   assert (matcut != 0);   
   const Material *mat =  matcut->GetMaterial();
   auto  elementVec=     mat->GetElementVector();
   auto  relNumAtomsVec= mat->GetRelativeNumberOfAtomsPerVolume();   
   using  std::cout;
   using  std::endl;

   // double kinE_inGeV= kinEnergy;
   kinEnergy *= 1000.0; // ( CLHEP::GeV / geant::GeV ) ; // Since, for now, VecPhys is using CLHEP units
   
   static bool firstCall= true;
   if( firstCall ) { 
     std::cout << " GammaModelWrapper - method ComputeMacroscopicXSection() called for model " << this->GetName() << std::endl;
     firstCall= false;
   }
   
   assert (mat != 0);
   // assert (particle == Gamma::Gamma() );

   // double xsec2 = 0.0;  // For cross-checking

   static const MaterialCuts* oldMatcut = nullptr; // For initial debugging only.
   bool debugPrint= 1; // (matcut != oldMatcut); 
   oldMatcut = matcut;
   
   if (kinEnergy>GetLowEnergyUsageLimit() && kinEnergy < GetHighEnergyUsageLimit())
   {
      double densityXnumAv= geant::kAvogadro * mat->GetDensity();

      //  Density of material is in internal [weight/ length/f$^3/f$] units.
      
      int numElements= mat->GetNumberOfElements();

      if( debugPrint )      
         cout << "Call for X-section of material " << mat->GetName()
              << "  #components= " << numElements << endl;

      for( int i=0; i < numElements; i++)
      {
         Element  *elem   =   elementVec[i];
         double   Zelement=   elem->GetZ();
         double  xsecPart= 0.0;

         //  X-section right from VecPhys - uses CLHEP units
         double xsecPerAtom= fVecPhysModel // ->GetG4CrossSection( Zelement, kinEnergy );
                             -> template CrossSectionKernel<vecCore::backend::Scalar> ( kinEnergy, Zelement );
                         // -> ComputeXSectionPerAtom( Zelement, matcut, kinEnergy );
#if 1
         xsecPerAtom /= CLHEP::barn;

         std::cout << " GMWrapper CXSPA call - model " << this->GetName()
                   << " Element= " << elem->GetName() << " Z= " << int(Zelement)
                   << " Ekin= " << std::scientific << kinEnergy / CLHEP::MeV << " MeV "
                   << " X-sec/atom = " << std::scientific
                   << xsecPerAtom << " barn "  << std::endl;

         // Convert to GeantV Units
         xsecPerAtom *= geant::barn;
#else
         // Unit correction - since, for now, VecPhys is using CLHEP units         
         const double barnCorrection= geant::barn / CLHEP::barn;
         xsecPerAtom *= barnCorrection;
#endif        
         //  Density is from Real-Physics - uses new GeantV units
         double numAtomsPerVol= relNumAtomsVec[i] * densityXnumAv / elem->GetA();

         // Both are now in GeantV units
         xsecPart = xsecPerAtom * numAtomsPerVol;
            
         xsec += xsecPart;
         if( debugPrint )
           cout << "   Z: " << Zelement << "  #atoms= " << numAtomsPerVol << "  xSecPart= " << xsecPart << endl;
      }
   }
 
   // constexpr double unitDensityGrPerCm3_geant = geant::gram / (geant::cm * geant::cm * geant::cm );
   // constexpr double unitDensityGrPerCm_clhep = CLHEP::gram / CLHEP::cm;   
   
   std::cout << " GMWrapper CXSPA call - model " << this->GetName()
             << " Material= " << mat->GetName() << " Ekin= " << std::scientific << kinEnergy
             << " X-sec = " << std::scientific
             << xsec        << " g / cm "
             << std::endl;
   std::cout << std::fixed;

   // xsec *= unitDensityGrPerCm_geant;
   
   return xsec;   
}

// Alternative is to use the VecPhys method:
//  VECCORE_ATT_HOST double EmModelBase<EmModel>::G4CrossSectionPerVolume(const vecgeom::Material *material, double energy)

} // namespace geantphysics

#endif // GammaModelWrapper_H
