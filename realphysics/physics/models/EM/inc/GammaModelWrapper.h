
#ifndef GammaModelWrapper_H
#define GammaModelWrapper_H

#include <string>
#include <vector>
#include <cassert>

#include "base/VecPhys.h"
// #include "base/SystemOfUnits.h"   //  The OLD units - clhep namespace

// #include "VecCore/Backend/Scalar.h"
// #include "GeantTaskData.h"

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
class GeantTaskData;   
class Material;
class MaterialCuts;
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
   * @brief CTR
   */
  GammaModelWrapper(const std::string&  name,
                    PhysModelType*      vecPhysModel )
      : EMModel( name )
                    // bool             isConversion = false ) // double minimumPrimaryEnergy= 0.0);
  {
     fVecPhysModel= vecPhysModel;
  }
  /**
   * @brief DTR
   */
  ~GammaModelWrapper() {}

//
// The following 5 virtual methods might be implemented by the derived electromagnetic models.
//
  // will set the physics parameter object that belongs to the regions(s) where this model is active
  /**
   * @brief The base class implementation Initialize() must be called in the first line because to set
   *        the PhysicsParameters member variable for the regions(s) where this model is active.
   *        This method is called when its EMPhysicsProcess, is initialized by the
   *        EmModelManager member of the EMPhysicsProcess.
   */
  void  Initialize() final;

  /**
   * @brief Method to compute stopping power relevant to energy loss models only.
   *
   */
  virtual double ComputeDEDX(const MaterialCuts * /*matcut*/, double /*kinenergy*/, const Particle * /*particle*/,
                             bool /*istotal=false*/ ) { return 1.0; }

  /**
   * @brief Method to compute macroscopic cross section for a given MaterialCuts, Partcile, kinetic energy.
   *
   * ALL DISCRETE MODELS (i.e. those that describes interaction happening at the post-step point) NEEDS TO IMPLEMENT
   * this method. This method is called at initialization to build the lambda tables through the corresponding
   * PhysicsProcess by the PhysicsManagerPerParticle object.
   *
   * @param[in] matcut      Pointer to the MaterialCuts object in which the macroscopic cross section must be computed.
   * @param[in] kinenergy   Kinetic energy of the Partcile at which the macroscopic cross section must be computed.
   * @param[in] particle    Pointer to the Partcile object for which the macroscopic cross section must be computed.
   * @return    Macroscopic cross section computed by the given electromagnetic model in internal [1/length] units for
   *            the given Particle, MaterialCuts/Material and particle kinetic energy combination.
   */
  virtual double ComputeMacroscopicXSection(const MaterialCuts * /*matcut*/, double /*kinenergy*/,
                                            const Particle * /*particle*/);

  /**
   * @brief Method to compute atomic cross section for a given Element, MaterialCuts, Partcile, kinetic energy.
   *
   * Method is called at initialization to build the lambda tables through the corresponding PhysicsProcess by
   *  - by the PhysicsManagerPerParticle object and 
   *  - by the EMElementSelector object if it was requested in the derived class Initialize() method.
   * 
   * Note, that the MaterialCut object, that corresponds to the material-region to which the current element belongs to
   * is also provided to include any material, or cut dependences. However, this infomation does not necessary to use by
   * each GammaModelWrapper-s (e.g. those models that describes discrete interaction that doesn't depend on production
   * threshold and material dependent properties won't use this extra information).
   *
   * @param[in] elem        Pointer to the Element object for which the atomic cross section must be computed.
   * @param[in] matcut      Pointer to the MaterialCuts object in which the atomic cross section must be computed.
   * @param[in] kinenergy   Kinetic energy of the Partcile at which the atomic cross section must be computed.
   * @param[in] particle    Pointer to the Partcile object for which the atomic cross section must be computed.
   * @return    Atomic cross section computed by the given electromagnetic model in internal [length^2] units for
   *            the given ELement, Particle, MaterialCuts/Material and particle kinetic energy combination.
   */
  double ComputeXSectionPerAtom(const Element * /*elem*/, const MaterialCuts * /*matcut*/, double /*kinenergy*/,
                                const Particle * /*particle*/) final;

  int    SampleSecondaries(LightTrack & /*track*/, std::vector<LightTrack> & /*sectracks*/,
                                   Geant::GeantTaskData * /*td*/) final;

  /**
   * @brief Method to obtain minimum primary particle kinetic energy at which the discrete part (if any) of the 
   *        interaction can happen i.e. kinematic minimum for the discrete part of the interaction.
   *
   * ALL DISCRETE MODELS (i.e. those that describes interaction happening at the post-step point) THAT HAS PRODUCTION
   * CUT DEPENDENCE NEEDS TO IMPLEMENT this method. It is used e.g. to build the energy grid of the target element
   * selector (if it was requested by the derived model) for providing run-time functionality to select target element
   * for the discrete interaction.
   *
   * @param[in] matcut      Pointer to the MaterialCuts object in which the discrete kinematic minimum is requested.
   * @param[in] particle    Pointer to the Partcile object for which the discrete kinematic minimum is requested.
   * @return    Minimum primary particle kinetic energy at which the discrete interaction can happen in interanl
   *            [energy] units for the given Particle, MaterialCuts/Material combination.
   */
  double MinimumPrimaryEnergy(const MaterialCuts * /*matcut*/, const Particle * /*part*/) const final
  { return  ( isConversion ? 2.0 * geant::kElectronMassC2 : 0.0 ); }

protected:
  // Some short cuts will be used in methods below, for speed - knowing 
  //    the type of primary and the implementation of the interactions
  GUTrack    ConvertLightToGUTrack(LightTrack &lightTrack) const;
  LightTrack ConvertGUtoLightTrack(GUTrack      &guTrack,
                                   LTrackStatus  status,  
                                   int           materialCutCoupleIndex,
                                   double        time= -1.0,
                                   double        weight = 1.0,                                   
                                   double        eDeposit = 0.0,
                                   ExtraInfo*    aExtraInfo = 0
                                   ) const ;

  // friend PhysModelType;
  
private:
  PhysModelType*  fVecPhysModel;  
  double          fModelminimumPrimaryEnergy= 0.0;

};


// template <class PhysModelType, bool isConversion>         
// void
//    GammaModelWrapper<PhysModelType,isConversion>:: METHOD
// { }

template <class PhysModelType, bool isConversion>         
void
   GammaModelWrapper<PhysModelType,isConversion>::Initialize()
{
  //  fVecPhysModel->template Initialize<vecCore::backend::Scalar>Initialization();
  fVecPhysModel->Initialization();
}

template <class PhysModelType, bool isConversion>         
GUTrack
   GammaModelWrapper<PhysModelType,isConversion>::ConvertLightToGUTrack(LightTrack &lightTrack) const
{
  GUTrack guTrack;

  guTrack.status= 1;                             // lightTrack.GetTrackStatus();
  guTrack.particleType= lightTrack.GetGVcode();
  guTrack.id= lightTrack.GetTrackIndex();       // counter
  guTrack.parentId= 0;                           // id of parent
  guTrack.proc= lightTrack.GetProcessIndex();    // index of physics process
  guTrack.x= 0.0;
  guTrack.y= 0.0;
  guTrack.z= 0.0;
  // double restMass= guTrack.GetMass();
  assert(guTrack.particleType == 22); // Gamma.ParticleType);
  double eKin= lightTrack.GetKinE();
  // double momentum= std::sqrt(eKin*(eKin+2*restMass));
  double momentum= eKin;  // Shortcut for gamma :   E = E_kin = p c   
  guTrack.px= momentum * lightTrack.GetDirX();
  guTrack.py= momentum * lightTrack.GetDirY();
  guTrack.pz= momentum * lightTrack.GetDirZ();
  guTrack.E= eKin;
  guTrack.q= 0;  // Shortcut for photon 
  guTrack.nint= lightTrack.GetNumOfInteractionLegthLeft();
  //  guTrack.lambda= lightTrack.GetIntLen();   // interaction length
  guTrack.s= lightTrack.GetStepLength();    // step length ??

  return guTrack;
}

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
  double aMass = (guTrack.id == 22) ? 0.0: geant::kElectronMassC2; 
  return LightTrack( 
             status,  // LTrackStatus::kAlive
             guTrack.particleType,  // aGVcode,
             guTrack.id,    //  aGTrackIndex,   => should be -1/null for new particles
             materialCutCoupleIndex,
             0, // aProcessIndex ==> this particle is new, no process is chosen for it
             0, // aTargetZ,
             0, // aTargetN, 
             guTrack.px, guTrack.py, guTrack.pz,
             guTrack.E, 
             aMass,
             time, 
             weight,
             guTrack.s,
             eDeposit,
             guTrack.nint,    // aNintLen
             guTrack.lambda,  // aIntLen,
             aExtraInfo
         );
}

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

// template <class PhysModelType, PhysModelType VecPhysModel, bool isConversion>            
template <class PhysModelType, bool isConversion>   
int
GammaModelWrapper<PhysModelType,isConversion>::
   SampleSecondaries(LightTrack & projectile, std::vector<LightTrack> &secondaryTracks,
                     Geant::GeantTaskData * /*td*/ )
{
  GUTrack  guSecondary;
  GUTrack  guProjectile= ConvertLightToGUTrack(projectile);
  int      matCutId=    projectile.GetMaterialCutCoupleIndex();
  ExtraInfo *extraInfo= projectile.GetExtraInfo();
  int  ZtargetElement= projectile.GetTargetZ();
  double    weight=    projectile.GetWeight();
  double    time=      projectile.GetTime();

  fVecPhysModel->template Interact<vecCore::backend::Scalar>( guProjectile, ZtargetElement, guSecondary);
  // *****************************

  // To DO:  if the models do not respect the production thresholds (they don't),
  //         we can (should?) apply the production threshold for ongoing gammas & electrons.
  //         Ongoing positrons should be cut only if the production threshold is above 2 * massElectron
  
  // status = guProjectile.status; // GetStatus();
  // if( status == GUTrack::kAlive ){ ...

  if( isConversion ) {
     // Deal with the HACK in Conversion - where the projectile photon is turned into a lepton!!
     LightTrack ongoingTrack = ConvertGUtoLightTrack( guProjectile,
                                                      LTrackStatus::kAlive,
                                                      // guProjectile.particleType, // PID  right now
                                                      matCutId,
                                                      time,
                                                      weight,
                                                      0.0,    // eDeposit = 0 for discrete process
                                                      extraInfo
        );
     secondaryTracks.push_back( ongoingTrack );  //  A copy is made into the vector !? 
     projectile.SetTrackStatus( LTrackStatus::kKill /*ed*/ );  // Defined in LightTrack.h
  } else {      
     projectile.SetKinE( guProjectile.E );
     // projectile.SetDirection( guProjectile.GetDirection() );
     double px=guProjectile.px, py=guProjectile.py, pz=guProjectile.pz;
     double invMomentum = 1.0 / std::sqrt( px*px + py*py + pz*pz ); // ( pMag2 );
     projectile.SetDirection( px * invMomentum, py * invMomentum, pz * invMomentum);
  }

  LightTrack  lightSecondary= ConvertGUtoLightTrack( guSecondary,
                                                     LTrackStatus::kAlive,
                                                     guSecondary.particleType,
                                                     matCutId,
                                                     time,
                                                     0.0,    // eDeposit = 0 for discrete process
                                                     extraInfo
     );
  secondaryTracks.push_back(lightSecondary);

  // void Interact(GUTrack &projectile, const int targetElement, GUTrack &secondary);
  return 2;
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

   assert (mat != 0);
   // assert (particle == Gamma::Gamma() );

   double xsec2 = 0.0;  // For cross-checking

   static const MaterialCuts* oldMatcut = nullptr; // For initial debugging only.
   bool debugPrint= (matcut != oldMatcut); 
   oldMatcut = matcut;
   
   if (kinEnergy>GetLowEnergyUsageLimit() && kinEnergy < GetHighEnergyUsageLimit()) {

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
         double numAtomsPerVol= relNumAtomsVec[i] * densityXnumAv / elem->GetA();

         xsecPart= numAtomsPerVol * fVecPhysModel // ->GetG4CrossSection( Zelement, kinEnergy );
                           -> template CrossSectionKernel<vecCore::backend::Scalar> ( kinEnergy, Zelement );
                         // -> ComputeXSectionPerAtom( Zelement, matcut, kinEnergy );
         xsec += xsecPart;
         if( debugPrint )
           cout << "   Z: " << Zelement << "  #atoms= " << numAtomsPerVol << "  xSecPart= " << xsecPart << endl;
         
         xsec2 += numAtomsPerVol *
                  fVecPhysModel ->
                        template CrossSectionKernel<vecCore::backend::Scalar> (kinEnergy, Zelement);
      }
   }
   assert( std::fabs(xsec - xsec2) < 0.5e5 * (xsec+xsec2) );
   
   xsec *= 100.0; // ( clhep::barn / geant::barn ) ; // Since, for now, VecPhys is using CLHEP units
   return xsec;   
}
  
} // namespace geantphysics

#endif // GammaModelWrapper_H
