//==================================================================
// @file NUDYCrossSection.h
// @author Abhijit Bhattacharyya
// @brief Interface of NUDY to GV for getting CS to use in HAD
//
// Like Hadronic CS we are initially showing elastic and inelastic
// Later it should cover thermal etc.
//======================================================================

#ifndef NUDY_CROSS_SECTION_H 
#define NUDY_CROSS_SECTION_H

#ifdef USE_ROOT
#include "NudyInterface.h"
#endif

// Forward declaration
namespace geantphysics {
  class LightTrack;
  class LightTrack_v;
  class HadronicCrossSection;
  class HadronicProcess;

  inline namespace GEANT_IMPL_NAMESPACE {
    class Isotope;
    class Material;
    class Element;
  }
}

namespace NudyPhysics {
  class NudyInterface;
}


namespace geantphysics {

  class NUDYCrossSection {
  public:
    // @brief Use NUDY Crossection if HadronicCrossSection method
    // HadronicCrossSection::IsApplilcable() returns false;
    // This could be done for elements having resonance and KE of "n" < 20 MeV.

    NUDYCrossSection();
    NUDYCrossSection(const std::string isoName, const int Nval, const int Zval,
		     const std::string interactionType, const double projectileKE
		     );

    virtual ~NUDYCrossSection();    // destructor

    //--------- Getters --------------
    std::string GetisoName() const { return fisoName; }           // returns name of the isotope
    int GetNval() const { return ftargetN; }                      // returns neucleon number of isotope
    int GetZval() const { return ftargetZ; }                      // returns atomic number of isotope
    double GetProjKE() const { return fprojEnergy; }              // returns KE of projectile

    virtual bool IsApplicable(const int projCode, const double projKE ); // this checks for neutron and KE < 20 MeV

    virtual double GetIsotopeCrossSection ( const int projectilecode, const double projectileKE,
					    const int Z, const int N
					    );
	

    //---------- Setters ----------------
    void SetisoName ( const std::string &isoname ) { fisoName = isoname; }
    void SetNval ( const int Nval ) { ftargetN = Nval; }
    void SetZval ( const int Zval ) { ftargetZ = Zval; }
    void SetProjKE ( const double projke ) { fprojEnergy = projke; }


  private:
    std::string fisoName;                     // Isotope Name for which CS required
    int ftargetN;                          // Nucleon number for the Isotope
    int ftargetZ;                          // Atomic number for the isotope 
    std::vector<int> fProjectileCodeVec;   // vector of GV code for allowed projectiles
    int fProjCode;                        // code for projectile as par code vector from HAD CS e.g. neutron = 3
    double fprojEnergy;                    // LE for the projectile
  };

}  // namespace

#endif

