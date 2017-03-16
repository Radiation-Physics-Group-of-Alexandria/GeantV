//
//  J. Apostolakis & S.Y. Jun,  March 2017
//
#include "ComptonProcess.h"
#include "ComptonKleinNishina.h"
#include "GammaModelWrapper.h"

namespace geantphysics {

ComptonProcess::ComptonProcess(const std::string &name)
: EMPhysicsProcess(name)
{
  // set process type to be an EM process.  Note: i.e. not relevant for E-loss
  SetType(ProcessType::kElectromagnetic);
  // Set it as a discrete process
  SetIsContinuous(false);
  SetIsDiscrete(true);
  // fill the list of particles that this process can be used to (note: models need to set either to be for e- or e+)
  AddToListParticlesAlloedToAssigned(Gamma::Definition());

  double lowE_inKeV = 1.0;
  double highE_inGeV = 10.0; 
  
  // ComptonKleinNishina(Random_t *states = 0, int threadId = -1);    
  auto comptonKNvec= new vecphys::ComptonKleinNishina();

  comptonKNvec->SetLowEnergyLimit ( lowE_inKeV  * CLHEP::keV);  // geant::keV);
  comptonKNvec->SetHighEnergyLimit( highE_inGeV * CLHEP::GeV);  // geant::GeV);
  
  EMModel *comptonModel=
     new GammaModelWrapper<vecphys::ComptonKleinNishina>(
        "Compton Klein Nishina model (Wrapped VecPhys model)",
        comptonKNvec);
  comptonModel->SetLowEnergyUsageLimit ( lowE_inKeV  * geant::keV);
  comptonModel->SetHighEnergyUsageLimit( highE_inGeV * geant::GeV);
  std::cout << "Registered Compton model for energy range: "
               << lowE_inKeV << " keV ( = "<< lowE_inKeV  * geant::keV << ")  to "
            << highE_inGeV << " GeV ( = " << highE_inGeV * geant::GeV << " ) " << std::endl;
  AddModel( comptonModel );
}

} // namespace geantphysics
