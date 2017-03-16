//
//  J. Apostolakis & S.Y. Jun,  March 2017
//
#include "ConversionProcess.h"
#include "ConversionBetheHeitler.h"
#include "GammaModelWrapper.h"

namespace geantphysics {

ConversionProcess::ConversionProcess(const std::string &name)
: EMPhysicsProcess(name) {
   using vecphys::VECPHYS_IMPL_NAMESPACE::ConversionBetheHeitler;

  // set process type to be an EM process.  Note: i.e. not relevant for E-loss
  SetType(ProcessType::kElectromagnetic);
  // Set it as a discrete process
  SetIsContinuous(false);
  SetIsDiscrete(true);
  // fill the list of particles that this process can be used to (note: models need to set either to be for e- or e+)
  AddToListParticlesAlloedToAssigned(Gamma::Definition());

  // ConversionModel(Random_t *states = 0, int threadId = -1);  
  auto vecphysConversionSG= new ConversionBetheHeitler(nullptr);
  double lowE_inMeV =    1.0;
  double highE_inGeV = 100.0;
  
  vecphysConversionSG->SetLowEnergyLimit (  lowE_inMeV * CLHEP::MeV);  // geant::MeV);
  vecphysConversionSG->SetHighEnergyLimit( highE_inGeV * CLHEP::GeV);  // geant::GeV);
  
  EMModel *conversionModel=
     new GammaModelWrapper<ConversionBetheHeitler, true>(
        "Photon Conversion process model (Wrapped VecPhys model)",
        vecphysConversionSG);
  conversionModel->SetLowEnergyUsageLimit (  lowE_inMeV * geant::MeV);
  conversionModel->SetHighEnergyUsageLimit( highE_inGeV * geant::GeV);

  std::cout << "Registered Conversion model for gammas in energy range: "
            << lowE_inMeV << " MeV ( = "<< lowE_inMeV  * geant::MeV << ")  to "
            << highE_inGeV << " GeV ( = " << highE_inGeV * geant::GeV << " ) " << std::endl;  
  
  AddModel( conversionModel );
}

} // namespace geantphysics
