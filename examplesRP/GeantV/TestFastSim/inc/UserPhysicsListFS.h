
#ifndef FS_USERPHYSICSLIST_H
#define FS_USERPHYSICSLIST_H

#include "PhysicsList.h"
// for the MSCSteppingAlgorithm enums
#include "MSCModel.h"

#include <string>


namespace userfsapplication {

/**
 * @brief User physics list for TeatFastSim.
 *
 * The physics list contains the available GeantV standard EM interactions and a user defined StepMaxProcess for e-/e+.
 * The step limit value of the this StepMaxProcess and the multiple Coulomb scattering process stepping algorithm type
 * are configurable from input arguments.
 *
 * @class   UserPhysicsListFS
 * @author  M Novak
 * @date    July 2017
 */

class UserPhysicsListFS : public geantphysics::PhysicsList {
public:
  // CTR
  UserPhysicsListFS(const std::string &name);
  // DTR
 ~UserPhysicsListFS();
  // interface method to assigne physics-process to particles
  virtual void Initialize();

  // public method to allow multiple scattering step limit configuration
  void SetMSCStepLimit(geantphysics::MSCSteppingAlgorithm stepping);
  // public method to allow configuration of the user-defined step-max process
  void SetStepMaxValue(double val);

private:
  geantphysics::MSCSteppingAlgorithm  fMSCSteppingAlgorithm;
  double                              fStepMaxValue;
};

}      //  namespace userfsapplication


#endif // USERPHYSICSLIST_H
