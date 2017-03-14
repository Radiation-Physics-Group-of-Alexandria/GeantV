
#ifndef PhotoElectricProcess_H
#define PhotoElectricProcess_H

#include "EMPhysicsProcess.h"

// #include "Electron.h"
// #include "Positron.h"
#include "Gamma.h"

#include <string>

namespace geantphysics {

class PhotoElectricProcess : public EMPhysicsProcess {
public:
  // CTR
  PhotoElectricProcess(const std::string &name = "PhotoElectricProcess");
};

} // namespace geantphysics

#endif
