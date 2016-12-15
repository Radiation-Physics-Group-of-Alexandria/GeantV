#ifndef GAMMAPHOTOELECTRICPROCESS_H
#define GAMMAPHOTOELECTRICPROCESS_H

#include "EMPhysicsProcess.h"

#include "Electron.h"
#include "Positron.h"
#include "Gamma.h"

#include <string>

namespace geantphysics {

class GammaPhotoElectricProcess : public EMPhysicsProcess {
public:
  // CTR
  GammaPhotoElectricProcess(const std::string &name = "gPhotoElectric");
};

} // namespace geantphysics

#endif
