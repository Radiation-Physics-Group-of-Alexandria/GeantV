#ifndef GAMMACOMPTONPROCESS_H
#define GAMMACOMPTONPROCESS_H

#include "EMPhysicsProcess.h"

#include "Electron.h"
#include "Positron.h"
#include "Gamma.h"

#include <string>

namespace geantphysics {

class GammaComptonProcess : public EMPhysicsProcess {
public:
  GammaComptonProcess(const std::string &name = "gCompton");
};

} // namespace geantphysics

#endif
