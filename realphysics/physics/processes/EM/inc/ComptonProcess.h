
#ifndef ComptonProcess_H
#define ComptonProcess_H

#include "EMPhysicsProcess.h"

// #include "Electron.h"
// #include "Positron.h"
#include "Gamma.h"

#include <string>

namespace geantphysics {

class ComptonProcess : public EMPhysicsProcess {
public:
  // CTR
  ComptonProcess(const std::string &name = "Compton");
};

} // namespace geantphysics

#endif
