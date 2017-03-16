
#ifndef ConversionProcess_H
#define ConversionProcess_H

#include "EMPhysicsProcess.h"

// #include "Electron.h"
// #include "Positron.h"
#include "Gamma.h"

#include <string>

namespace geantphysics {

class ConversionProcess : public EMPhysicsProcess {
public:
  // CTR
  ConversionProcess(const std::string &name = "Conversion");
};

} // namespace geantphysics

#endif
