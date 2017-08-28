#ifndef GEANTV_HEPMCGENERATORMULTFILES_H
#define GEANTV_HEPMCGENERATORMULTFILES_H

#include "HepMCGenerator.h"

namespace Geant {
inline namespace GEANT_IMPL_NAMESPACE {

class HepMCGeneratorMultFiles : public HepMCGenerator {
private:
  std::string filename;
  int currentOffset;
  int desiredOffset;

public:
  HepMCGeneratorMultFiles();
  void SetEventSource(const std::string &file, int offset);

  virtual GeantEventInfo NextEvent();
};

} // GEANT_IMPL_NAMESPACE
} // Geant

#endif // GEANTV_HEPMCGENERATORMULTFILES_H
