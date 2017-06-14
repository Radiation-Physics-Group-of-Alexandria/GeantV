#ifndef HepMCGenerator_h
#define HepMCGenerator_h

#include "TPartIndex.h"
#include "PrimaryGenerator.h"

#include "HepMC/Reader.h"
#include "HepMC/GenEvent.h"
#include "HepMC/Search/FindParticles.h"

namespace Geant {
inline namespace GEANT_IMPL_NAMESPACE {

class HepMCGenerator : public PrimaryGenerator {
protected:
  HepMC::Reader *input_file;
  HepMC::FindParticles *search;
  void LoadFile(const std::string &filename);
public:
  HepMCGenerator();
  HepMCGenerator(std::string &filename);
  virtual ~HepMCGenerator();

  // set one GeantTrack primary track properties
  virtual void InitPrimaryGenerator();
  virtual GeantEventInfo NextEvent();
  virtual void GetTrack(int n, Geant::GeantTrack &gtrack);
  // used from Geant4 test-complex to take one primary track
  void GetTrack(int n, double &tpx, double &tpy, double &tpz, double &te, double &x0, double &y0, double &z0, int &pdg);

private:
  HepMCGenerator(const HepMCGenerator &);            // no imp.
  HepMCGenerator &operator=(const HepMCGenerator &); // no imp.

};

} // GEANT_IMPL_NAMESPACE
} // Geant

#endif
